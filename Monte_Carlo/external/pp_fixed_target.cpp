#include "Pythia8/Pythia.h"
#include <cmath>
#include <iomanip>
#include <fstream>
#include <limits>
#include <iostream>
#include <string>

int main() {
    const double T_beam = 1000.0; // Incident proton kinetic energy [GeV].
    const int N_events = 10000; // Independent collisions at the same beam energy.
    Pythia8::Pythia pythia;
    const double m_p = pythia.particleData.m0(2212); // Proton mass [GeV].
    const double E_beam = T_beam + m_p;

    // PDG ID 2212 is a proton. Frame 2 accepts laboratory total energies;
    // setting the target energy to its rest mass makes it stationary.
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:frameType = 2");
    pythia.settings.parm("Beams:eA", E_beam);
    pythia.settings.parm("Beams:eB", m_p);

    // Generate inelastic pp events, including diffractive processes.
    // Keep default particle decays, including pi0 -> gamma gamma.
    pythia.readString("SoftQCD:inelastic = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 12345");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    if (!pythia.init()) return 1;
    // Run from external/. Each row describes one independent collision.
    std::ofstream output("pp_fixed_target_events.csv");
    if (!output) return 1;
    output << "event_id,T_beam_GeV,nucleon1_pdg,T_nucleon1_GeV,nucleon2_pdg,T_nucleon2_GeV\n";
    output << std::setprecision(16);
    const double pz_initial = std::sqrt(T_beam * (T_beam + 2.0 * m_p));
    double max_dE = 0.0;
    double max_dp = 0.0;
    int N_no_proton = 0;
    int N_failed_attempts = 0;
    for (int event_id = 0; event_id < N_events; ++event_id) {
        bool generated = false;
        for (int attempt = 0; attempt < 10; ++attempt) {
            if (pythia.next()) {
                generated = true;
                break;
            }
            ++N_failed_attempts;
        }
        if (!generated) {
            std::cerr << "Generation failed at event " << event_id << '\n';
            return 1;
        }
        Pythia8::Vec4 p_final;
        int N_protons = 0;
        int N_nucleons = 0;
        int nucleon1_pdg = 0;
        int nucleon2_pdg = 0;
        double T_nucleon1 = std::numeric_limits<double>::quiet_NaN();
        double T_nucleon2 = std::numeric_limits<double>::quiet_NaN();
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Pythia8::Particle& particle = pythia.event[i];
            // Only particles remaining after the enabled default decays.
            if (!particle.isFinal()) continue;
            p_final += particle.p();
            if (particle.id() == 2212) ++N_protons;
            // Keep the two highest LAB kinetic energies among final protons
            // and neutrons, excluding antinucleons. Use each species' mass.
            // These are energies before neutron beta decay; no decay is imposed.
            if (particle.id() == 2212 || particle.id() == 2112) {
                const double T = particle.e() - pythia.particleData.m0(particle.id());
                if (N_nucleons == 0 || T > T_nucleon1) {
                    T_nucleon2 = T_nucleon1;
                    nucleon2_pdg = nucleon1_pdg;
                    T_nucleon1 = T;
                    nucleon1_pdg = particle.id();
                } else if (N_nucleons == 1 || T > T_nucleon2) {
                    T_nucleon2 = T;
                    nucleon2_pdg = particle.id();
                }
                ++N_nucleons;
            }
        }
        if (N_protons == 0) ++N_no_proton;
        const double dE = p_final.e() - (E_beam + m_p);
        const double dp = std::sqrt(p_final.px() * p_final.px() + p_final.py() * p_final.py()
                                   + (p_final.pz() - pz_initial) * (p_final.pz() - pz_initial));
        if (std::abs(dE) > max_dE) max_dE = std::abs(dE);
        if (dp > max_dp) max_dp = dp;
        output << event_id << ',' << T_beam << ',' << nucleon1_pdg << ',' << T_nucleon1 << ','
               << nucleon2_pdg << ',' << T_nucleon2 << '\n';
    }
    output.close();
    if (!output) return 1;
    std::cout << std::defaultfloat << std::setprecision(10) << "\nGenerated " << N_events << " independent inelastic pp events\n"
              << "Projectile kinetic energy = " << T_beam << " GeV\n"
              << "sqrt(s) = " << std::sqrt(4.0 * m_p * m_p + 2.0 * m_p * T_beam) << " GeV\n"
              << "No-final-proton fraction = " << double(N_no_proton) / N_events << '\n'
              << "Failed generation attempts = " << N_failed_attempts << '\n'
              << "Maximum |energy residual| = " << max_dE << " GeV\n"
              << "Maximum momentum residual = " << max_dp << " GeV\n"
              << "Wrote pp_fixed_target_events.csv\n";
    return 0;
}
