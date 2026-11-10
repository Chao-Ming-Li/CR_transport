// Independent pp collisions, followed by all available hadron decays.
#include "Pythia8/Pythia.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    // Run parameters: energies are total energies unless marked as kinetic.
    int N_events = 10000;
    if (argc > 1) {
        N_events = std::stoi(argv[1]);
    }
    if (N_events <= 0) {
        return 1;
    }
    const double E_min = 10.0; // Projectile TOTAL energy [GeV].
    const double E_max = 1e7; // 10 PeV [GeV].
    const int N_bins = 140;
    const double log_energy_range = std::log(E_max / E_min);

    // Initialize the pp collision and decay models.
    Pythia8::Pythia pythia;
    const double m_p = pythia.particleData.m0(2212);
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    // Generate and decay in the CM frame to avoid numerical cancellation
    // in highly boosted, small-Q neutron beta decays. Boost to LAB afterward.
    pythia.readString("Beams:frameType = 1");
    pythia.settings.parm("Beams:eCM", std::sqrt(2.0 * m_p * (E_max + m_p)));
    pythia.readString("Beams:allowVariableEnergy = on");
    pythia.readString("SoftQCD:inelastic = on");
    pythia.readString("LowEnergyQCD:inelastic = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 12345");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    // Infinite-time, vacuum decay approximation: enable every hadron decay
    // table, including n -> p e- anti-nu_e and weak hyperon decays.
    // PYTHIA samples decay kinematics; no neutron-to-proton energy rescaling.
    // Antihadrons decay to antiparticles and are not counted as protons.
    pythia.readString("ParticleDecays:limitTau0 = off");
    pythia.readString("ParticleDecays:limitTau = off");
    pythia.readString("ParticleDecays:limitRadius = off");
    pythia.readString("ParticleDecays:limitCylinder = off");
    if (!pythia.init()) {
        return 1;
    }
    // nextId visits the particle table in PDG-code order. The library returns
    // a shared pointer; auto keeps its library-defined type without copying
    // the entry. mayDecay changes the decay flag in PYTHIA's particle table.
    for (int id = 1; id != 0; id = pythia.particleData.nextId(id)) {
        auto entry = pythia.particleData.particleDataEntryPtr(id);
        if (entry->isHadron() && entry->sizeChannels() > 0) {
            pythia.particleData.mayDecay(id, true);
        }
    }
    pythia.particleData.mayDecay(2112, true);

    // Run from external/. Overwrite outputs; all energies are LAB energies.
    std::ofstream events("pp_powerlaw_events.csv");
    std::ofstream protons("pp_powerlaw_protons.csv");
    std::ofstream spectrum("pp_powerlaw_spectrum.csv");
    if (!events || !protons || !spectrum) {
        return 1;
    }
    events << "event_id,E_beam_GeV,process_code,N_protons\n";
    protons << "event_id,E_beam_GeV,particle_index,mother_pid,E_proton_GeV,T_proton_GeV\n";
    spectrum << "E_low_GeV,E_high_GeV,N_injected,N_protons,dN_injected_dE_per_event,dN_protons_dE_per_event\n";
    events << std::setprecision(16);
    protons << std::setprecision(16);
    spectrum << std::setprecision(16);

    // Logarithmic energy bins and injection random-number generator.
    // vector<int>(N_bins, 0) creates N_bins counters initialized to zero.
    std::vector<int> injected(N_bins, 0);
    std::vector<int> final_protons(N_bins, 0);
    // Include recoil/decay protons below the minimum injection energy.
    const double log_bin_min = std::log(m_p);
    const double dlogE = std::log((E_max + m_p) / m_p) / N_bins;
    std::mt19937_64 rng(67890);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    int N_protons_total = 0;
    int N_neutron_daughters = 0;
    double max_energy_residual = 0.0;
    double max_relative_residual = 0.0;
    int N_failed_attempts = 0;

    // Independent collisions: each sampled projectile undergoes one reaction.
    for (int i = 0; i < N_events; ++i) {
        // Uniform ln(E) means dN/dE proportional to E^-1; equal event weights.
        const double E_beam = E_min * std::exp(log_energy_range * uniform(rng));
        const double E_CM = std::sqrt(2.0 * m_p * (E_beam + m_p));
        if (!pythia.setKinematics(E_CM)) {
            return 1;
        }
        bool generated = false;
        for (int attempt = 0; attempt < 10; ++attempt) {
            if (pythia.next()) {
                generated = true;
                break;
            }
            ++N_failed_attempts;
        }
        if (!generated) {
            std::cerr << "Failed at event " << i << ", E_beam = " << E_beam << '\n';
            return 1;
        }
        // Finish decay chains before transforming all four-momenta to LAB.
        if (!pythia.moreDecays()) {
            return 1;
        }
        const double p_beam = std::sqrt((E_beam - m_p) * (E_beam + m_p));
        const double beta_CM = p_beam / (E_beam + m_p);
        pythia.event.bst(0.0, 0.0, beta_CM);
        // The logarithmic coordinate is nonnegative, so conversion to int
        // selects the containing bin by discarding its fractional part.
        const int bin_initial = int((std::log(E_beam) - log_bin_min) / dlogE);
        ++injected[bin_initial];
        int N_protons = 0;
        double E_final_sum = 0.0;
        for (int j = 0; j < pythia.event.size(); ++j) {
            const Pythia8::Particle& particle = pythia.event[j];
            if (!particle.isFinal()) {
                continue;
            }
            E_final_sum += particle.e();
            if (particle.id() == 2112) {
                std::cerr << "Undecayed final neutron at event " << i << '\n';
                return 1;
            }
            if (particle.id() != 2212) {
                continue;
            }
            ++N_protons;
            const int mother_pid = pythia.event[particle.mother1()].id();
            if (mother_pid == 2112) {
                ++N_neutron_daughters;
            }
            protons << i << ',' << E_beam << ',' << j << ',' << mother_pid << ','
                    << particle.e() << ',' << particle.e() - m_p << '\n';
            const int bin = int((std::log(particle.e()) - log_bin_min) / dlogE);
            if (bin < 0 || bin >= N_bins) {
                std::cerr << "Proton outside spectrum bins\n";
                return 1;
            }
            ++final_protons[bin];
        }
        const double residual = std::abs(E_final_sum - E_beam - m_p);
        const double relative_residual = residual / (E_beam + m_p);
        if (residual > max_energy_residual) {
            max_energy_residual = residual;
        }
        if (relative_residual > max_relative_residual) {
            max_relative_residual = relative_residual;
        }
        N_protons_total += N_protons;
        events << i << ',' << E_beam << ',' << pythia.info.code() << ',' << N_protons << '\n';
    }

    // Inclusive spectrum per injected event: integral over dE is the mean
    // number of final protons, not a probability normalized to one.
    for (int b = 0; b < N_bins; ++b) {
        const double E_low = std::exp(log_bin_min + b * dlogE);
        const double E_high = std::exp(log_bin_min + (b + 1) * dlogE);
        const double normalization = N_events * (E_high - E_low);
        spectrum << E_low << ',' << E_high << ',' << injected[b] << ',' << final_protons[b] << ','
                 << injected[b] / normalization << ',' << final_protons[b] / normalization << '\n';
    }
    events.close();
    protons.close();
    spectrum.close();
    if (!events || !protons || !spectrum) {
        return 1;
    }
    std::cout << std::defaultfloat << std::setprecision(12)
              << "\nIndependent inelastic events: " << N_events
              << "\nFinal protons: " << N_protons_total
              << "\nProtons from neutron beta decays: " << N_neutron_daughters
              << "\nMaximum total-energy residual [GeV]: " << max_energy_residual
              << "\nMaximum relative energy residual: " << max_relative_residual
              << "\nFailed generation attempts: " << N_failed_attempts << '\n';
    return 0;
}
