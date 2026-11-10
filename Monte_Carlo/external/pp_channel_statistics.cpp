#include "Pythia8/Pythia.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    const int N_events = argc > 1 ? std::stoi(argv[1]) : 10000;
    if (N_events <= 0) return 1;
    const double T_beam = 1000.0; // Projectile kinetic energy [GeV].
    Pythia8::Pythia pythia;
    const double m_p = pythia.particleData.m0(2212);
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:frameType = 2");
    pythia.settings.parm("Beams:eA", T_beam + m_p);
    pythia.settings.parm("Beams:eB", m_p);
    // Fractions below are conditional on an INELASTIC collision, not all pp.
    pythia.readString("SoftQCD:inelastic = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 12345");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    if (!pythia.init()) return 1;

    std::ofstream events("pp_channel_events.csv");
    std::ofstream particles("pp_channel_particles.csv");
    if (!events || !particles) return 1;
    events << "event_id,process_code,process_name,T_beam_GeV,N_final\n";
    particles << "event_id,particle_index,pid,name,E_GeV,T_GeV,px_GeV,py_GeV,pz_GeV\n";
    events << std::setprecision(16);
    particles << std::setprecision(16);
    int N_failed = 0;
    for (int i = 0; i < N_events; ++i) {
        bool generated = false;
        for (int attempt = 0; attempt < 10; ++attempt) {
            if (pythia.next()) {
                generated = true;
                break;
            }
            ++N_failed;
        }
        if (!generated) {
            std::cerr << "Generation failed at event " << i << '\n';
            return 1;
        }
        int N_final = 0;
        for (int j = 0; j < pythia.event.size(); ++j) {
            const Pythia8::Particle& particle = pythia.event[j];
            // Final means after enabled PYTHIA decays, not infinite-time decay.
            if (!particle.isFinal()) continue;
            ++N_final;
            particles << i << ',' << j << ',' << particle.id() << ",\"" << particle.name() << "\","
                      << particle.e() << ',' << particle.e() - particle.m() << ','
                      << particle.px() << ',' << particle.py() << ',' << particle.pz() << '\n';
        }
        events << i << ',' << pythia.info.code() << ",\"" << pythia.info.name() << "\","
               << T_beam << ',' << N_final << '\n';
    }
    events.close();
    particles.close();
    if (!events || !particles) return 1;
    pythia.stat();
    std::cout << "\nSaved " << N_events << " events; failed generation attempts = " << N_failed << '\n';
    return 0;
}
