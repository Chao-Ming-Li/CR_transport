#include "Pythia8/Pythia.h"
#include <iostream>

int main() {
    const double T_beam = 1000.0; // Projectile kinetic energy [GeV].
    Pythia8::Pythia pythia;
    const double m_p = pythia.particleData.m0(2212);
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:frameType = 2");
    pythia.settings.parm("Beams:eA", T_beam + m_p);
    pythia.settings.parm("Beams:eB", m_p);
    pythia.readString("SoftQCD:inelastic = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 12345");
    // Print explicitly below to avoid duplicate automatic event listings.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    if (!pythia.init()) return 1;
    if (!pythia.next()) return 1;

    std::cout << "\nSingle inelastic pp event: projectile kinetic energy 1000 GeV, stationary proton target\n";
    pythia.info.list();
    std::cout << "\nProcess record:\n";
    pythia.process.list(true, true);
    std::cout << "\nComplete event history (all particles, including decayed and intermediate particles):\n";
    pythia.event.list(true, true);
    pythia.stat();
    return 0;
}
