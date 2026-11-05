#include "Pythia8/Pythia.h"
#include <iostream>

int main() {
    Pythia8::Pythia pythia;
    // A 1000 GeV kinetic-energy proton collides with a stationary proton.
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:frameType = 2");
    pythia.readString("Beams:eA = 1000.9382720813");
    pythia.readString("Beams:eB = 0.9382720813");
    pythia.readString("SoftQCD:inelastic = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 12345");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    if (!pythia.init()) return 1;

    int N_events = 0;
    int N_photons = 0;
    for (int i = 0; i < 10; ++i) {
        if (!pythia.next()) continue;
        ++N_events;
        for (int j = 0; j < pythia.event.size(); ++j) {
            if (pythia.event[j].isFinal() && pythia.event[j].id() == 22) {
                ++N_photons;
            }
        }
    }
    std::cout << "Smoke test: " << N_events << "/10 events, " << N_photons << " final photons\n";
    return N_events == 10 && N_photons > 0 ? 0 : 1;
}
