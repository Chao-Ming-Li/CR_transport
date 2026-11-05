# Local PYTHIA 8 installation

## Simple 1 TeV fixed-target example

`pp_fixed_target.cpp` generates 10000 independent inelastic collisions of a proton
with kinetic energy 1000 GeV on a stationary proton. Default decays are enabled.
Each event is written to `pp_fixed_target_events.csv` with exactly six columns:
`event_id,T_beam_GeV,nucleon1_pdg,T_nucleon1_GeV,nucleon2_pdg,T_nucleon2_GeV`.
All energies are laboratory kinetic energies in GeV. The CSV is overwritten each run.
From `external/`, run:

```sh
make
./pp_fixed_target
```

To save the terminal output, use `./pp_fixed_target > pp_fixed_target.log`.
`T_beam` at the beginning of the source sets the projectile kinetic energy;
the beam energy passed to PYTHIA includes the proton rest energy.
`N_events` sets the number of independent collisions.

The two retained nucleons are the highest-kinetic-energy final protons or
neutrons, excluding antinucleons. They are ordered by kinetic energy.
PDG 2212 denotes a proton and 2112 a neutron. Missing entries use PDG 0
and energy NaN. Neutron energies are before beta decay.

`pp_fixed_target_analysis.ipynb` preserves the previous executed proton-only
analysis. Its input cells require the old, wider CSV schema and cannot be
rerun directly on the new six-column output.

## Installation

Version: 8.318. Source: https://pythia.org/releases/pythia83/pythia8318.tgz

The source tree and local C++ libraries are in `external/pythia8318`.
Configured with `./configure --cxx=c++` and compiled with `make -j4`.
No optional external packages or shell startup changes are required.
Configuration and compilation logs are in that directory as
`configure-local.log` and `build-local.log`.

From the Monte_Carlo project directory, compile and run the small fixed-target
test with:

```sh
c++ external/pythia_smoke_test.cpp $(external/pythia8318/bin/pythia8-config --cxxflags --libs) -o external/pythia_smoke_test
./external/pythia_smoke_test
```

The test generates ten inelastic proton-proton events with a stationary target
and projectile kinetic energy 1000 GeV. It counts final photons after the
generator's default decays. This checks installation and event generation;
it does not validate the physical model across the full propagation energy range.

For another C++ source, use the same flags returned by `pythia8-config`.
They include the header path, library path and runtime library search path.
Keep this directory in place: the library also needs its XML data files.
The existing `pion_decay_mc.cpp` and project Makefile have not been changed.
