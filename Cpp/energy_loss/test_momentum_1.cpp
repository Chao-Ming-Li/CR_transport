#include "momentum_solver.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    const SpatialGrid grid(AxisGrid::linear(2, 0.0, 2.0), AxisGrid::linear(2, 0.0, 2.0));
    const auto Rg = RgGrid::logarithmic(513, 1.0, 100.0);
    Field2D vr(grid, Field2D::Location::RadialFace);
    Field2D vz(grid, Field2D::Location::VerticalFace);
    const cr_diffusion::DiffusionModel model(1e-7, 1.0 / 3.0, 30.0);
    Field3D q(grid, Rg);
    for (int k = 0; k < Rg.size(); ++k) {
        for (int i = 0; i < grid.nR(); ++i) {
            for (int j = 0; j < grid.nz(); ++j) {
                // q(k, i, j) = std::pow(Rg.center(k), -1.0); // power law spectrum
                q(k, i, j) = std::exp(-std::pow((Rg.center(k) - 20.0) / (std::sqrt(2.0) * 2.0), 2.0));
            }
        }
    } // Gaussian spectrum
    auto const initial = q;

    // Pure reacceleration from the linked spatial/rigidity diffusion model.
    // D0=1e-7 kpc^2/yr at 1 GV, delta=1/3, vA=30 km/s.
    for (int i = 0; i < vr.nR(); ++i) {
        for (int j = 0; j < vr.nz(); ++j) {
            vr(i, j) = 0.0;
        }
    }
    for (int i = 0; i < vz.nR(); ++i) {
        for (int j = 0; j < vz.nz(); ++j) {
            vz(i, j) = 0.0;
        }
    }

    cr_momentum::Solver cooling(grid, Rg, vr, vz, model, 1000.0, cr_momentum::Boundary::Outflow, cr_momentum::Boundary::Outflow);

    // 100 Myr with a 1000 yr timestep.
    for (int t = 0; t < 100000; ++t) {
        cooling.advance(q);
    }

    // Native double records: Rg [GV], initial dn/dRg, accelerated dn/dRg.
    // All spatial cells evolve identically; save the spectrum at (i,j)=(0,0).
    std::ofstream output("momentum_spectra_linked_t100Myr.bin", std::ios::binary | std::ios::trunc);

    for (int k = 0; k < Rg.size(); ++k) {
        const double record[] = {Rg.center(k), initial(k, 0, 0), q(k, 0, 0)};
        output.write(reinterpret_cast<const char *>(record), sizeof(record));
    }
    output.close();
}
