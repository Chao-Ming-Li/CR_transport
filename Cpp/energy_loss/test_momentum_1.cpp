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
    std::vector<double> Dpp(Rg.size() + 1, 0.0);
    Field3D q(grid, Rg);
    for (int k = 0; k < Rg.size(); ++k)
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                // q(k, i, j) = std::pow(Rg.center(k), -1.0); // power law spectrum
                q(k, i, j) = std::exp(-std::pow((Rg.center(k) - 20.0) / (std::sqrt(2.0) * 2.0), 2.0)); // Gaussian spectrum
    auto const initial = q;

    // test adiabatic cooling: vR=a*R, vz=a*z => cylindrical divergence=3*a, including R=0.
    // const double expansion = 0.1;
    // for (int i = 0; i < vr.nR(); ++i)
    //     for (int j = 0; j < vr.nz(); ++j)
    //         vr(i,j) = expansion * grid.R().face(i);
    // for (int i = 0; i < vz.nR(); ++i)
    //     for (int j = 0; j < vz.nz(); ++j)
    //         vz(i,j) = expansion * grid.z().face(j);

    // cr_momentum::Solver cooling(grid, Rg, vr, vz, Dpp, 0.01, cr_momentum::Boundary::Outflow, cr_momentum::Boundary::Outflow);

    // for (int t = 0; t < 1000; ++t)
    //     cooling.advance(q);
    // test stochastic diffusion: constant Dpp, vR=vz=0; zero-flux equilibrium q/Rg^2=constant.
    for (int i = 0; i < vr.nR(); ++i)
        for (int j = 0; j < vr.nz(); ++j)
            vr(i, j) = 0.0;
    for (int i = 0; i < vz.nR(); ++i)
        for (int j = 0; j < vz.nz(); ++j)
            vz(i, j) = 0.0;

    for (int k = 0; k < Dpp.size(); ++k)
        Dpp[k] = 1.0;
    cr_momentum::Solver cooling(grid, Rg, vr, vz, Dpp, 0.001, cr_momentum::Boundary::Outflow, cr_momentum::Boundary::Outflow);

    for (int t = 0; t < 100000; ++t)
        cooling.advance(q);

    // Native double records: Rg [GV], initial dn/dRg, accelerated dn/dRg.
    // All spatial cells evolve identically; save the spectrum at (i,j)=(0,0).
    std::ofstream output("momentum_spectra_stochastic_t100.bin", std::ios::binary | std::ios::trunc);

    for (int k = 0; k < Rg.size(); ++k) {
        const double record[] = {Rg.center(k), initial(k, 0, 0), q(k, 0, 0)};
        output.write(reinterpret_cast<const char *>(record), sizeof(record));
    }
    output.close();
}
