#include "momentum_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Driver for energy_loss.ipynb: cells, steps, time [yr], a [yr^-1],
// profile (peak/power/cutoff), boundary (outflow/closed), wind (both/radial/vertical).
// Initial values are true cell averages, not samples at geometric centers.
int main(int argc, char **argv) {
    if (argc != 8) return 1;
    const int cells = std::stoi(argv[1]), steps = std::stoi(argv[2]);
    const double time = std::stod(argv[3]), a = std::stod(argv[4]);
    const std::string profile = argv[5], wind = argv[7];
    const auto boundary = std::string(argv[6]) == "closed"
        ? cr_momentum::Boundary::ZeroFlux : cr_momentum::Boundary::Outflow;
    const SpatialGrid grid(AxisGrid::from_faces({0.0, 0.3, 2.0}),
                           AxisGrid::from_faces({-2.0, -0.4, 2.0}));
    const auto Rg = RgGrid::logarithmic(cells + 1, 1.0, 1e6);
    Field2D vr(grid, Field2D::Location::RadialFace);
    Field2D vz(grid, Field2D::Location::VerticalFace);
    // All three winds have cylindrical divergence 3*a, including the axis.
    const double ar = wind == "vertical" ? 0.0 : (wind == "radial" ? 1.5*a : a);
    const double az = wind == "radial" ? 0.0 : (wind == "vertical" ? 3.0*a : a);
    for (int i = 0; i < vr.nR(); ++i)
        for (int j = 0; j < vr.nz(); ++j) vr(i,j) = ar * grid.R().face(i);
    for (int i = 0; i < vz.nR(); ++i)
        for (int j = 0; j < vz.nz(); ++j) vz(i,j) = az * grid.z().face(j);
    const auto integral = [&](double lo, double hi) {
        if (profile == "power") return 1.0/lo - 1.0/hi;
        if (profile == "cutoff") return std::log(hi/lo);
        // Gaussian in ln(Rg): Rg*q peaks at Rg=1000 GV, sigma_lnRg=0.5.
        const double scale = std::sqrt(2.0)*0.5;
        return 0.5 * (std::erf(std::log(hi/1000.0)/scale)
                    - std::erf(std::log(lo/1000.0)/scale));
    };
    Field3D q(grid, Rg);
    for (int k = 0; k < cells; ++k)
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                q(k,i,j) = integral(Rg.face(k), Rg.face(k+1))/Rg.width(k);
    const auto initial = q;
    q(0,-1,-1) = 7.0;
    const std::vector<double> Dpp(cells + 1, 0.0);
    cr_momentum::Solver solver(grid, Rg, vr, vz, Dpp, time/steps, boundary, boundary);
    const auto outgoing = [&]() {
        if (boundary == cr_momentum::Boundary::ZeroFlux) return 0.0;
        return std::max(a, 0.0)*Rg.face(0)*q(0,0,0)
             + std::max(-a, 0.0)*Rg.face(cells)*q(cells-1,0,0);
    };
    double escaped = 0.0;
    for (int t = 0; t < steps; ++t) {
        const double old_flux = outgoing();
        solver.advance(q);
        escaped += 0.5*time/steps*(old_flux + outgoing());
    }
    double N0 = 0.0, N = 0.0, spatial_error = 0.0;
    for (int k = 0; k < cells; ++k) {
        N0 += initial(k,0,0)*Rg.width(k);
        N += q(k,0,0)*Rg.width(k);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                spatial_error = std::max(spatial_error, std::abs(q(k,i,j)-q(k,0,0)));
    }
    std::cout << std::setprecision(17);
    std::cout << "# " << N0 << ' ' << N << ' ' << escaped << ' '
              << spatial_error << ' ' << q(0,-1,-1) << '\n';
    std::cout << "# lo hi Rg initial final\n";
    for (int k = 0; k < cells; ++k)
        std::cout << Rg.face(k) << ' ' << Rg.face(k+1) << ' ' << Rg.center(k)
                  << ' ' << initial(k,0,0) << ' ' << q(k,0,0) << '\n';
}
