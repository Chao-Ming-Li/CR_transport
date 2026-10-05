#include "CR_advection.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

void check(bool ok) {
    if (!ok)
        throw std::runtime_error("Reusable solver test failed");
}
int main() {
    const Grid2D grid(AxisGrid::geometric(9, 0, 3, 1.1), AxisGrid::geometric(12, 0, 2, 1.12));
    Field2D density(grid), vr(grid, Field2D::Location::RadialFace),
        vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j)
            density(i, j) = std::exp(-std::pow(grid.radial_centroid(i) - 1, 2) -
                                     std::pow(grid.z().center(j) - 0.7, 2));
    for (double &v : vr.data)
        v = 0.1;
    for (double &v : vz.data)
        v = -0.1;
    cr_advection::Options options;
    options.limiter = cr_advection::Limiter::MC;
    cr_advection::Solver advection(grid, options);
    cr_diffusion::Solver diffusion(grid, 0.02, 0.001);
    Field2D reference = density;
    for (int step = 0; step < 20; ++step) {
        const double dt = step % 2 ? 0.001 : 0.002;
        advection.advance(density, vr, vz, dt);
        diffusion.advance(density);
        cr_advection::Solver fresh_advection(grid, options);
        fresh_advection.advance(reference, vr, vz, dt);
        cr_diffusion::Solver fresh_diffusion(grid, 0.02, 0.001);
        fresh_diffusion.advance(reference);
        for (std::size_t k = 0; k < density.data.size(); ++k)
            check(std::abs(density.data[k] - reference.data[k]) < 1e-14);
    }
    auto before = density.data;
    bool rejected = false;
    Field2D wrong(grid);
    try {
        advection.advance(density, wrong, vz, 0.001);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected && before == density.data);
    cr_diffusion::Solver disabled(grid, 0, 1);
    disabled.advance(density);
    check(before == density.data);
    std::cout << "Reusable solver tests passed\n";
}
