#include "cross_section.hpp"
#include "spatial_solver.hpp"
#include "momentum_solver.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
    void check(bool ok) {
        if (!ok) throw std::runtime_error("Coupled diffusion test failed");
    }
}

int main() {
    const double D0 = 1e-7, delta = 1.0 / 3.0, vA = 30.0;
    const cr_diffusion::DiffusionModel model(D0, delta, vA);
    constexpr double speed = 30.0 * physical_constants::km_s_to_kpc_yr;
    // delta=1/3 gives the analytic coefficient 108/385.
    for (double Rg : {1.0, 10.0, 1000.0}) {
        check(std::abs(model.spatial_coefficient(Rg) / (D0 * std::cbrt(Rg)) - 1.0) < 1e-14);
        check(std::abs(model.rigidity_coefficient(Rg) * model.spatial_coefficient(Rg) /
                       (Rg * Rg * speed * speed) - 108.0 / 385.0) < 1e-14);
        check(cr_diffusion::DiffusionModel(D0, delta, 0.0).rigidity_coefficient(Rg) == 0.0);
        check(std::abs(cr_diffusion::DiffusionModel(D0, delta, 60.0).rigidity_coefficient(Rg) /
                       model.rigidity_coefficient(Rg) - 4.0) < 1e-14);
    }
    const SpatialGrid grid(AxisGrid::geometric(8, 0.0, 2.0, 1.1), AxisGrid::linear(6, 0.0, 2.0));
    const auto Rg = RgGrid::logarithmic(9, 1.0, 1000.0);
    Field3D q(grid, Rg);
    for (int k = 0; k < Rg.size(); ++k) {
        for (int i = 0; i < grid.nR(); ++i) {
            for (int j = 0; j < grid.nz(); ++j) {
                q(k,i,j) = std::exp(-grid.radial_centroid(i) * grid.radial_centroid(i)
                                                  - grid.z().center(j) * grid.z().center(j)) / Rg.center(k);
            }
        }
    }
    auto reference = q;
    cr_diffusion::Solver spatial(grid, Rg, model, 1000.0);
    spatial.advance(q);
    Field2D plane(grid);
    for (int k = 0; k < Rg.size(); ++k) {
        auto values = reference.spatial_plane(k);
        std::copy(values.begin(), values.end(), plane.data.begin());
        cr_diffusion::Solver single(grid, D0 * std::pow(Rg.center(k), delta), 1000.0);
        single.advance(plane);
        check(std::equal(plane.data.begin(), plane.data.end(), q.spatial_plane(k).begin()));
    }
    Field2D vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    reference = q;
    double number0 = 0.0;
    for (int k = 0; k < Rg.size(); ++k) number0 += q(k,0,0) * Rg.width(k);
    cr_momentum::Solver momentum(grid, Rg, vr, vz, model, 1000.0);
    momentum.advance(q);
    double number = 0.0;
    for (int k = 0; k < Rg.size(); ++k) number += q(k,0,0) * Rg.width(k);
    check(std::abs(number / number0 - 1.0) < 1e-13);
    std::cout << "Coupled diffusion: normalization, Alfvén scaling, rigidity planes and momentum conservation passed.\n";
}
