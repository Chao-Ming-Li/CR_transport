#include "cross_section.hpp"
#include "nuclear_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
    void check(bool condition) {
        if (!condition) throw std::runtime_error("Nuclear solver test failed");
    }

    double reaction_rate(cross_section::Isotope parent, cross_section::Isotope daughter, double Rg, double nH) {
        const cross_section::Nuclide &n = cross_section::nuclide(parent);
        const double p = Rg * n.charge / n.mass_number;
        const double beta = p / std::hypot(p, physical_constants::nucleon_rest_energy_GeV);
        return nH * beta * physical_constants::speed_of_light_cm_s * physical_constants::seconds_per_year *
               cross_section::sigma[cross_section::index(parent)][cross_section::index(daughter)];
    }
}

int main() {
    using cross_section::Isotope;
    const SpatialGrid grid(AxisGrid::linear(2, 0.0, 2.0), AxisGrid::linear(2, 0.0, 2.0));
    const RgGrid Rg = RgGrid::logarithmic(41, 1.0, 100.0);
    Field2D gas(grid);
    for (double &value : gas.data) value = 1e6;
    const int k0 = 15;

    // Equal A/Z isolates a two-species analytic decay/production problem.
    const std::vector<Isotope> isotopes = {Isotope::O16, Isotope::C12};
    const double a = reaction_rate(isotopes[0], isotopes[0], Rg.center(k0), gas(0,0));
    const double b = reaction_rate(isotopes[0], isotopes[1], Rg.center(k0), gas(0,0));
    const double d = reaction_rate(isotopes[1], isotopes[1], Rg.center(k0), gas(0,0));
    const double exact_parent = std::exp(-a);
    const double exact_daughter = b * (std::exp(-a) - std::exp(-d)) / (d - a);
    double errors[2];
    for (int resolution = 0; resolution < 2; ++resolution) {
        const int steps = 100 * (resolution + 1);
        std::vector<Field3D> q(2, Field3D(grid, Rg));
        q[0](k0,0,0) = 1.0;
        q[0](k0,-1,-1) = 7.0;
        cr_nuclear::Solver solver(grid, Rg, gas, isotopes, 1.0 / steps, 1.0);
        for (int t = 0; t < steps; ++t) {
            solver.advance(q);
        }
        errors[resolution] = std::abs(q[0](k0,0,0) - exact_parent) + std::abs(q[1](k0,0,0) - exact_daughter);
        check(q[0](k0,-1,-1) == 7.0);
        for (const Field3D &field : q) {
            for (double value : field.data()) {
                check(value >= 0.0);
            }
        }
        if (resolution == 1) check(errors[resolution] < 1e-3);
    }
    check(errors[0] / errors[1] > 1.9 && errors[0] / errors[1] < 2.1);

    // C12 -> B11 shifts rigidity by alpha=1.1. Integrated production must
    // equal the parent-bin reaction count, with only out-of-domain mass lost.
    const std::vector<Isotope> carbon_boron = {Isotope::C12, Isotope::B11};
    for (int parent_bin : {k0, Rg.size() - 1}) {
        std::vector<Field3D> q(2, Field3D(grid, Rg));
        std::vector<Field3D> production(2, Field3D(grid, Rg)), loss(2, Field3D(grid, Rg));
        q[0](parent_bin,0,0) = 1.0;
        cr_nuclear::Solver solver(grid, Rg, gas, carbon_boron, 0.1, 1.0);
        solver.calculate_terms(q, production, loss);
        double produced = 0.0;
        for (int k = 0; k < Rg.size(); ++k) {
            produced += production[1](k,0,0) * Rg.width(k);
            if (Rg.face(k + 1) <= 1.1 * Rg.face(parent_bin) || Rg.face(k) >= 1.1 * Rg.face(parent_bin + 1))
                check(production[1](k,0,0) == 0.0);
        }
        const double retained_width = (std::min(100.0, 1.1 * Rg.face(parent_bin + 1)) - 1.1 * Rg.face(parent_bin)) / 1.1;
        const double rate = reaction_rate(carbon_boron[0], carbon_boron[1], Rg.center(parent_bin), gas(0,0));
        check(std::abs(produced / (rate * retained_width) - 1.0) < 1e-13);
        check(std::abs(loss[0](parent_bin,0,0) / reaction_rate(carbon_boron[0], carbon_boron[0], Rg.center(parent_bin), gas(0,0)) - 1.0) < 1e-13);
    }

    // Three-species chain, deliberately reversed input order.
    const std::vector<Isotope> chain = {Isotope::O16, Isotope::C12, Isotope::B11};
    const std::vector<Isotope> reversed = {Isotope::B11, Isotope::C12, Isotope::O16};
    std::vector<Field3D> q(3, Field3D(grid, Rg)), reverse_q(3, Field3D(grid, Rg));
    q[0](k0,0,0) = reverse_q[2](k0,0,0) = 1.0;
    cr_nuclear::Solver forward(grid, Rg, gas, chain, 0.1, 1.0), backward(grid, Rg, gas, reversed, 0.1, 1.0);
    forward.advance(q);
    backward.advance(reverse_q);
    for (int s = 0; s < 3; ++s) {
        // Different input order can change the roundoff order of source sums.
        for (std::size_t l = 0; l < q[s].data().size(); ++l) {
            check(std::abs(q[s].data()[l] - reverse_q[2-s].data()[l]) < 1e-14);
        }
    }
    check(q[2](k0,0,0) > 0.0);

    // No gas leaves the full network unchanged, including ghosts.
    for (double &value : gas.data) value = 0.0;
    const std::vector<Field3D> before = q;
    cr_nuclear::Solver vacuum(grid, Rg, gas, chain, 1e10, 1.0);
    vacuum.advance(q);
    for (int s = 0; s < 3; ++s) {
        check(std::equal(q[s].data().begin(), q[s].data().end(), before[s].data().begin()));
    }
    std::vector<Field3D> proton(1, Field3D(grid, Rg));
    proton[0](k0,0,0) = 1.0;
    for (double &value : gas.data) value = 1.0;
    cr_nuclear::Solver proton_solver(grid, Rg, gas, {Isotope::Proton}, 1e10, 1.0);
    proton_solver.advance(proton);
    check(proton[0](k0,0,0) == 1.0);
    std::cout << "Nuclear: analytic production/loss, first-order convergence, rigidity mapping, domain escape, chain ordering, vacuum and proton passed.\n";
}
