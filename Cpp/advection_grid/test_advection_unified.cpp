#include "CR_advection_unified.hpp"
#include "CR_advection.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cr_advection;

void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}

double mass(const Field2D &n, const Grid2D &g) {
    double total = 0;
    for (int i = 0; i < g.nR(); ++i)
        for (int j = 0; j < g.nz(); ++j)
            total += n(i, j) * g.volume(i, j);
    return total;
}

int main() {
    for (bool stretched : {false, true}) {
        Grid2D g(stretched ? AxisGrid::geometric(12, 0, 12, 1.08) : AxisGrid::linear(12, 0, 12),
                 stretched ? AxisGrid::geometric(14, 0, 14, 1.06) : AxisGrid::linear(14, 0, 14));
        Field2D initial(g), vr(g, Field2D::Location::RadialFace),
            vz(g, Field2D::Location::VerticalFace), temp(g);
        for (int i = 0; i < g.nR(); ++i)
            for (int j = 0; j < g.nz(); ++j)
                initial(i, j) = 1.0 + 0.2 * std::sin(0.4 * i) * std::cos(0.3 * j);
        // Both velocity signs, with closed boundaries for mass conservation.
        for (int i = 1; i < g.nR(); ++i)
            for (int j = 0; j < g.nz(); ++j)
                vr(i, j) = 0.2 * std::cos(0.4 * i + 0.3 * j);
        for (int i = 0; i < g.nR(); ++i)
            for (int j = 1; j < g.nz(); ++j)
                vz(i, j) = 0.2 * std::sin(0.5 * i - 0.2 * j);

        for (int method = 0; method < 4; ++method) {
            Options opt;
            auto old = initial, unified = initial;
            apply_boundary_conditions(old, g);
            if (method == 0) {
                opt.reconstruction = Reconstruction::Constant;
                opt.integrator = TimeIntegrator::Euler;
                advection_upwind(old, temp, vr, vz, 0.01, g);
            } else if (method == 1) {
                opt.limiter = Limiter::MC;
                opt.integrator = TimeIntegrator::Euler;
                opt.splitting = Splitting::RadialThenVertical;
                advection_PLM(old, temp, vr, vz, 0.01, g);
            } else if (method == 2) {
                opt.integrator = TimeIntegrator::Euler;
                opt.splitting = Splitting::RadialThenVertical;
                advection_TVD(old, temp, vr, vz, 0.01, g);
            } else {
                advance_TVD_SSPRK2(old, temp, vr, vz, 0.01, g);
            }
            Solver(g, opt).advance(unified, vr, vz, 0.01);
            double error = 0;
            for (std::size_t k = 0; k < old.data.size(); ++k)
                error = std::max(error, std::abs(old.data[k] - unified.data[k]));
            require(error < 1e-10, "Legacy comparison failed");
        }

        for (auto reconstruction : {Reconstruction::Constant, Reconstruction::PLM})
            for (auto limiter : {Limiter::Minmod, Limiter::VanLeer, Limiter::MC})
                for (auto integrator : {TimeIntegrator::Euler, TimeIntegrator::SSPRK2}) {
                    Options opt{reconstruction, limiter, integrator, Splitting::Unsplit};
                    auto n = initial;
                    Solver solver(g, opt);
                    for (int k = 0; k < 10; ++k)
                        solver.advance(n, vr, vz, 0.01);
                    require(std::abs(mass(n, g) / mass(initial, g) - 1) < 1e-13,
                            "Closed-boundary mass conservation failed");
                }

        auto n = initial;
        Solver solver(g);
        solver.advance(n, vr, vz, 0.0);
        auto expected = initial;
        apply_boundary_conditions(expected, g);
        require(n.data == expected.data, "Zero timestep changed physical density");
        bool rejected = false;
        try {
            solver.advance(n, vr, vz, -0.1);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Negative timestep accepted");
        require(n.data == expected.data, "Rejected step modified density");

        // Excessive combined outflow must fail without partially updating n.
        Field2D spike(g), outR(g, Field2D::Location::RadialFace),
            outZ(g, Field2D::Location::VerticalFace);
        spike(4, 4) = 1;
        outR(5, 4) = 1;
        outZ(4, 5) = 1;
        const auto before = spike.data;
        rejected = false;
        try {
            solver.advance(spike, outR, outZ, 10.0);
        } catch (const std::runtime_error &) {
            rejected = true;
        }
        require(rejected && spike.data == before, "Failed stage was not rejected atomically");
    }
    std::cout << "Unified advection: legacy agreement, all limiter/integrator combinations,\n"
                 "mass conservation, ghost refresh, and invalid-step checks passed.\n";
}
