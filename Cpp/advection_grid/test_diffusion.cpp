#include "CR_advection.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
    void check(bool condition, const char *message) {
        if (!condition)
            throw std::runtime_error(message);
    }

    double mass(const Field2D &density, const Grid2D &grid) {
        double total = 0.0;
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                total += density(i, j) * grid.volume(i, j);
        return total;
    }

    void test_boundaries() {
        const Grid2D grid(AxisGrid::geometric(12, 0.0, 3.0, 1.15),
                          AxisGrid::geometric(15, 0.0, 2.0, 1.1));
        Field2D density(grid);
        for (double &value : density.data)
            value = 1.0;
        const auto initial = density.data;
        cr_diffusion::Solver(grid, 0.0, 1.0).advance(density);
        check(density.data == initial, "Zero D must be a no-op");
        cr_diffusion::Solver solver(grid, 0.1, 0.001);
        double previous_mass = mass(density, grid);
        for (int step = 0; step < 50; ++step) {
            solver.advance(density);
            const double current_mass = mass(density, grid);
            check(current_mass < previous_mass, "Absorbing boundaries must remove mass");
            previous_mass = current_mass;
            for (int i = 0; i < grid.nR(); ++i)
                for (int j = 0; j < grid.nz(); ++j)
                    check(density(i, j) >= 0.0, "Small-step diffusion became negative");
        }
        check(density(grid.nR(), 2) == 0.0 && density(2, grid.nz()) == 0.0,
              "Outer ghost cells must absorb");
        check(density(-1, 2) == density(0, 2) && density(2, -1) == density(2, 0),
              "Lower ghost cells must reflect");
    }

    // Exact discrete eigenmode for reflecting lower / zero-ghost outer boundaries.
    double eigenmode_error(int steps) {
        const int nz = 16;
        const Grid2D grid(AxisGrid::linear(2, 0.0, 2.0), AxisGrid::linear(nz, 0.0, 1.0));
        const double D = 0.1, duration = 0.2, dt = duration / steps;
        // Radial operator / D = {{-9/4, 9/4}, {3/4, -9/4}}.
        const double radial_mode[] = {1.0, 1.0 / std::sqrt(3.0)};
        const double theta = std::numbers::pi / (2.0 * nz + 1.0);
        const double lambda_R = D * (-2.25 + std::sqrt(27.0) / 4.0);
        const double lambda_z = -4.0 * D * nz * nz * std::pow(std::sin(theta / 2.0), 2);
        Field2D density(grid);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < nz; ++j)
                density(i, j) = radial_mode[i] * std::cos(theta * (j + 0.5));
        cr_diffusion::Solver solver(grid, D, dt);
        for (int step = 0; step < steps; ++step)
            solver.advance(density);
        const auto amplification = [dt](double lambda) {
            return (1.0 + 0.5 * dt * lambda) / (1.0 - 0.5 * dt * lambda);
        };
        const double discrete = std::pow(amplification(lambda_R) * amplification(lambda_z), steps);
        const double exact = std::exp((lambda_R + lambda_z) * duration);
        double error = 0.0;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < nz; ++j) {
                const double mode = radial_mode[i] * std::cos(theta * (j + 0.5));
                check(std::abs(density(i, j) - discrete * mode) < 2e-13,
                      "Incorrect two-dimensional ADI amplification");
                error = std::max(error, std::abs(density(i, j) - exact * mode));
            }
        return error;
    }
} // namespace

int main() {
    test_boundaries();
    const double coarse = eigenmode_error(4), fine = eigenmode_error(8);
    check(coarse / fine > 3.9 && coarse / fine < 4.1, "Diffusion must be second order in time");
    std::cout << "Diffusion tests passed; temporal error ratio = " << coarse / fine << '\n';
}
