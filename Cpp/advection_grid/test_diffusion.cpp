#include "CR_advection.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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

    cr_diffusion::Options reflecting() {
        return {cr_diffusion::Boundary::Reflecting, cr_diffusion::Boundary::Reflecting};
    }

    void test_nonuniform_conservation() {
        const Grid2D grid(AxisGrid::geometric(12, 0.0, 3.0, 1.15),
                          AxisGrid::geometric(15, 0.0, 2.0, 1.1));
        Field2D density(grid);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                density(i, j) = (i + 2 * j) % 4 == 0 ? 1.0 : 0.0;
        const double initial_mass = mass(density, grid);
        cr_diffusion::Solver solver(grid, 0.1, 0.001, reflecting());
        for (int step = 0; step < 50; ++step)
            solver.advance(density);
        check(std::abs(mass(density, grid) / initial_mass - 1.0) < 2e-13,
              "Nonuniform closed-domain mass is not conserved");
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                check(density(i, j) >= 0.0, "Small-step diffusion became negative");

        for (double &value : density.data)
            value = 3.0;
        cr_diffusion::Solver(grid, 0.1, 2.0, reflecting()).advance(density);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                check(std::abs(density(i, j) - 3.0) < 1e-13,
                      "Reflecting diffusion changed a constant field");
    }

    // Vertical cosine is an exact discrete eigenmode on a uniform reflecting grid.
    // This independently checks the CN amplification and second-order time error.
    double eigenmode_error(int steps) {
        const int nz = 16;
        const Grid2D grid(AxisGrid::geometric(6, 0.0, 2.0, 1.2), AxisGrid::linear(nz, 0.0, 1.0));
        Field2D density(grid);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < nz; ++j)
                density(i, j) = 2.0 + std::cos(std::numbers::pi * (j + 0.5) / nz);
        const double D = 0.1, duration = 0.2, dt = duration / steps;
        const double lambda =
            -4.0 * D * nz * nz * std::pow(std::sin(std::numbers::pi / (2.0 * nz)), 2);
        cr_diffusion::Solver solver(grid, D, dt, reflecting());
        for (int step = 0; step < steps; ++step)
            solver.advance(density);
        const double amplification =
            std::pow((1.0 + 0.5 * dt * lambda) / (1.0 - 0.5 * dt * lambda), steps);
        double error = 0.0;
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < nz; ++j) {
                const double mode = std::cos(std::numbers::pi * (j + 0.5) / nz);
                check(std::abs(density(i, j) - (2.0 + amplification * mode)) < 2e-13,
                      "ADI does not match discrete CN eigenmode amplification");
                error = std::max(
                    error, std::abs(density(i, j) - (2.0 + std::exp(lambda * duration) * mode)));
            }
        return error;
    }

    void test_absorbing_and_validation() {
        const Grid2D grid(AxisGrid::linear(5, 0.0, 2.0), AxisGrid::linear(7, 0.0, 3.0));
        Field2D density(grid);
        for (double &value : density.data)
            value = 1.0;
        const Field2D initial = density;
        cr_diffusion::Solver(grid, 0.0, 1.0).advance(density);
        check(density.data == initial.data, "Zero D must be a no-op");
        cr_diffusion::Solver(grid, 1.0, 0.0).advance(density);
        check(density.data == initial.data, "Zero dt must be a no-op");
        for (double D : {-1.0, std::numeric_limits<double>::infinity()}) {
            bool rejected = false;
            try {
                cr_diffusion::Solver(grid, D, 0.001).advance(density);
            } catch (const std::invalid_argument &) {
                rejected = true;
            }
            check(rejected && density.data == initial.data, "Invalid D changed density");
        }
        const double initial_mass = mass(density, grid);
        cr_diffusion::Solver(grid, 0.1, 0.001).advance(density);
        check(mass(density, grid) < initial_mass, "Absorbing boundaries did not remove mass");
        check(density(grid.nR(), 2) == 0.0 && density(2, grid.nz()) == 0.0,
              "Absorbing ghost cells are not zero");
        check(density(-1, 2) == density(0, 2) && density(2, -1) == density(2, 0),
              "Lower ghost cells are not reflecting");
    }

    void test_two_dimensional_mode() {
        const int nz = 8;
        const Grid2D grid(AxisGrid::linear(2, 0.0, 2.0), AxisGrid::linear(nz, 0.0, 1.0));
        Field2D density(grid);
        const double D = 0.1, dt = 0.02;
        // For these two annuli, L_R has eigenvector (1, -1/3) and eigenvalue -3D.
        // Their volume centroids are 2/3 and 14/9, with spacing 8/9.
        const double radial_mode[] = {1.0, -1.0 / 3.0};
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < nz; ++j)
                density(i, j) = 2.0 + radial_mode[i] * std::cos(std::numbers::pi * (j + 0.5) / nz);
        const double lambda_R = -3.0 * D;
        const double lambda_z =
            -4.0 * D * nz * nz * std::pow(std::sin(std::numbers::pi / (2.0 * nz)), 2);
        const auto amplification = [dt](double lambda) {
            return (1.0 + 0.5 * dt * lambda) / (1.0 - 0.5 * dt * lambda);
        };
        cr_diffusion::Solver(grid, D, dt, reflecting()).advance(density);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < nz; ++j) {
                const double expected = 2.0 + amplification(lambda_R) * amplification(lambda_z) *
                                                  radial_mode[i] *
                                                  std::cos(std::numbers::pi * (j + 0.5) / nz);
                check(std::abs(density(i, j) - expected) < 2e-14,
                      "Two-dimensional ADI amplification is incorrect");
            }
    }
} // namespace

int main() {
    test_nonuniform_conservation();
    test_absorbing_and_validation();
    test_two_dimensional_mode();
    const double coarse = eigenmode_error(4), fine = eigenmode_error(8);
    check(coarse / fine > 3.9 && coarse / fine < 4.1,
          "Diffusion is not second-order accurate in time");
    std::cout << "Diffusion tests passed; temporal error ratio = " << coarse / fine << '\n';
}
