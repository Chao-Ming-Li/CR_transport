#include "CR_advection.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

namespace {
constexpr double center = 3.0;
constexpr double sigma = 0.4;
constexpr double velocity = 0.2;
constexpr double total_time = 10.0;

double gaussian_cell_average(double left, double right, double shift)
{
    const double scale = std::sqrt(2.0) * sigma;
    return sigma * std::sqrt(std::acos(-1.0) / 2.0) *
           (std::erf((right - center - shift) / scale) -
            std::erf((left - center - shift) / scale)) / (right - left);
}

struct Metrics { double mass, mean, variance, peak, relative_l1; };

Metrics measure(const Field2D& n, const Grid2D& grid, double shift)
{
    double mass = 0.0, first = 0.0, second = 0.0, peak = 0.0, l1 = 0.0, reference_mass = 0.0;
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) {
            const double z = grid.z().center(j), width = grid.z().width(j);
            const double volume = grid.volume(i,j), q = n(i,j);
            const double reference = gaussian_cell_average(grid.z().face(j), grid.z().face(j + 1), shift);
            mass += q * volume;
            first += q * volume * z;
            second += q * volume * (z * z + width * width / 12.0);
            peak = std::max(peak, q);
            l1 += std::abs(q - reference) * volume;
            reference_mass += reference * volume;
        }
    const double mean = first / mass;
    return {mass, mean, second / mass - mean * mean, peak, l1 / reference_mass};
}

void run(std::string_view name, const Grid2D& grid, const Field2D& initial,
         double dt, int steps)
{
    Field2D n = initial, temp(grid);
    Field2D vr(grid, Field2D::Location::RadialFace);
    Field2D vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 1; j <= grid.nz(); ++j) vz(i,j) = velocity;
    for (int step = 0; step < steps; ++step) {
        if (name == "upwind") advection_upwind(n, temp, vr, vz, dt, grid);
        else if (name == "TVD") advection_TVD(n, temp, vr, vz, dt, grid);
        else advection_PLM(n, temp, vr, vz, dt, grid);
    }
    const Metrics before = measure(initial, grid, 0.0);
    const Metrics after = measure(n, grid, velocity * total_time);
    const double effective_diffusion = (after.variance - before.variance) / (2.0 * total_time);
    std::cout << dt << ' ' << name << ' ' << after.mass / before.mass - 1.0 << ' '
              << after.mean - before.mean << ' ' << after.variance - before.variance << ' '
              << effective_diffusion << ' ' << after.peak << ' ' << after.relative_l1 << '\n';
}
} // namespace

int main()
{
    const Grid2D grid(AxisGrid::linear(2, 0.0, 1.0), AxisGrid::linear(160, 0.0, 10.0));
    Field2D initial(grid);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j)
            initial(i,j) = gaussian_cell_average(grid.z().face(j), grid.z().face(j + 1), 0.0);
    std::cout << std::setprecision(10);
    std::cout << "dt scheme mass_rel_change mean_shift variance_increase D_effective peak relative_L1\n";
    for (int steps : {100, 200})
        for (std::string_view name : {"upwind", "TVD", "PLM"})
            run(name, grid, initial, total_time / steps, steps);
}
