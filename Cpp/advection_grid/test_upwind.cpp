#include "CR_advection.hpp"

#include <cmath>
#include <stdexcept>

namespace {
void check(bool condition) { if (!condition) throw std::runtime_error("Upwind test failed"); }

double mass(const Field2D& n, const Grid2D& grid)
{
    double result = 0.0;
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) result += n(i,j) * grid.volume(i,j);
    return result;
}

void closed_grid_test()
{
    const Grid2D grid(AxisGrid::geometric(8, 0.0, 2.0, 1.2),
                      AxisGrid::geometric(9, 0.0, 3.0, 1.15));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace),
            vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) n(i,j) = ((i + 2 * j) % 4 == 0) ? 1.0 : 0.0;
    for (int i = 0; i <= grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) vr(i,j) = i == grid.nR() ? 0.0 : (i % 2 ? 0.3 : -0.3);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j <= grid.nz(); ++j) vz(i,j) = j == grid.nz() ? 0.0 : (j % 2 ? -0.25 : 0.25);
    const Field2D vr_before = vr, vz_before = vz;
    const double initial = mass(n, grid);
    bool cfl_rejected = false;
    try { advection_upwind(n, temp, vr, vz, 2.0, grid); }
    catch (const std::invalid_argument&) { cfl_rejected = true; }
    check(cfl_rejected);
    const double unchanged_mass = mass(n, grid);
    check(unchanged_mass == initial);
    for (int step = 0; step < 10; ++step)
        advection_upwind(n, temp, vr, vz, 0.02, grid);
    check(std::abs(mass(n, grid) - initial) < 1e-12 * initial);
    check(vr.data == vr_before.data && vz.data == vz_before.data);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) check(std::isfinite(n(i,j)) && n(i,j) >= 0.0);
}

void directional_and_boundary_test()
{
    const Grid2D grid(AxisGrid::linear(2, 0.0, 2.0), AxisGrid::linear(3, 0.0, 3.0));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace),
            vz(grid, Field2D::Location::VerticalFace);
    n(0,1) = 1.0;
    for (int i = 0; i < grid.nR(); ++i) {
        vz(i,0) = 7.0; // lower boundary is reflecting regardless of supplied velocity
        vz(i,1) = -0.2;
        vz(i,2) = 0.2;
        vz(i,3) = -0.2; // vacuum inflow at upper boundary
    }
    const double initial = mass(n, grid);
    advection_upwind(n, temp, vr, vz, 0.1, grid);
    check(n(0,0) > 0.0 && n(0,2) > 0.0); // both velocity signs are honored
    check(std::abs(mass(n, grid) - initial) < 1e-13 * initial);
    check(vz(0,0) == 7.0);
    // A single upper-face outflow removes exactly h * area * v * n.
    for (double& value : n.data) value = 0.0;
    n(0,2) = 1.0;
    for (double& value : vz.data) value = 0.0;
    vz(0,3) = 0.2;
    const double open_initial = mass(n, grid);
    advection_upwind(n, temp, vr, vz, 0.1, grid);
    check(std::abs(mass(n, grid) -
                   (open_initial - 0.1 * 0.2 * grid.vertical_face_area(0))) < 1e-13);
    bool rejected = false;
    try { advection_upwind(n, n, vr, vz, 0.1, grid); }
    catch (const std::invalid_argument&) { rejected = true; }
    check(rejected);
}
} // namespace

int main()
{
    closed_grid_test();
    directional_and_boundary_test();
}
