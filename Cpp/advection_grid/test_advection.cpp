#include "CR_advection.hpp"
#include "initialization.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void check(bool ok, const std::string& message)
{
    if (!ok) throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance, const std::string& message)
{
    check(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}
template<class F> void rejects(F action)
{
    bool rejected = false;
    try { action(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Invalid input was accepted");
}
double mass(const Field2D& n, const Grid2D& grid)
{
    double sum = 0.0;
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) sum += n(i,j) * grid.volume(i,j);
    return sum;
}
void test_zero_and_validation()
{
    const Grid2D grid(AxisGrid::linear(3, 0., 1.), AxisGrid::linear(5, 0., 1.));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) n(i,j) = 1. + i + j;
    const Field2D before = n;
    // Deliberately stale ghost cells must be replaced before reconstruction.
    n(-1, 0) = std::numeric_limits<double>::quiet_NaN();
    advection_TVD(n, temp, vr, vz, 100., grid);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) near(n(i,j), before(i,j), 0., "Zero velocity changed density");
    near(n(-1, 0), n(0, 0), 0., "Inner ghost not filled");
    near(n(-2, -2), n(1, 1), 0., "Corner not filled");
    near(n(grid.nR(), 1), 0., 0., "Outer ghost not zero");
    advection_TVD(n, temp, vr, vz, 0., grid);
    rejects([&] { advection_TVD(n, n, vr, vz, 1., grid); });
    rejects([&] { advection_TVD(n, temp, n, vz, 1., grid); });
    rejects([&] { advection_TVD(n, temp, vr, temp, 1., grid); });
    rejects([&] { advection_TVD(n, temp, vr, vz, -1., grid); });
    rejects([&] { advection_TVD(n, temp, vr, vz, std::numeric_limits<double>::infinity(), grid); });
    Field2D wrong(2, 2);
    rejects([&] { advection_TVD(n, wrong, vr, vz, 1., grid); });
    n(0,0) = -1.;
    rejects([&] { advection_TVD(n, temp, vr, vz, 1., grid); });
    n(0,0) = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { advection_TVD(n, temp, vr, vz, 1., grid); });
    n(0,0) = 1.;
    vz(1,5) = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { advection_TVD(n, temp, vr, vz, 1., grid); });
}
void test_closed_stretched_grid()
{
    const Grid2D grid(AxisGrid::geometric(14, 0., 2., 1.4),
                      AxisGrid::geometric(17, 0., 3., 1.3));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) n(i,j) = (i + j) % 3 == 0 ? 1. : 0.;
    // Nonzero supplied lower velocities must still produce reflecting fluxes.
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) {
            vr(i,j) = 0.2 * std::cos(i + j);
            vz(i,j) = 0.2 * std::sin(i - j + 1.);
        }
    const double initial = mass(n, grid);
    advection_TVD(n, temp, vr, vz, 0.6, grid); // Requires many CFL substeps.
    near(mass(n, grid), initial, initial * 2e-13, "Closed-grid mass conservation");
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j)
            check(std::isfinite(n(i,j)) && n(i,j) >= 0., "Negative density on stretched grid");
}
void test_open_boundary()
{
    const Grid2D grid(AxisGrid::linear(3, 0., 2.), AxisGrid::linear(30, 0., 3.));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i) {
        for (int j = 0; j < grid.nz(); ++j) n(i,j) = 1.;
        for (int j = 0; j <= grid.nz(); ++j) vz(i,j) = 0.5;
    }
    const double initial = mass(n, grid);
    advection_TVD(n, temp, vr, vz, 0.05, grid);
    const double area = grid.vertical_face_area(0) + grid.vertical_face_area(1) + grid.vertical_face_area(2);
    near(mass(n, grid), initial - 0.05 * 0.5 * area, 1e-12, "Boundary flux mass balance");
    // A vacuum cannot inject mass, even if the outer velocity points inward.
    std::fill(n.data.begin(), n.data.end(), 0.);
    for (double& v : vz.data) v = -0.5;
    advection_TVD(n, temp, vr, vz, 0.2, grid);
    near(mass(n, grid), 0., 0., "Vacuum inflow injected mass");
}

double temporal_error(int steps)
{
    const Grid2D grid(AxisGrid::linear(6, 0., 1.), AxisGrid::linear(7, 0., 1.));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i <= grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) vr(i,j) = 0.2 * grid.R().face(i);
    for (int i = 0; i < grid.nR(); ++i) {
        for (int j = 0; j <= grid.nz(); ++j) vz(i,j) = 0.2 * grid.z().face(j);
        for (int j = 0; j < grid.nz(); ++j) n(i,j) = 1.;
    }
    for (int step = 0; step < steps; ++step) advection_TVD(n, temp, vr, vz, 0.1 / steps, grid);
    // div(v) = 3*0.2: uniform density decays exponentially.
    const double exact = std::exp(-0.06);
    double error = 0.;
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) error = std::max(error, std::abs(n(i,j) - exact));
    return error;
}

double gaussian_average(double lower, double upper, double shift)
{
    return std::sqrt(std::acos(-1.)) / 20. *
        (std::erf(10. * (upper - shift)) - std::erf(10. * (lower - shift))) / (upper - lower);
}
double spatial_error(int cells, double velocity, bool stretched)
{
    const AxisGrid z = stretched ? AxisGrid::geometric(cells, -1., 1., std::exp(1. / cells))
                                : AxisGrid::linear(cells, -1., 1.);
    const Grid2D grid(AxisGrid::linear(2, 0., 1.), z);
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i) {
        for (int j = 0; j <= grid.nz(); ++j) vz(i,j) = velocity;
        for (int j = 0; j < grid.nz(); ++j)
            n(i,j) = gaussian_average(z.face(j), z.face(j + 1), 0.);
    }
    advection_TVD(n, temp, vr, vz, 0.2, grid);
    double error = 0.;
    for (int j = 0; j < grid.nz(); ++j) {
        check(n(0,j) >= 0., "Gaussian became negative");
        error += z.width(j) * std::abs(n(0,j) - gaussian_average(z.face(j), z.face(j + 1), velocity * 0.2));
    }
    return error;
}
double radial_error(int cells, double velocity)
{
    const Grid2D grid(AxisGrid::geometric(cells, 0., 4., std::exp(1. / cells)),
                      AxisGrid::linear(2, 0., 1.));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    auto primitive = [](double R) {
        return std::sqrt(std::acos(-1.)) / 10. * std::erf(10. * (R - 2.))
             - std::exp(-100. * (R - 2.) * (R - 2.)) / 200.;
    };
    auto average = [&](int i, double time) {
        const double a = grid.R().face(i), b = grid.R().face(i + 1);
        return (primitive(b - velocity * time) - primitive(a - velocity * time))
             / ((b - a) * (a + b) * 0.5);
    };
    for (int i = 0; i <= cells; ++i)
        for (int j = 0; j < 2; ++j) vr(i,j) = velocity;
    for (int i = 0; i < cells; ++i)
        for (int j = 0; j < 2; ++j) n(i,j) = std::max(0., average(i, 0.));
    advection_TVD(n, temp, vr, vz, 0.2, grid);
    double error = 0.;
    for (int i = 0; i < cells; ++i) {
        check(n(i,0) >= 0., "Radial Gaussian became negative");
        error += grid.volume(i,0) * std::abs(n(i,0) - average(i, 0.2));
    }
    return error;
}
void test_substeps()
{
    const Grid2D grid(AxisGrid::linear(2, 0., 1.), AxisGrid::linear(16, 0., 1.));
    Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 16; ++j) n(i,j) = j > 3 && j < 9 ? 1. : 0.;
        for (int j = 0; j <= 16; ++j) vz(i,j) = 1.;
    }
    Field2D reference = n;
    advection_TVD(n, temp, vr, vz, 1., grid);
    // max outgoing rate = 16; ceil(16 / 0.45) = 36 substeps.
    for (int k = 0; k < 36; ++k) advection_TVD(reference, temp, vr, vz, 1. / 36., grid);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 16; ++j)
            near(n(i,j), reference(i,j), 2e-14, "CFL substeps did not cover requested interval");
}
void test_tiny_density_and_tvd()
{
    const Grid2D grid(AxisGrid::linear(2, 0., 1.), AxisGrid::geometric(80, 0., 1., 1.02));
    Field2D n(grid), tiny(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 80; ++j) {
            n(i,j) = (j >= 25 && j < 45) ? 1. : 0.;
            tiny(i,j) = n(i,j) * 1e-100;
        }
        for (int j = 0; j <= 80; ++j) vz(i,j) = 0.5;
    }
    advection_TVD(n, temp, vr, vz, 0.1, grid);
    advection_TVD(tiny, temp, vr, vz, 0.1, grid);
    double variation = std::abs(n(0,0));
    for (int j = 0; j < 80; ++j) {
        check(n(0,j) >= 0. && n(0,j) <= 1. + 1e-14, "Top-hat overshoot");
        near(tiny(0,j) / 1e-100, n(0,j), 2e-14, "Density scale dependence");
        if (j) variation += std::abs(n(0,j) - n(0,j-1));
    }
    variation += std::abs(n(0,79));
    check(variation <= 2. + 1e-13, "1D total variation increased");
}

void test_long_transport()
{
    // Same Gaussian width, grid spacing, velocity and elapsed time as main.cpp.
    const Grid2D grid(AxisGrid::linear(2, 0., 1.), AxisGrid::linear(100, 0., 100.));
    for (double velocity : {-0.1, 0.1}) {
        Field2D n(grid), temp(grid), vz(grid, Field2D::Location::VerticalFace);
        const double start = velocity > 0.0 ? 30.0 : 70.0;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 100; ++j)
                n(i,j) = std::exp(-std::pow((grid.z().center(j) - start) / 3.0, 2));
            for (int j = 0; j <= 100; ++j) vz(i,j) = velocity;
        }
        const Field2D before_velocity = vz;
        const double initial_mass = mass(n, grid);
        for (int step = 0; step < 500; ++step) advection_TVD_Z(n, temp, vz, 1., grid);
        double error = 0., peak = 0., sum = 0., mean = 0., variance = 0.;
        for (int j = 0; j < 100; ++j) {
            const double z = grid.z().center(j);
            check(n(0,j) >= 0., "Long Gaussian transport became negative");
            error += std::abs(n(0,j) - std::exp(-std::pow((z - start - velocity * 500.) / 3., 2)));
            peak = std::max(peak, n(0,j));
            sum += n(0,j);
            mean += z * n(0,j);
        }
        mean /= sum;
        for (int j = 0; j < 100; ++j)
            variance += n(0,j) * std::pow(grid.z().center(j) - mean, 2) / sum;
        std::cout << "Long Gaussian (v=" << velocity << "): L1=" << error
                  << ", peak=" << peak << ", variance=" << variance << '\n';
        // Original implementation: L1=1.41792, peak=0.713125, variance=6.53811.
        check(error < 1.2 && peak > 0.76 && variance < 6.4, "Gaussian dissipation regression");
        near(mass(n, grid), initial_mass, initial_mass * 1e-10, "Long Gaussian mass balance");
        check(vz.data == before_velocity.data, "Advection modified supplied velocity");
        rejects([&] { advection_TVD_Z(n, n, vz, 1., grid); });
    }
}

void test_directional_consistency()
{
    const Grid2D grid(AxisGrid::geometric(12, 0., 2., 1.1), AxisGrid::geometric(15, 0., 3., 1.1));
    for (bool radial : {false, true}) {
        Field2D n(grid), temp(grid), vr(grid, Field2D::Location::RadialFace), vz(grid, Field2D::Location::VerticalFace);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j) n(i,j) = 1. + std::sin(i + j);
        for (double& v : (radial ? vr.data : vz.data)) v = -0.2;
        Field2D reference = n;
        advection_TVD(n, temp, vr, vz, 0.2, grid);
        if (radial) advection_TVD_R(reference, temp, vr, 0.2, grid);
        else advection_TVD_Z(reference, temp, vz, 0.2, grid);
        check(n.data == reference.data, "Directional and unsplit updates disagree for one active direction");
    }
}

} // namespace
int main()
{
    test_zero_and_validation();
    test_closed_stretched_grid();
    test_open_boundary();
    test_substeps();
    test_tiny_density_and_tvd();
    test_long_transport();
    test_directional_consistency();
    const double t1 = temporal_error(1), t2 = temporal_error(2), t4 = temporal_error(4);
    std::cout << "Temporal error ratios: " << t1/t2 << ", " << t2/t4 << '\n';
    check(t1/t2 > 3.8 && t2/t4 > 3.8, "Not second order in time");
    for (bool stretched : {false, true})
        for (double velocity : {-0.4, 0.4}) {
            const double coarse = spatial_error(80, velocity, stretched);
            const double fine = spatial_error(160, velocity, stretched);
            std::cout << "Spatial error ratio (stretched=" << stretched << ", v=" << velocity << "): " << coarse/fine << '\n';
            check(coarse/fine > 2.7, "Insufficient spatial convergence");
        }
    for (double velocity : {-0.4, 0.4}) {
        const double ratio = radial_error(160, velocity) / radial_error(320, velocity);
        std::cout << "Radial spatial error ratio (v=" << velocity << "): " << ratio << '\n';
        check(ratio > 2.7, "Insufficient radial convergence");
    }
    std::cout << "All advection tests passed.\n";
}
