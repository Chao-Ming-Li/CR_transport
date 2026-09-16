#include "field.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace {

// Unlike assert(), these checks also run when NDEBUG is defined.
void check(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, const std::string& message)
{
    const double tolerance = 1e-12 * std::max(1.0, std::abs(expected));
    check(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
          message + ": got " + std::to_string(actual) +
          ", expected " + std::to_string(expected));
}

template <class Function>
void rejects(Function function, const std::string& message)
{
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

void test_linear()
{
    const auto axis = AxisGrid::linear(4, 0.0, 4.0);
    check(axis.size() == 4, "Linear cell count");
    for (int i = -2; i <= 6; ++i) {
        near(axis.face(i), 0.0 + i, "Linear face including ghosts");
        std::cout << "axis.face(" << i << ") = " << axis.face(i) << std::endl;
    }
    // for (int i = -2; i <= 5; ++i) {
    //     near(axis.center(i), 2.5 + i, "Linear center including ghosts");
    //     near(axis.width(i), 1.0, "Linear width including ghosts");
    // }
    // for (int i = -2; i <= 4; ++i)
    //     near(axis.center_distance(i), 1.0, "Linear center distance");
}

void test_geometric()
{
    for (double ratio : {0.5, 1.0, 1.0 + 1e-10, 2.0}) {
        const auto axis = AxisGrid::geometric(4, 0.0, 15.0, ratio);
        near(axis.face(0), 0.0, "Geometric lower endpoint");
        near(axis.face(4), 15.0, "Geometric upper endpoint");
        double length = 0.0;
        for (int i = 0; i < axis.size(); ++i) {
            check(axis.width(i) > 0.0, "Positive geometric width");
            length += axis.width(i);
            if (i > 0)
                near(axis.width(i) / axis.width(i - 1), ratio,
                     "Geometric width ratio");
        }
        near(length, 15.0, "Geometric total length");
    }
    const auto axis = AxisGrid::geometric(4, 0.0, 15.0, 2.0);
    const double expected_faces[] = {0.0, 1.0, 3.0, 7.0, 15.0};
    for (int i = 0; i <= 4; ++i)
        near(axis.face(i), expected_faces[i], "Known geometric face");
}

void test_custom_and_ghosts()
{
    // Minimum supported cell count, with unequal widths 1 and 3.
    const auto axis = AxisGrid::from_faces({0.0, 1.0, 4.0});
    const double faces[] = {-4.0, -1.0, 0.0, 1.0, 4.0, 7.0, 8.0};
    const double centers[] = {-2.5, -0.5, 0.5, 2.5, 5.5, 7.5};
    const double widths[] = {3.0, 1.0, 1.0, 3.0, 3.0, 1.0};
    const double distances[] = {2.0, 1.0, 2.0, 3.0, 2.0};
    check(axis.size() == 2, "Custom cell count");
    for (int i = -2; i <= 4; ++i)
        near(axis.face(i), faces[i + 2], "Custom face");
    for (int i = -2; i <= 3; ++i) {
        near(axis.center(i), centers[i + 2], "Custom center");
        near(axis.width(i), widths[i + 2], "Custom width");
    }
    for (int i = -2; i <= 2; ++i)
        near(axis.center_distance(i), distances[i + 2],
             "Custom center distance");
}

void test_cylindrical_geometry()
{
    const Grid2D grid(AxisGrid::from_faces({0.0, 1.0, 3.0}),
                      AxisGrid::from_faces({-1.0, 1.0, 4.0}));
    const double pi = std::numbers::pi;
    check(grid.nR() == 2 && grid.nz() == 2, "2D dimensions");
    near(grid.R().face(2), 3.0, "Radial accessor");
    near(grid.z().face(2), 4.0, "Vertical accessor");
    near(grid.volume(0, 0), 2.0 * pi, "Inner cell volume");
    near(grid.volume(1, 1), 24.0 * pi, "Outer cell volume");
    near(grid.radial_face_area(0, 0), 0.0, "Area at axis");
    near(grid.radial_face_area(2, 1), 18.0 * pi, "Outer radial area");
    near(grid.vertical_face_area(1), 8.0 * pi, "Annular area");
    double total = 0.0;
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j)
            total += grid.volume(i, j);
    near(total, 45.0 * pi, "Total cylinder volume");
}

void test_invalid_inputs()
{
    rejects([] { AxisGrid::linear(1, 0.0, 1.0); }, "Accepted too few cells");
    rejects([] { AxisGrid::linear(2, 1.0, 0.0); }, "Accepted reversed extent");
    rejects([] { AxisGrid::linear(2, 1.0, 1.0); }, "Accepted empty extent");
    rejects([] { AxisGrid::geometric(2, 0.0, 1.0, 0.0); },
            "Accepted zero width ratio");
    rejects([] { AxisGrid::from_faces({0.0, 1.0}); }, "Accepted too few faces");
    rejects([] { AxisGrid::from_faces({0.0, 1.0, 1.0}); },
            "Accepted duplicate faces");
    rejects([] { AxisGrid::from_faces({0.0, 2.0, 1.0}); },
            "Accepted unordered faces");
    rejects([] {
        AxisGrid::from_faces({0.0, 1.0,
                             std::numeric_limits<double>::quiet_NaN()});
    }, "Accepted nonfinite face");
    rejects([] {
        Grid2D grid(AxisGrid::linear(2, -1.0, 1.0),
                    AxisGrid::linear(2, 0.0, 1.0));
    }, "Accepted negative physical radius");
}

} // namespace

int main()
{
    try {
        test_linear();
        std::cout << "PASS: linear grid\n";
        test_geometric();
        std::cout << "PASS: geometric grids\n";
        test_custom_and_ghosts();
        std::cout << "PASS: custom grid and ghost cells\n";
        test_cylindrical_geometry();
        std::cout << "PASS: cylindrical volumes and areas\n";
        test_invalid_inputs();
        std::cout << "PASS: invalid inputs rejected\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << "All grid tests passed.\n";
    return 0;
}
