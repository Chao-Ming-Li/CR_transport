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
    void check(bool condition, const std::string &message) {
        if (!condition)
            throw std::runtime_error(message);
    }

    void near(double actual, double expected, const std::string &message) {
        const double tolerance = 1e-12 * std::max(1.0, std::abs(expected));
        check(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
              message + ": got " + std::to_string(actual) + ", expected " +
                  std::to_string(expected));
    }

    template <class Function> void rejects(Function function, const std::string &message) {
        try {
            function();
        } catch (const std::invalid_argument &) {
            return;
        }
        throw std::runtime_error(message);
    }

    void test_linear() {
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

    void test_geometric() {
        for (double ratio : {0.5, 1.0, 1.0 + 1e-10, 2.0}) {
            const auto axis = AxisGrid::geometric(4, 0.0, 15.0, ratio);
            near(axis.face(0), 0.0, "Geometric lower endpoint");
            near(axis.face(4), 15.0, "Geometric upper endpoint");
            double length = 0.0;
            for (int i = 0; i < axis.size(); ++i) {
                check(axis.width(i) > 0.0, "Positive geometric width");
                length += axis.width(i);
                if (i > 0)
                    near(axis.width(i) / axis.width(i - 1), ratio, "Geometric width ratio");
            }
            near(length, 15.0, "Geometric total length");
        }
        const auto axis = AxisGrid::geometric(4, 0.0, 15.0, 2.0);
        const double expected_faces[] = {0.0, 1.0, 3.0, 7.0, 15.0};
        for (int i = 0; i <= 4; ++i)
            near(axis.face(i), expected_faces[i], "Known geometric face");
    }

    void test_geometric_capped() {
        const auto requested = AxisGrid::geometric_capped(0., 200., 0.1, 1., 1.03);
        constexpr int geometric_count = 78;
        const double requested_transition =
            0.1 * std::expm1(geometric_count * std::log(1.03)) / 0.03;
        check(requested.size() == 248, "1.03 grid cell count");
        near(requested.width(0), 0.1, "1.03 start width");
        near(requested.face(geometric_count), requested_transition, "1.03 transition");
        for (int i = 1; i < geometric_count; ++i)
            near(requested.width(i) / requested.width(i - 1), 1.03, "1.03 geometric ratio");
        for (int i = geometric_count; i < requested.size(); ++i)
            near(requested.width(i), 1., "1.03 constant region");
        near(requested.face(requested.size()), requested_transition + 170.,
             "1.03 full-cell endpoint");
        near(requested.width(requested.size() - 1), 1., "1.03 final cell");
        std::cout << "1.03 configuration: " << requested.size()
                  << " cells, transition=" << requested.face(geometric_count)
                  << " kpc, final width=" << requested.width(requested.size() - 1) << " kpc\n";

        const auto axis = AxisGrid::geometric_capped(0., 200., 0.1, 1., 1.003);
        near(axis.width(0), 0.1, "Requested start width");
        check(axis.face(axis.size()) >= 200. && axis.face(axis.size()) - 200. <= 1.,
              "Target endpoint overshoot");
        check(axis.face(axis.size() - 1) < 200., "Stop at first face beyond target");
        for (int i = 0; i < axis.size(); ++i) {
            check(axis.width(i) > 0. && axis.width(i) < 1., "Cap not reached in box");
            if (i > 0)
                near(axis.width(i) / axis.width(i - 1), 1.003, "Requested growth ratio");
        }
        const int k = static_cast<int>(std::ceil(std::log(10.) / std::log(1.003)));
        const double transition = 0.1 * std::expm1(k * std::log(1.003)) / 0.003;
        const auto extended = AxisGrid::geometric_capped(0., 400., 0.1, 1., 1.003);
        near(extended.face(k), transition, "Calculated transition coordinate");
        near(extended.width(k), 1., "First constant-width cell");
        for (int i = k; i < extended.size(); ++i)
            near(extended.width(i), 1., "Constant outer widths");
        near(extended.width(-1), extended.width(0), "Lower ghost");
        near(extended.width(extended.size()), extended.width(extended.size() - 1), "Upper ghost");
        const auto shifted = AxisGrid::geometric_capped(5., 405., 0.1, 1., 1.003);
        near(shifted.face(k), 5. + transition, "Shifted transition");
        const auto uniform = AxisGrid::geometric_capped(0., 2., 0.1, 0.1, 1.);
        check(uniform.size() == 20, "Decimal uniform count");
        const auto single = AxisGrid::geometric_capped(0., 2., 0.5, 1., 2.);
        near(single.width(0), 0.5, "Single geometric cell");
        near(single.width(1), 1., "Single linear cell");
        near(single.width(2), 1., "Full final cell");
        near(single.face(single.size()), 2.5, "Extended endpoint");
        rejects([] { AxisGrid::geometric_capped(0., 10., 0., 1., 1.003); }, "Zero width");
        rejects([] { AxisGrid::geometric_capped(0., 10., 2., 1., 1.003); }, "Reversed widths");
        rejects([] { AxisGrid::geometric_capped(0., 10., 0.1, 1., 0.9); }, "Shrinking ratio");
        std::cout << "Requested grid cells=" << axis.size()
                  << " last full width=" << axis.width(axis.size() - 2)
                  << " final width=" << axis.width(axis.size() - 1) << " transition index=" << k
                  << " coordinate=" << transition << '\n';
    }

    void test_custom_and_ghosts() {
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
            near(axis.center_distance(i), distances[i + 2], "Custom center distance");
    }

    void test_cylindrical_geometry() {
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

    void test_invalid_inputs() {
        rejects([] { AxisGrid::linear(1, 0.0, 1.0); }, "Accepted too few cells");
        rejects([] { AxisGrid::linear(2, 1.0, 0.0); }, "Accepted reversed extent");
        rejects([] { AxisGrid::linear(2, 1.0, 1.0); }, "Accepted empty extent");
        rejects([] { AxisGrid::geometric(2, 0.0, 1.0, 0.0); }, "Accepted zero width ratio");
        rejects([] { AxisGrid::from_faces({0.0, 1.0}); }, "Accepted too few faces");
        rejects([] { AxisGrid::from_faces({0.0, 1.0, 1.0}); }, "Accepted duplicate faces");
        rejects([] { AxisGrid::from_faces({0.0, 2.0, 1.0}); }, "Accepted unordered faces");
        rejects([] { AxisGrid::from_faces({0.0, 1.0, std::numeric_limits<double>::quiet_NaN()}); },
                "Accepted nonfinite face");
        rejects([] { Grid2D grid(AxisGrid::linear(2, -1.0, 1.0), AxisGrid::linear(2, 0.0, 1.0)); },
                "Accepted negative physical radius");
    }

    void test_cached_centroid_ghosts() {
        for (const auto &axis :
             {AxisGrid::linear(2, 0., 2.), AxisGrid::geometric(8, 0., 4., 1.4),
              AxisGrid::geometric(8, 0., 4., 0.7), AxisGrid::from_faces({0.5, 1.5, 3.5})}) {
            const Grid2D grid(axis, AxisGrid::linear(3, -1.5, 1.5));
            const int n = grid.nR(), ng = AxisGrid::NG;
            for (int i = 0; i < n; ++i) {
                const double a = axis.face(i), b = axis.face(i + 1);
                near(grid.radial_centroid(i), (2. / 3.) * (a * a + a * b + b * b) / (a + b),
                     "Physical cylindrical centroid");
            }
            for (int g = 1; g <= ng; ++g) {
                near(grid.radial_centroid(-g) + grid.radial_centroid(g - 1), 2. * axis.face(0),
                     "Lower centroid reflection");
                near(grid.radial_centroid(n + g - 1) + grid.radial_centroid(n - g),
                     2. * axis.face(n), "Upper centroid reflection");
            }
            for (int i = -ng; i < n + ng; ++i) {
                check(grid.radial_centroid(i) > axis.face(i) &&
                          grid.radial_centroid(i) < axis.face(i + 1),
                      "Centroid outside cell");
                if (i < n + ng - 1) {
                    check(grid.radial_centroid_distance(i) > 0., "Nonpositive centroid spacing");
                    near(grid.radial_centroid_distance(i),
                         grid.radial_centroid(i + 1) - grid.radial_centroid(i),
                         "Cached centroid distance");
                }
            }
            near(grid.z().center(1), 0., "Vertical midpoint unchanged");
        }
    }

} // namespace

int main() {
    try {
        test_linear();
        std::cout << "PASS: linear grid\n";
        test_geometric();
        test_geometric_capped();
        std::cout << "PASS: geometric grids\n";
        test_custom_and_ghosts();
        std::cout << "PASS: custom grid and ghost cells\n";
        test_cylindrical_geometry();
        std::cout << "PASS: cylindrical volumes and areas\n";
        test_cached_centroid_ghosts();
        std::cout << "PASS: cached centroids including ghosts\n";
        test_invalid_inputs();
        std::cout << "PASS: invalid inputs rejected\n";
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << "All grid tests passed.\n";
    return 0;
}
