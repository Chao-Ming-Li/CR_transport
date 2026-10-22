#include "field.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace {
    void check(bool condition) {
        if (!condition)
            throw std::runtime_error("Grid/field check failed");
    }

    bool close(double a, double b) {
        return std::abs(a - b) <= 1e-12 * std::max({1.0, std::abs(a), std::abs(b)});
    }
}

int main() {
    const Grid2D grid(AxisGrid::geometric_capped(0.0, 200.0, 0.25, 1.0, 1.03),
                      AxisGrid::linear(3, 0.0, 6.0));
    check(grid.R().face(grid.nR()) == 200.0);
    check(grid.R().width(grid.nR() - 1) <= 1.0);
    check(close(grid.radial_centroid(0), 1.0 / 6.0));
    check(close(grid.radial_centroid(-1), -grid.radial_centroid(0)));
    check(close(grid.z().center(-1), -grid.z().center(0)));
    double volume = 0.0;
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j)
            volume += grid.volume(i, j);
    check(close(volume, std::numbers::pi * 200.0 * 200.0 * 6.0));

    const auto rigidity = RigidityGrid::logarithmic(61, 1.0, 1e6);
    check(rigidity.size() == 60);
    check(rigidity.face(0) == 1.0 && rigidity.face(60) == 1e6);
    double total_width = 0.0;
    for (int k = 0; k < rigidity.size(); ++k) {
        check(rigidity.face(k) < rigidity.center(k));
        check(close(rigidity.center(k), std::sqrt(rigidity.face(k) * rigidity.face(k + 1))));
        check(rigidity.center(k) < rigidity.face(k + 1));
        check(close(rigidity.log_width(k), std::log(10.0) / 10.0));
        total_width += rigidity.width(k);
    }
    check(close(total_width, rigidity.face(60) - rigidity.face(0)));

    Field3D density(grid, rigidity);
    check(density.data().size() == 60 * density.spatial_size());
    for (int k = 0; k < density.nrigidity(); ++k)
        for (int i = -2; i < grid.nR() + 2; ++i)
            for (int j = -2; j < grid.nz() + 2; ++j)
                density(k, i, j) = k * 1000000.0 + (i + 2) * 100.0 + j + 2;
    const Field3D &read_only = density;
    for (int k = 0; k < density.nrigidity(); ++k) {
        check(&density(k, -2, -2) == density.spatial_plane(k).data());
        check(close(read_only(k, grid.nR() + 1, grid.nz() + 1),
                    k * 1000000.0 + (grid.nR() + 3) * 100.0 + grid.nz() + 3));
        check(read_only.spatial_plane(k).size() == density.spatial_size());
    }
    check(&density(0, 0, 1) == &density(0, 0, 0) + 1);
    density.spatial_plane(1)[0] = -5.0;
    check(density(1, -2, -2) == -5.0);
    Field2D radial_velocity(grid, Field2D::Location::RadialFace);
    check(radial_velocity.nR() == grid.nR() + 1);
    std::cout << "Spatial geometry, rigidity bins and Field3D layout passed.\n";
}
