#include "field.hpp"
#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace {
    // Validate that the grid extent is valid for the given number of cells.
    void validate_extent(int n, double lower, double upper)
    {
        if (n < AxisGrid::NG) {
            throw std::invalid_argument(
                "AxisGrid requires at least two physical cells");
        }

        if (!std::isfinite(lower) ||
            !std::isfinite(upper) ||
            !std::isfinite(upper - lower) ||
            upper <= lower) {
            throw std::invalid_argument("Invalid grid extent");
        }
    }

} // namespace

AxisGrid AxisGrid::linear(int n, double lower, double upper)
{
    validate_extent(n, lower, upper);

    std::vector<double> faces(static_cast<std::size_t>(n) + 1);

    for (int i = 0; i <= n; ++i) {
        faces[i] = lower + (upper - lower) * (static_cast<double>(i) / n);
    }

    faces.front() = lower;
    faces.back() = upper;

    return AxisGrid(std::move(faces));
}

AxisGrid AxisGrid::geometric(int n, double lower, double upper, double ratio)
{
    validate_extent(n, lower, upper);

    if (!std::isfinite(ratio) || ratio <= 0.0) {
        throw std::invalid_argument(
            "Geometric width ratio must be positive and finite");
    }

    if (ratio == 1.0) {
        return linear(n, lower, upper);
    }

    const double a = std::log(ratio);
    std::vector<double> faces(static_cast<std::size_t>(n) + 1);

    for (int i = 0; i <= n; ++i) {
        double fraction;

        if (a > 0.0) {
            // Equivalent to (ratio^i - 1)/(ratio^n - 1),
            // arranged to avoid overflowing positive exponentials.
            fraction =
                std::exp((static_cast<double>(i) - n) * a) *
                (-std::expm1(-static_cast<double>(i) * a)) /
                (-std::expm1(-static_cast<double>(n) * a));
        } else {
            fraction =
                std::expm1(static_cast<double>(i) * a) /
                std::expm1(static_cast<double>(n) * a);
        }

        faces[i] = lower + (upper - lower) * fraction;
    }

    faces.front() = lower;
    faces.back() = upper;

    // Constructor rejects grids whose cells collapse under
    // floating-point rounding because stretching is too extreme.
    return AxisGrid(std::move(faces));
}

AxisGrid AxisGrid::from_faces(std::vector<double> physical_faces)
{
    return AxisGrid(std::move(physical_faces));
}

AxisGrid::AxisGrid(std::vector<double> physical_faces)
    : n_(0)
{
    if (physical_faces.size() < static_cast<std::size_t>(NG + 1)) {
        throw std::invalid_argument("Too few grid faces");
    }

    if (physical_faces.size() > 
        static_cast<std::size_t>(std::numeric_limits<int>::max() - 2 * NG)) {
        throw std::invalid_argument("Too many grid faces");
    }

    n_ = static_cast<int>(physical_faces.size()) - 1;

    for (int i = 0; i <= n_; ++i) {
        if (!std::isfinite(physical_faces[i])) {
            throw std::invalid_argument("Nonfinite grid face");
        }

        if (i > 0 && physical_faces[i] <= physical_faces[i - 1]) {
            throw std::invalid_argument(
                "Grid faces must be strictly increasing");
        }
    }

    faces_.resize(n_ + 1 + 2 * NG);
    centers_.resize(n_ + 2 * NG);
    widths_.resize(n_ + 2 * NG);
    center_distances_.resize(n_ - 1 + 2 * NG);

    for (int i = 0; i <= n_; ++i) {
        faces_[i + NG] = physical_faces[i];
    }

    // Reflect geometry across each domain boundary.
    for (int g = 1; g <= NG; ++g) {
        faces_[NG - g] =
            physical_faces.front() -
            (physical_faces[g] - physical_faces.front());

        faces_[n_ + NG + g] =
            physical_faces.back() +
            (physical_faces.back() - physical_faces[n_ - g]);
    }

    for (int i = -NG; i <= n_ + NG; ++i) {
        if (!std::isfinite(face(i))) {
            throw std::invalid_argument("Nonfinite ghost face");
        }
    }

    for (int i = -NG; i < n_ + NG; ++i) {
        const double dx = face(i + 1) - face(i);

        if (!std::isfinite(dx) || dx <= 0.0) {
            throw std::invalid_argument("Invalid cell width");
        }

        widths_[i + NG] = dx;
        centers_[i + NG] = face(i) + 0.5 * dx;
    }

    for (int i = -NG; i < n_ + NG - 1; ++i) {
        const double distance = center(i + 1) - center(i);

        if (!std::isfinite(distance) || distance <= 0.0) {
            throw std::invalid_argument("Invalid center spacing");
        }

        center_distances_[i + NG] = distance;
    }
}

double AxisGrid::face(int i) const
{
    assert(i >= -NG && i <= n_ + NG);
    return faces_[i + NG];
}

double AxisGrid::center(int i) const
{
    assert(i >= -NG && i < n_ + NG);
    return centers_[i + NG];
}

double AxisGrid::width(int i) const
{
    assert(i >= -NG && i < n_ + NG);
    return widths_[i + NG];
}

double AxisGrid::center_distance(int i) const
{
    assert(i >= -NG && i < n_ + NG - 1);
    return center_distances_[i + NG];
}

Grid2D::Grid2D(AxisGrid radial, AxisGrid vertical)
    : radial_(std::move(radial)),
      vertical_(std::move(vertical))
{
    if (radial_.face(0) < 0.0) {
        throw std::invalid_argument(
            "Physical radial coordinates must be nonnegative");
    }
}

double Grid2D::volume(int i, int j) const
{
    assert(i >= 0 && i < nR());
    assert(j >= 0 && j < nz());

    return vertical_face_area(i) * vertical_.width(j);
}

double Grid2D::radial_face_area(int i, int j) const
{
    assert(i >= 0 && i <= nR());
    assert(j >= 0 && j < nz());

    return 2.0 * std::numbers::pi *
           radial_.face(i) * vertical_.width(j);
}

double Grid2D::vertical_face_area(int i) const
{
    assert(i >= 0 && i < nR());

    const double left = radial_.face(i);
    const double right = radial_.face(i + 1);

    // Factored form of pi * (right^2 - left^2).
    return std::numbers::pi * (right - left) * (right + left);
}