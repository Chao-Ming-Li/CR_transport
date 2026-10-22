#include "field.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>


AxisGrid AxisGrid::linear(int n, double lower, double upper) {

    std::vector<double> faces(static_cast<std::size_t>(n) + 1);

    for (int i = 0; i <= n; ++i) {
        faces[i] = lower + (upper - lower) * (static_cast<double>(i) / n);
    }

    faces.front() = lower;
    faces.back() = upper;

    return AxisGrid(std::move(faces));
}

AxisGrid AxisGrid::geometric(int n, double lower, double upper, double ratio) {


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
            fraction = std::exp((static_cast<double>(i) - n) * a) *
                       (-std::expm1(-static_cast<double>(i) * a)) /
                       (-std::expm1(-static_cast<double>(n) * a));
        } else {
            fraction =
                std::expm1(static_cast<double>(i) * a) / std::expm1(static_cast<double>(n) * a);
        }

        faces[i] = lower + (upper - lower) * fraction;
    }

    faces.front() = lower;
    faces.back() = upper;

    // Constructor rejects grids whose cells collapse under
    // floating-point rounding because stretching is too extreme.
    return AxisGrid(std::move(faces));
}

AxisGrid AxisGrid::geometric_capped(double lower, double upper, double first_width,
                                    double max_width, double ratio) {
    std::vector<double> faces{lower};
    double width = first_width;
    // Clip only the final cell so the physical boundary is exactly upper.
    while (faces.back() < upper) {
        const double next = std::min(upper, faces.back() + width);
        faces.push_back(next);
        // Cap before multiplying to avoid overflow for a large ratio.
        width = width >= max_width / ratio ? max_width : width * ratio;
    }
    return AxisGrid(std::move(faces));
}

AxisGrid AxisGrid::from_faces(std::vector<double> physical_faces) {
    return AxisGrid(std::move(physical_faces));
}

AxisGrid::AxisGrid(std::vector<double> physical_faces) : n_(0) {
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
            throw std::invalid_argument("Grid faces must be strictly increasing");
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
        faces_[NG - g] = physical_faces.front() - (physical_faces[g] - physical_faces.front());

        faces_[n_ + NG + g] =
            physical_faces.back() + (physical_faces.back() - physical_faces[n_ - g]);
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

double AxisGrid::face(int i) const {
    return faces_[i + NG];
}

double AxisGrid::center(int i) const {
    return centers_[i + NG];
}

double AxisGrid::width(int i) const {
    return widths_[i + NG];
}

double AxisGrid::center_distance(int i) const {
    return center_distances_[i + NG];
}

SpatialGrid::SpatialGrid(AxisGrid radial, AxisGrid vertical)
    : radial_(std::move(radial)), vertical_(std::move(vertical)) {
    if (radial_.face(0) < 0.0) {
        throw std::invalid_argument("Physical radial coordinates must be nonnegative");
    }
    constexpr int ng = AxisGrid::NG;
    radial_centroid_.resize(nR() + 2 * ng);
    radial_centroid_distance_.resize(nR() - 1 + 2 * ng);

    for (int i = 0; i < nR(); ++i) {
        const double center = radial_.center(i);
        const double width = radial_.width(i);
        // Cylindrical volume weighting; valid for nonuniform cells too.
        radial_centroid_[i + ng] = center + width * (width / center) / 12.0;
    }
    // Extend the reflected ghost geometry used by AxisGrid. These are
    // reflected coordinates, not centroids of additional physical annuli.
    // This remains well defined even if an inner ghost cell crosses R=0.
    const double lower = radial_.face(0), upper = radial_.face(nR());
    for (int g = 1; g <= ng; ++g) {
        radial_centroid_[ng - g] = lower - (radial_centroid_[ng + g - 1] - lower);
        radial_centroid_[ng + nR() + g - 1] = upper + (upper - radial_centroid_[ng + nR() - g]);
    }
    for (int i = -ng; i < nR() + ng; ++i) {
        const double centroid = radial_centroid_[i + ng];
        if (!std::isfinite(centroid) || centroid <= radial_.face(i) ||
            centroid >= radial_.face(i + 1)) {
            throw std::invalid_argument("Invalid radial centroid");
        }
    }
    for (int i = -ng; i < nR() + ng - 1; ++i) {
        const double distance = radial_centroid_[i + ng + 1] - radial_centroid_[i + ng];
        if (!std::isfinite(distance) || distance <= 0.0) {
            throw std::invalid_argument("Invalid radial centroid spacing");
        }
        radial_centroid_distance_[i + ng] = distance;
    }
}

double SpatialGrid::volume(int i, int j) const {

    return vertical_face_area(i) * vertical_.width(j);
}

double SpatialGrid::radial_face_area(int i, int j) const {

    return 2.0 * std::numbers::pi * radial_.face(i) * vertical_.width(j);
}

double SpatialGrid::vertical_face_area(int i) const {

    const double left = radial_.face(i);
    const double right = radial_.face(i + 1);

    // Factored form of pi * (right^2 - left^2).
    return std::numbers::pi * (right - left) * (right + left);
}

double SpatialGrid::radial_centroid(int i) const {
    return radial_centroid_[i + AxisGrid::NG];
}

double SpatialGrid::radial_centroid_distance(int i) const {
    return radial_centroid_distance_[i + AxisGrid::NG];
}

RgGrid::RgGrid(std::vector<double> centers, std::vector<double> faces)
    : centers_(std::move(centers)), faces_(std::move(faces)) {}

RgGrid RgGrid::logarithmic(int n_faces, double lower, double upper) {
    if (n_faces < 2 || !std::isfinite(lower) || !std::isfinite(upper) ||
        lower <= 0.0 || upper <= lower)
        throw std::invalid_argument("Invalid logarithmic rigidity grid");

    const double log_lower = std::log(lower);
    const double dlog = (std::log(upper) - log_lower) / (n_faces - 1);
    std::vector<double> faces(n_faces), centers(n_faces - 1);
    for (int k = 0; k < n_faces; ++k)
        faces[k] = std::exp(log_lower + k * dlog);
    faces.front() = lower;
    faces.back() = upper;
    for (int k = 0; k < n_faces - 1; ++k) {
        centers[k] = std::exp(0.5 * (std::log(faces[k]) + std::log(faces[k + 1])));
        if (!(faces[k] < centers[k] && centers[k] < faces[k + 1]))
            throw std::invalid_argument("Rigidity cell collapsed under rounding");
    }
    return RgGrid(std::move(centers), std::move(faces));
}

double RgGrid::center(int k) const {
    return centers_[k];
}

double RgGrid::face(int k) const {
    return faces_[k];
}

double RgGrid::width(int k) const {
    return faces_[k + 1] - faces_[k];
}

double RgGrid::center_distance(int k) const {
    return centers_[k + 1] - centers_[k];
}

double RgGrid::log_width(int k) const {
    return std::log(faces_[k + 1]) - std::log(faces_[k]);
}

Field3D::Field3D(const SpatialGrid &grid, const RgGrid &rigidity)
    : nR_(grid.nR()), nz_(grid.nz()), nrigidity_(rigidity.size()), spatial_size_(0) {
    const std::size_t rows = static_cast<std::size_t>(nR_) + 2 * NG;
    const std::size_t cols = static_cast<std::size_t>(nz_) + 2 * NG;
    if (rows > data_.max_size() / cols)
        throw std::length_error("Spatial plane exceeds storage limits");
    spatial_size_ = rows * cols;
    if (static_cast<std::size_t>(nrigidity_) > data_.max_size() / spatial_size_)
        throw std::length_error("Field3D exceeds storage limits");
    data_.resize(static_cast<std::size_t>(nrigidity_) * spatial_size_, 0.0);
}

std::size_t Field3D::index(int k, int i, int j) const {
    return static_cast<std::size_t>(k) * spatial_size_ +
           static_cast<std::size_t>(static_cast<long long>(i) + NG) *
               (static_cast<std::size_t>(nz_) + 2 * NG) +
           static_cast<std::size_t>(static_cast<long long>(j) + NG);
}

std::span<double> Field3D::spatial_plane(int k) {
    return std::span<double>(data_).subspan(static_cast<std::size_t>(k) * spatial_size_, spatial_size_);
}

std::span<const double> Field3D::spatial_plane(int k) const {
    return std::span<const double>(data_).subspan(static_cast<std::size_t>(k) * spatial_size_, spatial_size_);
}
