#pragma once
#include <vector>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <limits>

class AxisGrid {
public:
    static constexpr int NG = 2; // Number of ghost cells on each side.

    // Uniform cell widths.
    static AxisGrid linear(int n, double lower, double upper);

    // Successive physical cell widths satisfy:
    // width(i + 1) / width(i) = ratio.
    //
    // ratio == 1 gives a uniform grid.
    static AxisGrid geometric(int n, double lower, double upper, double ratio);

    // Arbitrary strictly increasing physical faces.
    // Number of cells is physical_faces.size() - 1.
    static AxisGrid from_faces(std::vector<double> physical_faces);

    int size() const noexcept { return n_; }

    double face(int i) const;
    double center(int i) const;
    double width(int i) const;

    // center(i + 1) - center(i)
    double center_distance(int i) const;

private:
    explicit AxisGrid(std::vector<double> physical_faces);

    int n_;
    std::vector<double> faces_;
    std::vector<double> centers_;
    std::vector<double> widths_;
    std::vector<double> center_distances_;
};

class Grid2D {
    public:
        Grid2D(AxisGrid radial, AxisGrid vertical);

        int nR() const noexcept { return radial_.size(); }
        int nz() const noexcept { return vertical_.size(); }

        const AxisGrid& R() const noexcept { return radial_; }
        const AxisGrid& z() const noexcept { return vertical_; }

        // Full axisymmetric annular-cell volume, including 2*pi.
        double volume(int i, int j) const;

        // Surface areas for conservative flux calculations.
        double radial_face_area(int i, int j) const;
        double vertical_face_area(int i) const;

        // Cached cylindrical volume centroids in physical cells; ghost coordinates
        // reflect physical centroids across the boundary faces.
        // Valid indices: -AxisGrid::NG <= i < nR() + AxisGrid::NG.
        double radial_centroid(int i) const;
        // centroid(i + 1) - centroid(i), including ghosts.
        // Valid indices: -AxisGrid::NG <= i < nR() + AxisGrid::NG - 1.
        double radial_centroid_distance(int i) const;


    private:
        AxisGrid radial_;
        AxisGrid vertical_;
        std::vector<double> radial_centroid_;
        std::vector<double> radial_centroid_distance_;

};

class Field2D
{
public:
    enum class Location { Centroid, RadialFace, VerticalFace };

private:
    Location location_;
    int nR_;
    int nz_;

    std::size_t index(int i, int j) const
    {
        assert(i >= -NG && static_cast<long long>(i) < static_cast<long long>(nR_) + NG);
        assert(j >= -NG && static_cast<long long>(j) < static_cast<long long>(nz_) + NG);
        return (static_cast<std::size_t>(static_cast<long long>(i) + NG)) * (static_cast<std::size_t>(nz_) + 2 * NG)
             + static_cast<std::size_t>(static_cast<long long>(j) + NG);
    }

public:
    static constexpr int NG = AxisGrid::NG;
    std::vector<double> data;

    explicit Field2D(const Grid2D& grid, Location location = Location::Centroid)
        : Field2D(grid.nR(), grid.nz(), location) {}

    // Arguments are grid CELL counts; stored dimensions count physical values.
    Field2D(int nR, int nz, Location location = Location::Centroid)
        : location_(location), nR_(nR), nz_(nz)
    {
        if (nR < NG || nz < NG) {
            throw std::invalid_argument("Field2D requires at least two cells per axis");
        }
        if (nR > std::numeric_limits<int>::max() - 2 * NG - 1 ||
            nz > std::numeric_limits<int>::max() - 2 * NG - 1) {
            throw std::length_error("Field2D dimensions exceed index limits");
        }
        if (location == Location::RadialFace) ++nR_;
        if (location == Location::VerticalFace) ++nz_;
        const std::size_t rows = static_cast<std::size_t>(nR_) + 2 * NG;
        const std::size_t cols = static_cast<std::size_t>(nz_) + 2 * NG;
        if (rows > data.max_size() / cols) {
            throw std::length_error("Field2D dimensions exceed storage limits");
        }
        data.resize(rows * cols);
    }

    // Physical sample counts (excluding ghosts), not necessarily cell counts.
    int nR() const noexcept { return nR_; }
    int nz() const noexcept { return nz_; }

    Location location() const noexcept { return location_; }
    bool matches(const Grid2D& grid, Location expected) const noexcept
    {
        return location_ == expected &&
               nR_ == grid.nR() + (expected == Location::RadialFace) &&
               nz_ == grid.nz() + (expected == Location::VerticalFace);
    }

    double& operator()(int i, int j) { return data[index(i, j)]; }
    const double& operator()(int i, int j) const { return data[index(i, j)]; }
};
