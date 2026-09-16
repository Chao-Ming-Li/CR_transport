#include "field.hpp"
#include "CR_advection.hpp"  
#include "initialization.hpp"
#include <algorithm>
#include <limits>
#include <fstream>
#include <iostream>
#include <cmath>   

namespace {
void validate_legacy_field(const Field2D& field, const Grid2D& grid,
                         Field2D::Location location = Field2D::Location::Center)
{
    if (!field.matches(grid, location)) {
        throw std::invalid_argument("Field dimensions and location must match the grid");
    }
}
} // namespace

void apply_boundary_conditions(Field2D& ndis, const Grid2D& grid) {
    validate_legacy_field(ndis, grid);
    for (int i = -Field2D::NG; i < ndis.nR() + Field2D::NG; i++) {
        ndis(i, ndis.nz()) = 0.0;     // ghost cell in upper absorbing boundary in z
        ndis(i, ndis.nz() + 1) = 0.0; // ghost cell in upper absorbing boundary in z
        ndis(i, -1) = ndis(i, 0);     // ghost cell in lower reflecting boundary in z
        ndis(i, -2) = ndis(i, 1);     // ghost cell in lower reflecting boundary in z
    }
    for (int j = -Field2D::NG; j < ndis.nz() + Field2D::NG; j++) {
        ndis(ndis.nR(), j) = 0.0;     // ghost cell in outer absorbing boundary in R
        ndis(ndis.nR() + 1, j) = 0.0; // ghost cell in outer absorbing boundary in R
        ndis(-1, j) = ndis(0, j);     // ghost cell in inner reflecting boundary in R
        ndis(-2, j) = ndis(1, j);     // ghost cell in inner reflecting boundary in R
    }
}

void write_array_to_bin(const std::string& filename, Field2D& arr, const std::size_t size) {
    if (size > arr.data.size()) {
        throw std::invalid_argument("Requested output exceeds field storage");
    }
    // std::ios::binary 表示以二进制模式打开
    // std::ios::app    追加模式; std::ios::trunc    覆写模式
    std::ofstream file(filename, std::ios::binary | std::ios::app);
    
    if (!file) {
        // 使用 std::cerr 代替 perror
        std::cerr << "File opening failed: " << filename << std::endl;
        return;
    }
    
    // reinterpret_cast 将 double 指针转换为 char 指针，符合 write 接口
    file.write(reinterpret_cast<const char*>(arr.data.data()), size * sizeof(double));
    
    // file 在离开作用域时会自动关闭，不需要显式调用 file.close()
}

namespace {
// Reconstruct about the volume centroid, since n stores cell averages.
double reconstruction_center(const AxisGrid& axis, int k, bool radial)
{
    const double c = axis.center(k);
    return radial ? c + axis.width(k) * (axis.width(k) / c) / 12.0 : c;
}

double face_state(const Field2D& n, const AxisGrid& axis, int k, int transverse,
                  bool radial, double face)
{
    auto value = [&](int cell) { return radial ? n(cell, transverse) : n(transverse, cell); };
    // Vacuum outside the upper boundary; constant reconstruction at domain edges.
    if (k < 0 || k >= axis.size()) return 0.0;
    const double q = value(k);
    if (k == 0 || k == axis.size() - 1) return q;
    const double x = reconstruction_center(axis, k, radial);
    const double dl = q - value(k - 1), dr = value(k + 1) - q;
    const bool extremum = dl == 0.0 || dr == 0.0 || std::signbit(dl) != std::signbit(dr);
    if (extremum && (k < 2 || k + 2 >= axis.size())) return q;
    auto variance = [&](int cell) {
        const double w = axis.width(cell);
        const double offset = reconstruction_center(axis, cell, radial) - axis.center(cell);
        return w * w / 12.0 - offset * offset;
    };
    // Fit a quadratic to three volume-weighted cell averages. This gives
    // (-q[k-1] + 5*q[k] + 2*q[k+1])/6 at the right face on a uniform z grid.
    const double xm = reconstruction_center(axis, k - 1, radial) - x;
    const double xp = reconstruction_center(axis, k + 1, radial) - x;
    const double vm = xm * xm + variance(k - 1) - variance(k);
    const double vp = xp * xp + variance(k + 1) - variance(k);
    const double curvature = (dr / xp + dl / xm) / (vp / xp - vm / xm);
    const double slope = (dr - curvature * vp) / xp;
    const double dx = face - x;
    const double increment = slope * dx + curvature * (dx * dx - variance(k));
    if (k >= 2 && k + 2 < axis.size()) {
        // Preserve resolved smooth extrema only when the curvature agrees on
        // three neighboring stencils. Discontinuities and unresolved peaks
        // revert to constant reconstruction. No absolute smoothness threshold.
        auto second = [&](int cell) {
            const double xl = reconstruction_center(axis, cell, radial) - reconstruction_center(axis, cell - 1, radial);
            const double xr = reconstruction_center(axis, cell + 1, radial) - reconstruction_center(axis, cell, radial);
            return ((value(cell + 1) - value(cell)) / xr -
                    (value(cell) - value(cell - 1)) / xl) / (xl + xr);
        };
        const double cm = second(k - 1), cc = second(k), cp = second(k + 1);
        if (cm != 0.0 && cc != 0.0 && cp != 0.0 &&
            std::signbit(cm) == std::signbit(cc) && std::signbit(cp) == std::signbit(cc) &&
            std::max({std::abs(cm), std::abs(cc), std::abs(cp)}) <=
                2.0 * std::min({std::abs(cm), std::abs(cc), std::abs(cp)}))
            return std::clamp(q + increment, 0.0, 2.0 * q);
    }
    if (extremum) return q;
    // Koren-type face limiting: no new extrema and states in [0, 2*q].
    // Limiting uses differences, so even arbitrarily small densities work.
    const double sign = std::copysign(1.0, dl) * std::copysign(1.0, dx);
    const double limited = sign * std::max(0.0, std::min({sign * increment,
        std::abs(dl), std::abs(dr), q}));
    return std::clamp(q + limited, 0.0, 2.0 * q);
}

void advance(Field2D& n, Field2D& temp, const Field2D* vr, const Field2D* vz,
             double dt, const Grid2D& grid)
{
    validate_legacy_field(n, grid);
    validate_legacy_field(temp, grid);
    if (&n == &temp || !std::isfinite(dt) || dt < 0.0)
        throw std::invalid_argument("Advection requires distinct fields and finite nonnegative dt");
    if (vr) validate_legacy_field(*vr, grid, Field2D::Location::RadialFace);
    if (vz) validate_legacy_field(*vz, grid, Field2D::Location::VerticalFace);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j)
            if (!std::isfinite(n(i,j)) || n(i,j) < 0.0)
                throw std::invalid_argument("Density must be finite and nonnegative");
    if (vr)
        for (int i = 0; i <= grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                if (!std::isfinite((*vr)(i,j))) throw std::invalid_argument("Nonfinite radial velocity");
    if (vz)
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j <= grid.nz(); ++j)
                if (!std::isfinite((*vz)(i,j))) throw std::invalid_argument("Nonfinite vertical velocity");

    // Outgoing face area * speed / cell volume, summed over both directions.
    // Face states <= 2*n, so CFL <= 1/2 preserves positivity at each Euler stage.
    // Reflecting lower boundaries have zero flux regardless of supplied velocity.
    double max_rate = 0.0;
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) {
            double rate = 0.0;
            if (vr) rate += (grid.R().face(i + 1) * std::max(0.0, (*vr)(i + 1,j)) +
                (i ? grid.R().face(i) * std::max(0.0, -(*vr)(i,j)) : 0.0)) /
                (grid.R().center(i) * grid.R().width(i));
            if (vz) rate += (std::max(0.0, (*vz)(i,j + 1)) +
                (j ? std::max(0.0, -(*vz)(i,j)) : 0.0)) / grid.z().width(j);
            if (!std::isfinite(rate)) throw std::invalid_argument("Nonfinite cell crossing rate");
            max_rate = std::max(max_rate, rate);
        }
    const double count = std::max(1.0, std::ceil(dt * max_rate / 0.45));
    if (!std::isfinite(count) || count >= static_cast<double>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Requested interval requires too many CFL substeps");
    apply_boundary_conditions(n, grid);
    if (dt == 0.0 || max_rate == 0.0) return;
    const int steps = static_cast<int>(count);
    const double h = dt / steps;
    Field2D fr(grid, Field2D::Location::RadialFace);
    Field2D fz(grid, Field2D::Location::VerticalFace);
    auto fluxes = [&](const Field2D& state) {
        // Compute each shared face flux once, ensuring conservative cancellation.
        if (vr)
            for (int i = 1; i <= grid.nR(); ++i)
                for (int j = 0; j < grid.nz(); ++j) {
                    const double v = (*vr)(i,j);
                    fr(i,j) = v == 0.0 ? 0.0 : grid.R().face(i) * v *
                        face_state(state, grid.R(), v >= 0.0 ? i - 1 : i, j, true, grid.R().face(i));
                }
        if (vz)
            for (int i = 0; i < grid.nR(); ++i)
                for (int j = 1; j <= grid.nz(); ++j) {
                    const double v = (*vz)(i,j);
                    fz(i,j) = v == 0.0 ? 0.0 : v * face_state(state, grid.z(), v >= 0.0 ? j - 1 : j, i, false, grid.z().face(j));
                }
    };
    auto increment = [&](int i, int j) {
        return h * ((fr(i,j) - fr(i + 1,j)) / (grid.R().center(i) * grid.R().width(i)) +
                    (fz(i,j) - fz(i,j + 1)) / grid.z().width(j));
    };
    for (int step = 0; step < steps; ++step) {
        // SSPRK(2,2): u1 = u + h L(u); u_new = (u + u1 + h L(u1))/2.
        fluxes(n);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j) temp(i,j) = n(i,j) + increment(i,j);
        apply_boundary_conditions(temp, grid);
        fluxes(temp);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j) {
                n(i,j) = 0.5 * n(i,j) + 0.5 * (temp(i,j) + increment(i,j));
                if (!std::isfinite(n(i,j)) || n(i,j) < 0.0)
                    throw std::runtime_error("Invalid density after advection step");
            }
        apply_boundary_conditions(n, grid);
    }
}
} // namespace

void advection_TVD(Field2D& n, Field2D& temp, Field2D& vr, Field2D& vz, double dt, const Grid2D& grid)
{
    advance(n, temp, &vr, &vz, dt, grid);
}

void advection_TVD_R(Field2D& n, Field2D& temp, Field2D& vr, double dt, const Grid2D& grid)
{
    advance(n, temp, &vr, nullptr, dt, grid);
}

void advection_TVD_Z(Field2D& n, Field2D& temp, Field2D& vz, double dt, const Grid2D& grid)
{
    advance(n, temp, nullptr, &vz, dt, grid);
}


void solve_advection_equation(Field2D& ndis_C, Field2D& ndis_H, Field2D& vR, Field2D& vZ, const Grid2D& grid) {
    validate_legacy_field(ndis_C, grid);
    validate_legacy_field(ndis_H, grid);
    validate_legacy_field(vR, grid, Field2D::Location::RadialFace);
    validate_legacy_field(vZ, grid, Field2D::Location::VerticalFace);

    Field2D temp(grid);
    // Unsplit second-order update includes both flux divergences at each stage.
    for (int step = 0; step < NT; ++step)
        advection_TVD(ndis_C, temp, vR, vZ, DT, grid);

    write_array_to_bin("ndis_C_t_ssad.bin", ndis_C, ndis_C.data.size());

}

