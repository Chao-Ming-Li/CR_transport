#include "CR_advection.hpp"
#include "CR_transport_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cr_advection {
namespace {

double limited_slope(double left, double center, double right, Limiter limiter)
{
    if (!((left > 0 && right > 0) || (left < 0 && right < 0))) return 0.0;
    const double small = std::min(std::abs(left), std::abs(right));
    const double large = std::max(std::abs(left), std::abs(right));
    switch (limiter) {
        case Limiter::Minmod: return std::copysign(small, left);
        case Limiter::VanLeer:
            return std::copysign(small * (2.0 / (1.0 + small / large)), left);
        case Limiter::MC:
            return std::copysign(std::min(2.0 * small, std::abs(center)), left);
    }
    throw std::invalid_argument("Unknown slope limiter");
}

// Reconstruct from the upwind cell's volume centroid to the requested face.
// The same implementation handles R and z, including two ghost-cell layers.
double face_value(const Field2D& n, const Grid2D& g, const Options& opt,
                  bool radial, int face, int transverse, double velocity)
{
    const int cell = velocity >= 0 ? face - 1 : face;
    const auto value = [&](int k) {
        return radial ? n(k, transverse) : n(transverse, k);
    };
    if (opt.reconstruction == Reconstruction::Constant) return value(cell);
    const auto coordinate = [&](int k) {
        return radial ? g.radial_centroid(k) : g.z().center(k);
    };
    const double dl = coordinate(cell) - coordinate(cell - 1);
    const double dr = coordinate(cell + 1) - coordinate(cell);
    const double sl = (value(cell) - value(cell - 1)) / dl;
    const double sr = (value(cell + 1) - value(cell)) / dr;
    const double sc = (value(cell + 1) - value(cell - 1)) / (dl + dr);
    const double xf = radial ? g.R().face(face) : g.z().face(face);
    return value(cell) + limited_slope(sl, sc, sr, opt.limiter)
                         * (xf - coordinate(cell));
}

// Return L(n), not a complete Euler update. Compute each shared face once.
void spatial_operator(Field2D& n, Field2D& rhs, Field2D& fr, Field2D& fz,
                      const Field2D& vr, const Field2D& vz, const Grid2D& g,
                      const Options& opt, bool radial, bool vertical)
{
    apply_boundary_conditions(n, g);
    if (radial) {
        for (int i = 1; i <= g.nR(); ++i)
            for (int j = 0; j < g.nz(); ++j)
                fr(i, j) = g.R().face(i) * vr(i, j)
                    * face_value(n, g, opt, true, i, j, vr(i, j));
    }
    if (vertical) {
        for (int i = 0; i < g.nR(); ++i)
            for (int j = 1; j <= g.nz(); ++j)
                fz(i, j) = vz(i, j)
                    * face_value(n, g, opt, false, j, i, vz(i, j));
    }
    for (int i = 0; i < g.nR(); ++i)
        for (int j = 0; j < g.nz(); ++j) {
            rhs(i, j) = 0.0;
            if (radial) rhs(i, j) -= (fr(i + 1, j) - fr(i, j))
                / (g.R().center(i) * g.R().width(i));
            if (vertical) rhs(i, j) -= (fz(i, j + 1) - fz(i, j))
                / g.z().width(j);
        }
}

} // namespace

void advance(Field2D& density, const Field2D& vR, const Field2D& vZ,
             double dt, const Grid2D& grid, const Options& opt)
{

    Field2D stage = density;
    Field2D rhs(grid), fr(grid, Field2D::Location::RadialFace),
        fz(grid, Field2D::Location::VerticalFace);
    const auto euler = [&](bool radial, bool vertical) {
        spatial_operator(stage, rhs, fr, fz, vR, vZ, grid, opt, radial, vertical);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                stage(i, j) += dt * rhs(i, j);
    };
    if (opt.splitting == Splitting::RadialThenVertical) {
        euler(true, false);
        euler(false, true);
    } else {
        euler(true, true);
        if (opt.integrator == TimeIntegrator::SSPRK2) {
            euler(true, true);
            for (int i = 0; i < grid.nR(); ++i)
                for (int j = 0; j < grid.nz(); ++j)
                    stage(i, j) = 0.5 * density(i, j) + 0.5 * stage(i, j);
        }
    }
    apply_boundary_conditions(stage, grid);
    density = std::move(stage);
}

} // namespace cr_advection
