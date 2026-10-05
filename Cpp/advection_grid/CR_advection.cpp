#include "CR_advection.hpp"
#include "CR_transport_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>


namespace cr_advection {

    Solver::Solver(const Grid2D &grid, const Options &options)
        : grid_(grid), options_(options), stage_(grid), rhs_(grid),
          fr_(grid, Field2D::Location::RadialFace), fz_(grid, Field2D::Location::VerticalFace) {}

    namespace {

        double limited_slope(double left, double center, double right, Limiter limiter) {
            if (!((left > 0 && right > 0) || (left < 0 && right < 0)))
                return 0.0;
            const double small = std::min(std::abs(left), std::abs(right));
            const double large = std::max(std::abs(left), std::abs(right));
            switch (limiter) {
                case Limiter::Minmod:
                    return std::copysign(small, left);
                case Limiter::VanLeer:
                    return std::copysign(small * (2.0 / (1.0 + small / large)), left);
                case Limiter::MC:
                    return std::copysign(std::min(2.0 * small, std::abs(center)), left);
            }
            // Solver construction has already validated the limiter.
            throw std::logic_error("Unexpected slope limiter");
        }

        // Reconstruct from the upwind cell's volume centroid to the requested face.
        // The same implementation handles R and z, including two ghost-cell layers.
        double face_value(const Field2D &n, const Grid2D &g, const Options &opt, bool radial,
                          int face, int transverse, double velocity) {
            const int cell = velocity >= 0 ? face - 1 : face;
            const auto value = [&](int k) { return radial ? n(k, transverse) : n(transverse, k); };
            if (opt.reconstruction == Reconstruction::Constant)
                return value(
                    cell); // upwind method, the face value is the same as the upwind cell value
            const auto coordinate = [&](int k) {
                return radial ? g.radial_centroid(k) : g.z().center(k);
            };
            const double dl = coordinate(cell) - coordinate(cell - 1);
            const double dr = coordinate(cell + 1) - coordinate(cell);
            const double sl = (value(cell) - value(cell - 1)) / dl;
            const double sr = (value(cell + 1) - value(cell)) / dr;
            const double sc = (value(cell + 1) - value(cell - 1)) / (dl + dr);
            const double xf = radial ? g.R().face(face) : g.z().face(face); // face coordinate
            const double reconstructed = value(cell) + limited_slope(sl, sc, sr, opt.limiter) *
                    (xf - coordinate(cell)); // linear reconstruction from the cell centroid to the face

            // Keep the face state between its adjacent cell averages. Near a steep
            // drop, roundoff in the coordinates can otherwise make it slightly negative.
            const double neighbor = value(cell + (velocity >= 0 ? 1 : -1));
            return std::clamp(reconstructed, std::min(value(cell), neighbor), std::max(value(cell),
                         neighbor)); // std::clamp(x, lower, upper) keeps x within [lower, upper]
        }

        // Return L(n), not a complete Euler update. Compute each shared face once.
        void spatial_operator(Field2D &n, Field2D &rhs, Field2D &fr, Field2D &fz, const Field2D &vr,
                              const Field2D &vz, const Grid2D &g, const Options &opt, bool radial,
                              bool vertical) {
            apply_boundary_conditions(n, g);
            if (radial) {
                for (int i = 1; i <= g.nR(); ++i) // radial faces are indexed 1..nR, excluding Rf[0] because flux is zero at the inner boundary, including the Rf[nR] face.
                    for (int j = 0; j < g.nz(); ++j) // vertical indexes of radial faces are 0..nz-1
                        fr(i, j) = g.R().face(i) * vr(i, j) * face_value(n, g, opt, true, i, j, vr(i, j));
            }
            if (vertical) {
                for (int i = 0; i < g.nR(); ++i)
                    for (int j = 1; j <= g.nz(); ++j)
                        fz(i, j) = vz(i, j) * face_value(n, g, opt, false, j, i, vz(i, j));
            }
            for (int i = 0; i < g.nR(); ++i)
                for (int j = 0; j < g.nz(); ++j) {
                    rhs(i, j) = 0.0;
                    if (radial)
                        rhs(i, j) -= (fr(i + 1, j) - fr(i, j)) / (g.R().center(i) * g.R().width(i));
                    if (vertical)
                        rhs(i, j) -= (fz(i, j + 1) - fz(i, j)) / g.z().width(j);
                }
        }

    } // namespace

    void Solver::advance(Field2D &density, const Field2D &vR, const Field2D &vZ, double dt) {

        const auto &grid = grid_;
        const auto &opt = options_;
        auto &stage = stage_;
        auto &rhs = rhs_;
        auto &fr = fr_;
        auto &fz = fz_;
        stage.data = density.data;
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
        density.data.swap(stage.data);
    }

} // namespace cr_advection

namespace cr_diffusion {
    namespace {

        // L(n)_k = lower_k (n_{k-1}-n_k) + upper_k (n_{k+1}-n_k).
        // Shared face conductances divided by cell volume ensure conservation.
        struct AxisCoefficients {
            std::vector<double> lower, upper;
            explicit AxisCoefficients(int n) : lower(n, 0.0), upper(n, 0.0) {}
        };

        AxisCoefficients radial_coefficients(const Grid2D &grid, double D, Boundary outer) {
            AxisCoefficients coefficients(grid.nR());
            for (int i = 0; i < grid.nR(); ++i) {
                const double volume = grid.R().center(i) * grid.R().width(i);
                if (i > 0)
                    coefficients.lower[i] =
                        D * grid.R().face(i) / (volume * grid.radial_centroid_distance(i - 1));
                if (i + 1 < grid.nR() || outer == Boundary::Absorbing)
                    coefficients.upper[i] =
                        D * grid.R().face(i + 1) / (volume * grid.radial_centroid_distance(i));
            }
            return coefficients;
        }

        AxisCoefficients vertical_coefficients(const Grid2D &grid, double D, Boundary outer) {
            AxisCoefficients coefficients(grid.nz());
            for (int j = 0; j < grid.nz(); ++j) {
                if (j > 0)
                    coefficients.lower[j] =
                        D / (grid.z().width(j) * grid.z().center_distance(j - 1));
                if (j + 1 < grid.nz() || outer == Boundary::Absorbing)
                    coefficients.upper[j] = D / (grid.z().width(j) * grid.z().center_distance(j));
            }
            return coefficients;
        }

    } // namespace

    Solver::Axis::Axis(int n) : lower(n), upper(n), a(n), cp(n), inverse_pivot(n) {}

    // Thomas algorithm factors the tridiagonal matrix for each axis. The lower and upper
    // coefficients are the same for every timestep, so factor once and reuse. The diagonal
    // is 1 + dt/2 * (lower + upper). The half_dt argument is dt/2, not dt, because the ADI method uses two half-steps per timestep. The factorization is in-place, overwriting the lower and upper coefficients with the Thomas factors.
    void Solver::Axis::factor(double half_dt) {
        for (std::size_t k = 0; k < lower.size(); ++k) {
            a[k] = -half_dt * lower[k];
            const double b = 1.0 + half_dt * (lower[k] + upper[k]);
            const double pivot = b - (k > 0 ? a[k] * cp[k - 1] : 0.0);
            inverse_pivot[k] = 1.0 / pivot;
            cp[k] = k + 1 < lower.size() ? -half_dt * upper[k] * inverse_pivot[k] : 0.0;
        }
    }

    void Solver::Axis::solve(std::vector<double> &line) const {
        // Forward substitution followed by back substitution, entirely in-place.
        line[0] *= inverse_pivot[0];
        for (std::size_t k = 1; k < lower.size(); ++k)
            line[k] = (line[k] - a[k] * line[k - 1]) * inverse_pivot[k];
        for (int k = static_cast<int>(lower.size()) - 2; k >= 0; --k)
            line[k] -= cp[k] * line[k + 1];
    }

    Solver::Solver(const Grid2D &grid, double D, double dt, const Options &options)
        : grid_(grid), options_(options), D_(D), half_dt_(0.5 * dt), radial_(grid.nR()),
          vertical_(grid.nz()), initial_(grid), intermediate_(grid), result_(grid),
          line_(std::max(grid.nR(), grid.nz())) {

        auto radial = radial_coefficients(grid_, D, options.outer_R);
        auto vertical = vertical_coefficients(grid_, D, options.upper_z);
        radial_.lower = std::move(radial.lower);
        radial_.upper = std::move(radial.upper);
        vertical_.lower = std::move(vertical.lower);
        vertical_.upper = std::move(vertical.upper);
        radial_.factor(half_dt_);
        vertical_.factor(half_dt_);
    }

    void Solver::advance(Field2D &density) {
        if (D_ == 0.0)   return;

        const auto &grid = grid_;
        const auto &options = options_;
        const auto &radial = radial_;
        const auto &vertical = vertical_;
        const double half_dt = half_dt_;
        auto &initial = initial_;
        auto &intermediate = intermediate_;
        auto &result = result_;
        initial.data = density.data;
        apply_boundary_conditions(initial, grid);
        // Step 1: (I - dt/2 L_R) intermediate = (I + dt/2 L_z) initial.
        for (int j = 0; j < grid.nz(); ++j) {
            for (int i = 0; i < grid.nR(); ++i)
                line_[i] = initial(i, j) + half_dt * (vertical.lower[j] * (initial(i, j - 1) - initial(i, j)) + vertical.upper[j] * (initial(i, j + 1) - initial(i, j)));
            radial.solve(line_);
            for (int i = 0; i < grid.nR(); ++i)
                intermediate(i, j) = line_[i];
        }
        apply_boundary_conditions(intermediate, grid);

        // Step 2: (I - dt/2 L_z) result = (I + dt/2 L_R) intermediate.
        for (int i = 0; i < grid.nR(); ++i) {
            for (int j = 0; j < grid.nz(); ++j)
                line_[j] = intermediate(i, j) + half_dt * (radial.lower[i] * (intermediate(i - 1, j) - intermediate(i, j)) + radial.upper[i] * (intermediate(i + 1, j) - intermediate(i, j)));
            vertical.solve(line_);
            for (int j = 0; j < grid.nz(); ++j) {
                result(i, j) = line_[j];
            }
        }
        apply_boundary_conditions(result, grid);

        density.data.swap(result.data);
    }

} // namespace cr_diffusion
