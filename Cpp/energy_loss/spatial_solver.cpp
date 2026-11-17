#include "spatial_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

    void apply_boundary_conditions(Field2D &ndis) {
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

} // namespace

namespace cr_advection {

    Solver::Workspace::Workspace(const SpatialGrid &grid)
        : stage(grid), rhs(grid), fr(grid, Field2D::Location::RadialFace),
          fz(grid, Field2D::Location::VerticalFace), plane(grid) {}


    Solver::Solver(const SpatialGrid &grid, const Field2D &vR, const Field2D &vZ, double dt, const Options &options)
        : grid_(grid), options_(options), vR_(vR), vZ_(vZ), dt_(dt) {
        work_.emplace_back(grid);
        if (!vR.matches(grid, Field2D::Location::RadialFace) || !vZ.matches(grid, Field2D::Location::VerticalFace))
            throw std::invalid_argument("Advection velocities must match the grid");
        if (!std::isfinite(dt) || dt < 0.0)
            throw std::invalid_argument("Advection timestep must be finite and nonnegative");
        for (double v : vR_.data) {
            if (!std::isfinite(v))
                throw std::invalid_argument("Radial velocity must be finite");
        }
        for (double v : vZ_.data) {
            if (!std::isfinite(v))
                throw std::invalid_argument("Vertical velocity must be finite");
        }
    }

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
            throw std::logic_error("Unexpected slope limiter");
        }

        // Reconstruct from the upwind cell's volume centroid to the requested face.
        // The same implementation handles R and z, including two ghost-cell layers.
        double face_value(const Field2D &n, const SpatialGrid &g, const Options &opt, bool radial, int face, int transverse, double velocity) {
            const int cell = velocity >= 0 ? face - 1 : face;
            const double center_value = radial ? n(cell, transverse) : n(transverse, cell);
            if (opt.reconstruction == Reconstruction::Constant)
                return center_value;

            const double left_value = radial ? n(cell - 1, transverse) : n(transverse, cell - 1);
            const double right_value = radial ? n(cell + 1, transverse) : n(transverse, cell + 1);
            const double xc = radial ? g.radial_centroid(cell) : g.z().center(cell);
            const double xl = radial ? g.radial_centroid(cell - 1) : g.z().center(cell - 1);
            const double xr = radial ? g.radial_centroid(cell + 1) : g.z().center(cell + 1);
            const double xf = radial ? g.R().face(face) : g.z().face(face);
            const double dl = xc - xl;
            const double dr = xr - xc;
            const double sl = (center_value - left_value) / dl;
            const double sr = (right_value - center_value) / dr;
            const double sc = (right_value - left_value) / (dl + dr);
            const double reconstructed = center_value + limited_slope(sl, sc, sr, opt.limiter) * (xf - xc);

            // Limit the face state to adjacent cell averages to avoid overshoots.
            const double neighbor = velocity >= 0 ? right_value : left_value;
            return std::clamp(reconstructed, std::min(center_value, neighbor), std::max(center_value, neighbor));
        }

        // Return L(n), not a complete Euler update. Compute each shared face once.
        void spatial_operator(Field2D &n, Field2D &rhs, Field2D &fr, Field2D &fz, const Field2D &vr, const Field2D &vz, const SpatialGrid &g, const Options &opt, bool radial, bool vertical) {
            apply_boundary_conditions(n);
            if (radial) {
                // Radial faces: 1..nR, including the outer face Rf[nR].
                // Exclude Rf[0], where the inner boundary flux is zero.
                // Vertical indexes of radial faces: 0..nz-1.
                for (int i = 1; i <= g.nR(); ++i) {
                    for (int j = 0; j < g.nz(); ++j) {
                        fr(i, j) = g.R().face(i) * vr(i, j) * face_value(n, g, opt, true, i, j, vr(i, j));
                    }
                }
            }
            if (vertical) {
                for (int i = 0; i < g.nR(); ++i) {
                    for (int j = 1; j <= g.nz(); ++j) {
                        fz(i, j) = vz(i, j) * face_value(n, g, opt, false, j, i, vz(i, j));
                    }
                }
            }
            for (int i = 0; i < g.nR(); ++i) {
                for (int j = 0; j < g.nz(); ++j) {
                    rhs(i, j) = 0.0;
                    if (radial)
                        rhs(i, j) -= (fr(i + 1, j) - fr(i, j)) / (g.R().center(i) * g.R().width(i));
                    if (vertical)
                        rhs(i, j) -= (fz(i, j + 1) - fz(i, j)) / g.z().width(j);
                }
            }
        }

    } // namespace

    // Euler stage: n <- n + dt * L(n), for the selected directions.
    void Solver::euler_step(Workspace &work, bool radial, bool vertical) {
        spatial_operator(work.stage, work.rhs, work.fr, work.fz, vR_, vZ_, grid_, options_, radial, vertical);
        for (int i = 0; i < grid_.nR(); ++i) {
            for (int j = 0; j < grid_.nz(); ++j) {
                work.stage(i, j) += dt_ * work.rhs(i, j);
            }
        }
    }

    void Solver::advance(Field2D &density) {
        if (!density.matches(grid_, Field2D::Location::Centroid))
            throw std::invalid_argument("Advection density must match the grid");
        if (dt_ == 0.0)
            return;
        advance_plane(density, work_.front());
    }

    void Solver::advance_plane(Field2D &density, Workspace &work) {
        const SpatialGrid &grid = grid_;
        const Options &opt = options_;
        Field2D &stage = work.stage;
        stage.data = density.data;
        if (opt.splitting == Splitting::RadialThenVertical) {
            euler_step(work, true, false);
            euler_step(work, false, true);
        } else {
            euler_step(work, true, true);
            if (opt.integrator == TimeIntegrator::SSPRK2) {
                euler_step(work, true, true);
                for (int i = 0; i < grid.nR(); ++i) {
                    for (int j = 0; j < grid.nz(); ++j) {
                        stage(i, j) = 0.5 * density(i, j) + 0.5 * stage(i, j);
                    }
                }
            }
        }
        apply_boundary_conditions(stage);
        density.data.swap(stage.data);
    }

    void Solver::advance(Field3D &density) {
        if (density.nR() != grid_.nR() || density.nz() != grid_.nz())
            throw std::invalid_argument("Advection density must match the grid");
        if (dt_ == 0.0)
            return;
        // Include ghosts in each copy; the 2D solver refreshes them in place.
#ifdef _OPENMP
        const int threads = std::min(omp_get_max_threads(), density.nrigidity());
        while (work_.size() < static_cast<std::size_t>(threads)) {
            work_.emplace_back(grid_);
        }
        #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
        for (int k = 0; k < density.nrigidity(); ++k) {
#ifdef _OPENMP
            Workspace &work = work_[omp_get_thread_num()];
#else
            Workspace &work = work_.front();
#endif
            auto values = density.spatial_plane(k);
            std::copy(values.begin(), values.end(), work.plane.data.begin());
            advance_plane(work.plane, work);
            std::copy(work.plane.data.begin(), work.plane.data.end(), values.begin());
        }
    }

} // namespace cr_advection

namespace cr_diffusion {

    Solver::Axis::Axis(int n) : lower(n), upper(n), a(n), cp(n), inverse_pivot(n) {}

    // L(n)_k = lower_k*(n_(k-1)-n_k) + upper_k*(n_(k+1)-n_k).
    // Shared face conductances divided by cell volume ensure conservation.
    void Solver::initialize_axis(Axis &axis, double D, bool radial) {
        if (radial) {
            for (int i = 0; i < grid_.nR(); ++i) {
                // Annular volume divided by 2*pi*dz: (R_right^2-R_left^2)/2.
                const double radial_volume = grid_.R().center(i) * grid_.R().width(i);
                if (i > 0)
                    axis.lower[i] = D * grid_.R().face(i) / (radial_volume * grid_.radial_centroid_distance(i - 1));
                axis.upper[i] = D * grid_.R().face(i + 1) / (radial_volume * grid_.radial_centroid_distance(i));
            }
        } else {
            for (int j = 0; j < grid_.nz(); ++j) {
                if (j > 0)
                    axis.lower[j] = D / (grid_.z().width(j) * grid_.z().center_distance(j - 1));
                axis.upper[j] = D / (grid_.z().width(j) * grid_.z().center_distance(j));
            }
        }
        axis.factor(half_dt_);
    }

    // Thomas algorithm factors the tridiagonal matrix for each axis. The lower and upper
    // coefficients are the same for every timestep, so factor once and reuse. The diagonal
    // is 1 + dt/2 * (lower + upper). Each ADI stage uses dt/2 in its
    // implicit and explicit directional operators. Retain lower and upper
    // separately because they are also needed for the explicit terms.
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
        for (std::size_t k = 1; k < lower.size(); ++k) {
            line[k] = (line[k] - a[k] * line[k - 1]) * inverse_pivot[k];
        }
        for (int k = static_cast<int>(lower.size()) - 2; k >= 0; --k) {
            line[k] -= cp[k] * line[k + 1];
        }
    }

    Solver::Workspace::Workspace(const SpatialGrid &grid)
        : initial(grid), intermediate(grid), result(grid), plane(grid),
          line(std::max(grid.nR(), grid.nz())) {}

    Solver::Solver(const SpatialGrid &grid, double D, double dt)
        : grid_(grid), D_(D), half_dt_(0.5 * dt), radial_(grid.nR()), vertical_(grid.nz()) {
        work_.emplace_back(grid);
        initialize_axis(radial_, D, true);
        initialize_axis(vertical_, D, false);
    }

    Solver::Solver(const SpatialGrid &grid, const RgGrid &Rg, const DiffusionModel &model, double dt)
        : Solver(grid, 0.0, dt) {
        rigidity_radial_.reserve(Rg.size());
        rigidity_vertical_.reserve(Rg.size());
        for (int k = 0; k < Rg.size(); ++k) {
            const double D = model.spatial_coefficient(Rg.center(k));
            rigidity_radial_.push_back(Axis(grid.nR()));
            rigidity_vertical_.push_back(Axis(grid.nz()));
            initialize_axis(rigidity_radial_[k], D, true);
            initialize_axis(rigidity_vertical_[k], D, false);
        }
    }

    void Solver::advance(Field2D &density) {
        if (!rigidity_radial_.empty())
            throw std::invalid_argument("Rigidity-dependent diffusion requires Field3D");
        if (D_ == 0.0 || half_dt_ == 0.0)
            return;
        advance_plane(density, radial_, vertical_, work_.front());
    }

    void Solver::advance_plane(Field2D &density, const Axis &radial, const Axis &vertical, Workspace &work) {
        const SpatialGrid &grid = grid_;
        const double half_dt = half_dt_;
        Field2D &initial = work.initial;
        Field2D &intermediate = work.intermediate;
        Field2D &result = work.result;
        initial.data = density.data;
        apply_boundary_conditions(initial);
        // Step 1: (I - dt/2 L_R) intermediate = (I + dt/2 L_z) initial.
        for (int j = 0; j < grid.nz(); ++j) {
            for (int i = 0; i < grid.nR(); ++i) {
                const double vertical_diffusion =
                    vertical.lower[j] * (initial(i, j - 1) - initial(i, j))
                    + vertical.upper[j] * (initial(i, j + 1) - initial(i, j));
                work.line[i] = initial(i, j) + half_dt * vertical_diffusion;
            }
            radial.solve(work.line);
            for (int i = 0; i < grid.nR(); ++i) {
                intermediate(i, j) = work.line[i];
            }
        }
        apply_boundary_conditions(intermediate);

        // Step 2: (I - dt/2 L_z) result = (I + dt/2 L_R) intermediate.
        for (int i = 0; i < grid.nR(); ++i) {
            for (int j = 0; j < grid.nz(); ++j) {
                const double radial_diffusion =
                    radial.lower[i] * (intermediate(i - 1, j) - intermediate(i, j))
                    + radial.upper[i] * (intermediate(i + 1, j) - intermediate(i, j));
                work.line[j] = intermediate(i, j) + half_dt * radial_diffusion;
            }
            vertical.solve(work.line);
            for (int j = 0; j < grid.nz(); ++j) {
                result(i, j) = work.line[j];
            }
        }
        apply_boundary_conditions(result);

        density.data.swap(result.data);
    }

    void Solver::advance(Field3D &density) {
        if (density.nR() != grid_.nR() || density.nz() != grid_.nz() ||
            (!rigidity_radial_.empty() && density.nrigidity() != static_cast<int>(rigidity_radial_.size())))
            throw std::invalid_argument("Diffusion density must match construction grids");
        if ((D_ == 0.0 && rigidity_radial_.empty()) || half_dt_ == 0.0)
            return;
#ifdef _OPENMP
        const int threads = std::min(omp_get_max_threads(), density.nrigidity());
        while (work_.size() < static_cast<std::size_t>(threads)) {
            work_.emplace_back(grid_);
        }
        #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
        for (int k = 0; k < density.nrigidity(); ++k) {
#ifdef _OPENMP
            Workspace &work = work_[omp_get_thread_num()];
#else
            Workspace &work = work_.front();
#endif
            auto values = density.spatial_plane(k);
            std::copy(values.begin(), values.end(), work.plane.data.begin());
            if (rigidity_radial_.empty())
                advance_plane(work.plane, radial_, vertical_, work);
            else
                advance_plane(work.plane, rigidity_radial_[k], rigidity_vertical_[k], work);
            std::copy(work.plane.data.begin(), work.plane.data.end(), values.begin());
        }
    }

} // namespace cr_diffusion
