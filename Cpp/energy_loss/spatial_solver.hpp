#pragma once

#include "field.hpp"
#include "initialization.hpp"

// Ported from advection_grid/CR_advection.hpp. Spatial units: kpc, yr.
// Field3D overloads evolve independent rigidity planes in parallel with OpenMP.
// Caller supplies matching grids; instances are sequential.

namespace cr_advection {

    enum class Reconstruction { Constant, PLM };
    enum class Limiter { Minmod, VanLeer, MC };
    enum class TimeIntegrator { Euler, SSPRK2 };
    enum class Splitting { Unsplit, RadialThenVertical };

    struct Options {
        Reconstruction reconstruction = Reconstruction::PLM;
        Limiter limiter = Limiter::VanLeer;
        TimeIntegrator integrator = TimeIntegrator::SSPRK2;
        Splitting splitting = Splitting::Unsplit;
    };

    // Conservative cylindrical transport of a nonnegative cell-average density.
    // Velocities are face fields copied at construction, together with a fixed dt.
    // Lower faces are reflecting (zero flux); outer ghost cells are zero,
    // as in the legacy solver.
    // Caller supplies a stable dt. No density clipping or automatic timestep
    // selection is performed. Input dimensions and timestep are validated.
    // Constant reconstruction ignores limiter. RadialThenVertical supports Euler
    // only and reproduces the legacy sequential directional update.
    // VanLeer uses sign-aware slopes without the legacy absolute 1e-12 offset.
    // Owns grid/velocity snapshots and reusable stage/RHS/flux storage.
    // Construct a new solver when the grid, velocities, timestep or options change.
    // Each instance is intended for sequential use; construct separate instances
    // for concurrent calls.
    class Solver {
    public:
        Solver(const SpatialGrid &grid, const Field2D &vR, const Field2D &vZ, double dt, const Options &options = {});
        void advance(Field2D &density);
        void advance(Field3D &density);

    private:
        struct Workspace {
            explicit Workspace(const SpatialGrid &grid);
            Field2D stage, rhs, fr, fz, plane;
        };
        void advance_plane(Field2D &density, Workspace &work);
        void euler_step(Workspace &work, bool radial, bool vertical);
        SpatialGrid grid_;
        Options options_;
        Field2D vR_, vZ_;
        double dt_;
        std::vector<Workspace> work_;
    };

} // namespace cr_advection

namespace cr_diffusion {

    // Spatially uniform isotropic diffusion: dn/dt = (1/R)d_R(R D d_R n) + d_z(D d_z n).
    // Peaceman–Rachford Crank–Nicolson ADI: implicit R/explicit z, then
    // implicit z/explicit R. Lower faces have zero flux. Absorbing outer
    // boundaries use zero ghost-cell density, matching the advection convention.
    // D is in kpc^2/yr; dt is in yr. Caller supplies finite, nonnegative values.
    // The model overload uses D(Rg) at each rigidity center.
    // Second-order in time for fixed D. CN is not positivity preserving at large
    // timesteps; signed results are retained without clipping. No source/loss terms.
    // Fixed reflecting lower and absorbing outer boundaries.
    // Fixed grids, coefficients and dt. Construct a new solver when these change.
    // Shares cached Thomas factors across threads and reuses separate workspaces.
    class Solver {
    public:
        Solver(const SpatialGrid &grid, double D, double dt);
        // Cache separate ADI factors at each rigidity center. Field3D only.
        Solver(const SpatialGrid &grid, const RgGrid &Rg, const DiffusionModel &model, double dt);
        void advance(Field2D &density);
        void advance(Field3D &density);

    private:
        struct Axis {
            std::vector<double> lower, upper, a, cp, inverse_pivot;
            explicit Axis(int n);
            void factor(double half_dt);
            void solve(std::vector<double> &line) const;
        };
        SpatialGrid grid_;
        double D_, half_dt_;
        Axis radial_, vertical_;
        std::vector<Axis> rigidity_radial_, rigidity_vertical_;
        void initialize_axis(Axis &axis, double D, bool radial);
        struct Workspace {
            explicit Workspace(const SpatialGrid &grid);
            Field2D initial, intermediate, result, plane;
            std::vector<double> line;
        };
        void advance_plane(Field2D &density, const Axis &radial, const Axis &vertical, Workspace &work);
        std::vector<Workspace> work_;
    };

} // namespace cr_diffusion
