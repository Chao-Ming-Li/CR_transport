#pragma once

#include "field.hpp"

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
    // Velocities are face fields held fixed during the step. Lower faces are
    // reflecting (zero flux); outer ghost cells are zero, as in the legacy solver.
    // Caller supplies a stable dt. No density clipping or automatic timestep
    // selection is performed. Input dimensions and timestep are validated.
    // Constant reconstruction ignores limiter. RadialThenVertical supports Euler
    // only and reproduces the legacy sequential directional update.
    // VanLeer uses sign-aware slopes without the legacy absolute 1e-12 offset.
    // Owns a grid snapshot and reusable stage/RHS/flux storage. Each instance is
    // intended for sequential use; construct separate instances for concurrent calls.
    class Solver {
        public:
            explicit Solver(const Grid2D &grid, const Options &options = {});
            void advance(Field2D &density, const Field2D &vR, const Field2D &vZ, double dt);

        private:
            Grid2D grid_;
            Options options_;
            Field2D stage_, rhs_, fr_, fz_;
    };

} // namespace cr_advection

namespace cr_diffusion {

    enum class Boundary { Reflecting, Absorbing };

    struct Options {
        Boundary outer_R = Boundary::Absorbing;
        Boundary upper_z = Boundary::Absorbing;
    };

    // Constant isotropic diffusion: dn/dt = (1/R)d_R(R D d_R n) + d_z(D d_z n).
    // Peaceman–Rachford Crank–Nicolson ADI: implicit R/explicit z, then
    // implicit z/explicit R. Lower faces have zero flux. Absorbing outer
    // boundaries use zero ghost-cell density, matching the advection convention.
    // D has units length^2/time; dt has units time. Both must be finite and >= 0.
    // Second-order in time for fixed D. CN is not positivity preserving at large
    // timesteps; signed results are retained without clipping. No source/loss terms.
    // Fixed grid, D, dt and boundaries. Construct a new solver when these change.
    // Owns reusable fields and cached Thomas factors; no timestep allocations.
    class Solver {
    public:
        Solver(const Grid2D &grid, double D, double dt, const Options &options = {});
        void advance(Field2D &density);

    private:
        struct Axis {
            std::vector<double> lower, upper, a, cp, inverse_pivot;
            explicit Axis(int n);
            void factor(double half_dt);
            void solve(std::vector<double> &line) const;
        };
        Grid2D grid_;
        Options options_;
        double D_, half_dt_;
        Axis radial_, vertical_;
        Field2D initial_, intermediate_, result_;
        std::vector<double> line_;
    };

} // namespace cr_diffusion
