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
// Caller supplies a stable dt; negative/nonfinite stage densities throw without
// modifying density. No clipping or automatic timestep selection is performed.
// Constant reconstruction ignores limiter. RadialThenVertical supports Euler
// only and reproduces the legacy sequential directional update.
// VanLeer uses sign-aware slopes without the legacy absolute 1e-12 offset.
void advance(Field2D& density, const Field2D& vR, const Field2D& vZ,
             double dt, const Grid2D& grid, const Options& options = {});

} // namespace cr_advection
