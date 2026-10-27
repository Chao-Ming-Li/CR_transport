#pragma once

#include "field.hpp"
#include <span>
#include <vector>

namespace cr_momentum {

    enum class Boundary {
        ZeroFlux,
        Outflow // Adiabatic outflow only; no inflow, zero diffusive flux.
    };

    // Evolves the cell-average q = dn/dRg at fixed spatial positions:
    // dq/dt = -dJ/dRg,
    // J = -(div(v)/3)*Rg*q - Rg^2*Dpp*d(q/Rg^2)/dRg.
    // Units: kpc, yr, GV; Dpp is in GV^2/yr, NOT D_RgRg.
    // Upwind adiabatic flux, two-point diffusion, Crank-Nicolson in time.
    // CN does not guarantee nonnegative densities for arbitrary timesteps.
    // Face coefficients are spatially uniform and copied at construction.
    // Zero diffusion is supplied as an array of zeros of length Rg.size()+1.
    // Grid, wind, dt and coefficients are fixed snapshots. Reconstruct when
    // they change. Reuses line storage; an instance is for sequential calls.
    class Solver {
    public:
        Solver(const SpatialGrid &grid, const RgGrid &Rg,
               const Field2D &vR, const Field2D &vZ,
               std::span<const double> Dpp, double dt,
               Boundary lower = Boundary::ZeroFlux,
               Boundary upper = Boundary::ZeroFlux);
        // Updates physical spatial cells only; spatial ghosts are untouched.
        void advance(Field3D &q);

    private:
        int nR_, nz_, n_;
        double dt_;
        Boundary lower_boundary_, upper_boundary_;
        std::vector<double> divergence_, faces_, widths_, diff_left_, diff_right_;
        std::vector<double> a_, b_, c_, rhs_;
    };

} // namespace cr_momentum
