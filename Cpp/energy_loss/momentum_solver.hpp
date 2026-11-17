#pragma once

#include "field.hpp"
#include "initialization.hpp"

namespace cr_momentum {

    enum class Boundary {
        ZeroFlux,
        Outflow // Adiabatic outflow only; no inflow, zero diffusive flux.
    };

    // Evolves the cell-average q = dn/dRg at fixed spatial positions:
    // dq/dt = -dJ/dRg,
    // J = -(div(v)/3)*Rg*q - Rg^2*Dpp*d(q/Rg^2)/dRg.
    // Units: kpc, yr, GV; Dpp denotes D_RgRg in GV^2/yr.
    // Upwind adiabatic flux, two-point diffusion, Crank-Nicolson in time.
    // CN does not guarantee nonnegative densities for arbitrary timesteps.
    // Face coefficients are spatially uniform and computed from the model at construction.
    // Set vA=0 in the model to disable rigidity diffusion.
    // Grid, wind, dt and coefficients are fixed snapshots. Reconstruct when
    // they change. Reuses line storage; an instance is for sequential calls.
    class Solver {
    public:
        Solver(const SpatialGrid &grid, const RgGrid &Rg, const Field2D &vR, const Field2D &vZ, const cr_diffusion::DiffusionModel &model, double dt, Boundary lower = Boundary::ZeroFlux, Boundary upper = Boundary::ZeroFlux);
        // Updates physical spatial cells only; spatial ghosts are untouched.
        void advance(Field3D &q);

    private:
        int nR_, nz_, n_;
        double dt_;
        Boundary lower_boundary_, upper_boundary_;
        std::vector<double> divergence_; // [R][z], in yr^-1.
        std::vector<double> faces_, widths_; // Rigidity faces and linear bin widths, GV.
        std::vector<double> diff_left_, diff_right_; // Interior diffusive face-flux coefficients.
        struct Workspace {
            explicit Workspace(int n) : a(n), b(n), c(n), rhs(n) {}
            std::vector<double> a, b, c, rhs;
        };
        std::vector<Workspace> work_;
    };

} // namespace cr_momentum
