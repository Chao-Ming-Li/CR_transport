#include "momentum_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cr_momentum {

    Solver::Solver(const SpatialGrid &grid, const RgGrid &Rg, const Field2D &vR, const Field2D &vZ, std::span<const double> Dpp, double dt, Boundary lower, Boundary upper) : nR_(grid.nR()), nz_(grid.nz()), n_(Rg.size()), dt_(dt), lower_boundary_(lower), upper_boundary_(upper), divergence_(static_cast<std::size_t>(nR_) * nz_), faces_(n_ + 1), widths_(n_), diff_left_(n_ + 1, 0.0), diff_right_(n_ + 1, 0.0), a_(n_), b_(n_), c_(n_), rhs_(n_) {
        if (!vR.matches(grid, Field2D::Location::RadialFace) || !vZ.matches(grid, Field2D::Location::VerticalFace))
            throw std::invalid_argument("Momentum wind fields must match spatial faces");
        if (!std::isfinite(dt) || dt < 0.0 || Dpp.size() != faces_.size())
            throw std::invalid_argument("Invalid momentum timestep or diffusion size");
        for (int k = 0; k <= n_; ++k) {
            faces_[k] = Rg.face(k);
        }
        for (int k = 0; k < n_; ++k)
            widths_[k] = Rg.width(k);
        for (int k = 1; k < n_; ++k) {
            const double conductance = faces_[k] * faces_[k] * Dpp[k] / Rg.center_distance(k - 1);
            diff_left_[k] = conductance / (Rg.center(k - 1) * Rg.center(k - 1));
            diff_right_[k] = conductance / (Rg.center(k) * Rg.center(k));
        }
        // Conservative cylindrical divergence; use annular volume, not the
        // volume-centroid coordinate in the radial denominator.
        for (int i = 0; i < nR_; ++i)
            for (int j = 0; j < nz_; ++j) {
                const double div =
                    (grid.radial_face_area(i + 1, j) * vR(i + 1, j) -
                     grid.radial_face_area(i, j) * vR(i, j) +
                     grid.vertical_face_area(i) * (vZ(i, j + 1) - vZ(i, j))) /
                    grid.volume(i, j);
                divergence_[static_cast<std::size_t>(i) * nz_ + j] = div;
            }
    }

    void Solver::advance(Field3D &q) {
        const double half_dt = 0.5 * dt_;
        // Caller supplies a density matching the construction grids.
        for (int i = 0; i < nR_; ++i)
            for (int j = 0; j < nz_; ++j) {
                const double rate = -divergence_[static_cast<std::size_t>(i) * nz_ + j] / 3.0;
                std::fill(a_.begin(), a_.end(), 0.0);
                std::fill(b_.begin(), b_.end(), 1.0);
                std::fill(c_.begin(), c_.end(), 0.0);
                for (int k = 0; k < n_; ++k)
                    rhs_[k] = q(k, i, j);
                // Shared face flux J_k = left*q_(k-1) + right*q_k.
                // CN: (I - dt*L/2) q_new = (I + dt*L/2) q_old.
                // Use each shared flux with opposite signs in adjoining cells.
                for (int k = 1; k < n_; ++k) {
                    const double velocity = rate * faces_[k];
                    const double left = std::max(velocity, 0.0) + diff_left_[k];
                    const double right = std::min(velocity, 0.0) - diff_right_[k];
                    const double old_flux = left * q(k - 1, i, j) + right * q(k, i, j);
                    b_[k - 1] += half_dt * left / widths_[k - 1];
                    c_[k - 1] += half_dt * right / widths_[k - 1];
                    a_[k] -= half_dt * left / widths_[k];
                    b_[k] -= half_dt * right / widths_[k];
                    rhs_[k - 1] -= half_dt * old_flux / widths_[k - 1];
                    rhs_[k] += half_dt * old_flux / widths_[k];
                }
                // Boundary contribution to Lq = -(J_right-J_left)/width.
                // Only the outgoing adiabatic velocity contributes; diffusion
                // and incoming flux remain zero. Apply the same CN rule at
                // both ends: b -= dt*L/2, rhs += dt*L*q_old/2.
                const auto add_outflow = [&](int k, double outward_velocity) {
                    const double loss = half_dt * outward_velocity / widths_[k];
                    b_[k] += loss;
                    rhs_[k] -= loss * q(k, i, j);
                };
                if (lower_boundary_ == Boundary::Outflow)
                    add_outflow(0, -std::min(rate * faces_[0], 0.0));
                if (upper_boundary_ == Boundary::Outflow)
                    add_outflow(n_ - 1, std::max(rate * faces_[n_], 0.0));

                // Thomas elimination and back substitution, no allocations.
                for (int k = 1; k < n_; ++k) {
                    const double factor = a_[k] / b_[k - 1];
                    b_[k] -= factor * c_[k - 1];
                    rhs_[k] -= factor * rhs_[k - 1];
                }
                rhs_[n_ - 1] /= b_[n_ - 1];
                for (int k = n_ - 2; k >= 0; --k)
                    rhs_[k] = (rhs_[k] - c_[k] * rhs_[k + 1]) / b_[k];
                for (int k = 0; k < n_; ++k)
                    q(k, i, j) = rhs_[k];
            }
    }

} // namespace cr_momentum
