#include "momentum_solver.hpp"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cr_momentum {

    Solver::Solver(const SpatialGrid &grid, const RgGrid &Rg, const Field2D &vR, const Field2D &vZ, const cr_diffusion::DiffusionModel &model, double dt, Boundary lower, Boundary upper)
        : nR_(grid.nR()), nz_(grid.nz()), n_(Rg.size()), dt_(dt), lower_boundary_(lower), upper_boundary_(upper),
          divergence_(static_cast<std::size_t>(nR_) * nz_), faces_(n_ + 1), widths_(n_), diff_left_(n_ + 1, 0.0), diff_right_(n_ + 1, 0.0) {
        work_.emplace_back(n_);
        for (int k = 0; k <= n_; ++k) {
            faces_[k] = Rg.face(k);
        }
        for (int k = 0; k < n_; ++k) {
            widths_[k] = Rg.width(k);
        }

        // Interior diffusive flux: J_k = diff_left[k]*q[k-1] - diff_right[k]*q[k].
        for (int k = 1; k < n_; ++k) {
            const double conductance = faces_[k] * faces_[k] * model.rigidity_coefficient(faces_[k]) / Rg.center_distance(k - 1);
            diff_left_[k] = conductance / (Rg.center(k - 1) * Rg.center(k - 1));
            diff_right_[k] = conductance / (Rg.center(k) * Rg.center(k));
        }
        // Conservative cylindrical divergence; use annular volume, not the
        // volume-centroid coordinate in the radial denominator.
        for (int i = 0; i < nR_; ++i) {
            for (int j = 0; j < nz_; ++j) {
                const double div = (
                    grid.radial_face_area(i + 1, j) * vR(i + 1, j)
                    - grid.radial_face_area(i, j) * vR(i, j)
                    + grid.vertical_face_area(i) * (vZ(i, j + 1) - vZ(i, j))
                ) / grid.volume(i, j);
                divergence_[static_cast<std::size_t>(i) * nz_ + j] = div;
            }
        }
    }

    void Solver::advance(Field3D &q) {
        const double half_dt = 0.5 * dt_;
        // Caller supplies a density matching the construction grids.
#ifdef _OPENMP
        const int threads = std::min(omp_get_max_threads(), nR_);
        while (work_.size() < static_cast<std::size_t>(threads)) {
            work_.emplace_back(n_);
        }
        #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
        for (int i = 0; i < nR_; ++i) {
#ifdef _OPENMP
            Workspace &work = work_[omp_get_thread_num()];
#else
            Workspace &work = work_.front();
#endif
            auto &a = work.a;
            auto &b = work.b;
            auto &c = work.c;
            auto &rhs = work.rhs;
            for (int j = 0; j < nz_; ++j) {
                // dRg/dt = adiabatic_rate*Rg: expansion cools, compression heats.
                const double adiabatic_rate = -divergence_[static_cast<std::size_t>(i) * nz_ + j] / 3.0;
                // Build the tridiagonal CN system for this spatial cell.
                std::fill(a.begin(), a.end(), 0.0);
                std::fill(b.begin(), b.end(), 1.0);
                std::fill(c.begin(), c.end(), 0.0);
                for (int k = 0; k < n_; ++k) {
                    rhs[k] = q(k, i, j);
                }

                // Shared face flux J_k = left*q_(k-1) + right*q_k.
                // CN: (I - dt*L/2) q_new = (I + dt*L/2) q_old.
                // Use each shared flux with opposite signs in adjoining cells.
                for (int k = 1; k < n_; ++k) {
                    const double velocity = adiabatic_rate * faces_[k];
                    const double left = std::max(velocity, 0.0) + diff_left_[k];
                    const double right = std::min(velocity, 0.0) - diff_right_[k];
                    const double old_flux = left * q(k - 1, i, j) + right * q(k, i, j);
                    b[k - 1] += half_dt * left / widths_[k - 1];
                    c[k - 1] += half_dt * right / widths_[k - 1];
                    a[k] -= half_dt * left / widths_[k];
                    b[k] -= half_dt * right / widths_[k];
                    rhs[k - 1] -= half_dt * old_flux / widths_[k - 1];
                    rhs[k] += half_dt * old_flux / widths_[k];
                }
                // Boundary contribution to Lq = -(J_right-J_left)/width.
                // Only the outgoing adiabatic velocity contributes; diffusion
                // and incoming flux remain zero. Apply the same CN rule at
                // both ends: b -= dt*L/2, rhs += dt*L*q_old/2.
                if (lower_boundary_ == Boundary::Outflow) {
                    const double loss = half_dt * (-std::min(adiabatic_rate * faces_[0], 0.0)) / widths_[0];
                    b[0] += loss;
                    rhs[0] -= loss * q(0, i, j);
                }
                if (upper_boundary_ == Boundary::Outflow) {
                    const double loss = half_dt * std::max(adiabatic_rate * faces_[n_], 0.0) / widths_[n_ - 1];
                    b[n_ - 1] += loss;
                    rhs[n_ - 1] -= loss * q(n_ - 1, i, j);
                }

                // Thomas elimination and back substitution, no allocations.
                for (int k = 1; k < n_; ++k) {
                    const double factor = a[k] / b[k - 1];
                    b[k] -= factor * c[k - 1];
                    rhs[k] -= factor * rhs[k - 1];
                }
                rhs[n_ - 1] /= b[n_ - 1];
                for (int k = n_ - 2; k >= 0; --k) {
                    rhs[k] = (rhs[k] - c[k] * rhs[k + 1]) / b[k];
                }
                for (int k = 0; k < n_; ++k) {
                    q(k, i, j) = rhs[k];
                }
            }
        }
    }

} // namespace cr_momentum
