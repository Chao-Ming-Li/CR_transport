#include "cross_section.hpp"
#include "nuclear_solver.hpp"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cr_nuclear {

    Solver::Solver(const SpatialGrid &grid, const RgGrid &Rg, const Field2D &gas, const std::vector<cross_section::Isotope> &isotopes, double dt, double helium_enhancement)
        : nR_(grid.nR()), nz_(grid.nz()), nRg_(Rg.size()), count_(static_cast<int>(isotopes.size())), dt_(dt), gas_(gas),
          source_(grid, Rg), order_(count_), destruction_(count_ * nRg_), transfers_(count_) {

        // s,d index the selected fields; id,daughter_id index the cross-section table.
        // A selected position need not equal the corresponding table index.
        for (int s = 0; s < count_; ++s) {
            const int id = static_cast<int>(cross_section::index(isotopes[s]));
            order_[s] = s;

            const cross_section::Nuclide &parent = cross_section::nuclide(isotopes[s]);
            // Cache beta*c*sigma for each parent bin; multiply by nH during evolution.
            for (int k = 0; k < nRg_; ++k) {
                const double p_per_nucleon = Rg.center(k) * parent.charge / parent.mass_number;
                const double mass = physical_constants::nucleon_rest_energy_GeV;
                const double beta = p_per_nucleon / std::sqrt(p_per_nucleon * p_per_nucleon + mass * mass);
                // Effective speed in cm/yr includes the optional helium enhancement.
                const double effective_speed = beta * physical_constants::speed_of_light_cm_s * physical_constants::seconds_per_year * helium_enhancement;
                destruction_[s * nRg_ + k] = effective_speed * cross_section::sigma[id][id];

                // Map each parent interval to daughter rigidity at fixed energy/nucleon.
                for (int d = 0; d < count_; ++d) {
                    const int daughter_id = static_cast<int>(cross_section::index(isotopes[d]));
                    if (id == daughter_id || cross_section::sigma[id][daughter_id] == 0.0)
                        continue;
                    const cross_section::Nuclide &daughter = cross_section::nuclide(isotopes[d]);
                    const double parent_A_over_Z = static_cast<double>(parent.mass_number) / parent.charge;
                    const double daughter_A_over_Z = static_cast<double>(daughter.mass_number) / daughter.charge;
                    const double alpha = daughter_A_over_Z / parent_A_over_Z;
                    const double mapped_lower = alpha * Rg.face(k);
                    const double mapped_upper = alpha * Rg.face(k + 1);
                    for (int l = 0; l < nRg_; ++l) {
                        const double overlap = std::min(mapped_upper, Rg.face(l + 1)) - std::max(mapped_lower, Rg.face(l));
                        if (overlap > 0.0) {
                            // q_d = q_p/alpha; overlap/width averages over the daughter bin.
                            // Only in-domain overlaps are retained, without renormalization.
                            const double rate = effective_speed * cross_section::sigma[id][daughter_id] * overlap / (alpha * Rg.width(l));
                            transfers_[d].push_back({s, k, l, rate});
                        }
                    }
                }
            }
        }

        // sigma is upper triangular in table order. This also accepts a caller's
        // arbitrary subset/order without changing the ordering of output fields.
        for (int s = 0; s < count_; ++s) {
            for (int t = s + 1; t < count_; ++t) {
                if (cross_section::index(isotopes[order_[t]]) < cross_section::index(isotopes[order_[s]]))
                    std::swap(order_[s], order_[t]);
            }
        }
    }

    void Solver::add_production(int daughter, const std::vector<Field3D> &density, Field3D &source) const {
        std::fill(source.data().begin(), source.data().end(), 0.0);
        for (const Transfer &transfer : transfers_[daughter]) {
            for (int i = 0; i < nR_; ++i) {
                for (int j = 0; j < nz_; ++j) {
                    source(transfer.daughter_bin, i, j) += gas_(i, j) * transfer.rate * density[transfer.parent](transfer.parent_bin, i, j);
                }
            }
        }
    }

    void Solver::calculate_terms(const std::vector<Field3D> &density, std::vector<Field3D> &production, std::vector<Field3D> &loss) const {
        for (int s = 0; s < count_; ++s) {
            add_production(s, density, production[s]);
            std::fill(loss[s].data().begin(), loss[s].data().end(), 0.0);
            for (int k = 0; k < nRg_; ++k) {
                for (int i = 0; i < nR_; ++i) {
                    for (int j = 0; j < nz_; ++j) {
                        loss[s](k, i, j) = gas_(i, j) * destruction_[s * nRg_ + k] * density[s](k, i, j);
                    }
                }
            }
        }
    }

    void Solver::advance(std::vector<Field3D> &density) {
        if (dt_ == 0.0)
            return;
        // Each radial row advances its complete parent-to-daughter chain.
        // Neither the reaction nor the rigidity remap couples spatial cells.
#ifdef _OPENMP
        const int threads = std::min(omp_get_max_threads(), nR_);
        #pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
        for (int i = 0; i < nR_; ++i) {
            for (const int s : order_) {
                for (int k = 0; k < nRg_; ++k) {
                    for (int j = 0; j < nz_; ++j) {
                        source_(k, i, j) = 0.0;
                    }
                }
                for (const Transfer &transfer : transfers_[s]) {
                    for (int j = 0; j < nz_; ++j) {
                        source_(transfer.daughter_bin, i, j) += gas_(i, j) * transfer.rate * density[transfer.parent](transfer.parent_bin, i, j);
                    }
                }
                for (int k = 0; k < nRg_; ++k) {
                    for (int j = 0; j < nz_; ++j) {
                        const double rate = gas_(i, j) * destruction_[s * nRg_ + k];
                        density[s](k, i, j) = (density[s](k, i, j) + dt_ * source_(k, i, j)) / (1.0 + dt_ * rate);
                    }
                }
            }
        }
    }

} // namespace cr_nuclear
