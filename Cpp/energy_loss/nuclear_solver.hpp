#pragma once

#include "field.hpp"
#include "cross_section.hpp"

namespace cr_nuclear {

    // One Field3D per selected isotope, in the same order as isotopes.
    // q=dn/dRg; gas is nH [cm^-3], dt [yr]. Cross sections come from sigma.
    // Straight-ahead approximation: kinetic energy per nucleon is unchanged,
    // using a common nucleon mass, hence Rg_d = (A_d/Z_d)/(A_p/Z_p)*Rg_p.
    // Piecewise-constant parent bins are conservatively mapped to daughter bins.
    // Products outside the rigidity domain escape; no renormalization is applied.
    // Backward Euler, parents before daughters: first-order time accuracy.
    // No radioactive decay, ionization, or proton reactions are included.
    class Solver {
    public:
        Solver(const SpatialGrid &grid, const RgGrid &Rg, const Field2D &gas, const std::vector<cross_section::Isotope> &isotopes, double dt, double helium_enhancement);
        void advance(std::vector<Field3D> &density);
        // Instantaneous production and destruction (both nonnegative for q>=0),
        // in q/yr. Output fields must be separate from density and each other.
        void calculate_terms(const std::vector<Field3D> &density, std::vector<Field3D> &production, std::vector<Field3D> &loss) const;

    private:
        struct Transfer {
            int parent; // Position in the selected isotope list.
            int parent_bin, daughter_bin;
            double rate; // cm^3/yr, including the conservative bin mapping.
        };
        int nR_, nz_, nRg_, count_;
        double dt_;
        Field2D gas_;
        Field3D source_;
        std::vector<int> order_; // Selected field positions sorted in parent-to-daughter order.
        std::vector<double> destruction_; // [species][rigidity], cm^3/yr.
        std::vector<std::vector<Transfer>> transfers_; // Grouped by daughter.
        void add_production(int daughter, const std::vector<Field3D> &density, Field3D &source) const;
    };

} // namespace cr_nuclear
