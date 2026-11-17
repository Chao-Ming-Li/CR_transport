#pragma once

#include "field.hpp"

void initialize_CR_source(Field3D &q, const SpatialGrid &grid, const RgGrid &Rg, double source_R, double source_z, double source_width, double source_index);
// vR must be RadialFace and vZ must be VerticalFace.
void initialize_wind_velocity(Field2D &vR, Field2D &vZ, double wind_km_per_s);
void initialize_gas(Field2D &ndis_H, const double nH, const SpatialGrid &grid);

namespace cr_diffusion {

    // Dxx = D0 * (Rg/GV)^delta, without a beta factor.
    // Quasilinear reacceleration with turbulence factor w=1:
    // D_RgRg * Dxx = 4 * Rg^2 * vA^2 / [3*delta*(4-delta^2)*(4-delta)].
    // Reference: https://arxiv.org/abs/1108.1023, equation (2).
    // For a fixed charge, transforming p to rigidity cancels the charge factors.
    class DiffusionModel {
    public:
        // D0: kpc^2/yr at 1 GV; vA: km/s. Use 0 < delta < 2.
        DiffusionModel(double D0, double delta, double vA);
        double spatial_coefficient(double Rg) const;
        double rigidity_coefficient(double Rg) const;
        std::vector<double> rigidity_face_coefficients(const RgGrid &Rg) const;

    private:
        double D0_, delta_, reacceleration_;
    };

} // namespace cr_diffusion
