#include "cross_section.hpp"
#include "initialization.hpp"
#include <cmath>

void initialize_CR_source(Field3D &q, const SpatialGrid &grid, const RgGrid &Rg, double source_R, double source_z, double source_width, double source_index) {
    for (int k = 0; k < q.nrigidity(); ++k) {
        // Evaluate the example power-law source spectrum dn/dRg ~ Rg^-source_index at bin centers.
        const double dN_dRg = std::pow(Rg.center(k), -source_index);
        for (int i = 0; i < grid.nR(); ++i) {
            const double R1 = grid.radial_centroid(i) - source_R;
            for (int j = 0; j < grid.nz(); ++j) {
                const double Z1 = grid.z().center(j) - source_z;
                q(k, i, j) = std::exp(-R1 * R1 / (source_width * source_width) - Z1 * Z1 / (source_width * source_width)) * dN_dRg;
            }
        }
    }
}

void initialize_wind_velocity(Field2D &vR, Field2D &vZ, double wind_km_per_s) {
    // Constant +z wind: convert km/s to kpc/yr. The lower face
    // has zero flux to match the reflecting z=0 boundary.
    // Thus div(v)=v_z/dz only in the bottom cell and is zero above it.
    const double v_z = wind_km_per_s * physical_constants::km_s_to_kpc_yr;
    for (int i = 0; i < vR.nR(); ++i) {
        for (int j = 0; j < vR.nz(); ++j) {
            vR(i, j) = 0.0;
        }
    }
    for (int i = 0; i < vZ.nR(); ++i) {
        for (int j = 0; j < vZ.nz(); ++j) {
            vZ(i, j) = j == 0 ? 0.0 : v_z;
        }
    }
}

void initialize_gas(Field2D &ndis_H, const double nH, const SpatialGrid &grid) {
    // Uniform cell-centered gas density, independent of position.
    for (int i = 0; i < grid.nR(); ++i) {
        for (int j = 0; j < grid.nz(); ++j) {
            ndis_H(i, j) = nH;
        }
    }
}

namespace cr_diffusion {

    DiffusionModel::DiffusionModel(double D0, double delta, double vA) : D0_(D0), delta_(delta) {
        const double speed = vA * physical_constants::km_s_to_kpc_yr;
        reacceleration_ = 4.0 * speed * speed / (3.0 * delta * (4.0 - delta * delta) * (4.0 - delta));
    }

    double DiffusionModel::spatial_coefficient(double Rg) const {
        return D0_ * std::pow(Rg, delta_);
    }

    double DiffusionModel::rigidity_coefficient(double Rg) const {
        return reacceleration_ * Rg * Rg / spatial_coefficient(Rg);
    }

    std::vector<double> DiffusionModel::rigidity_face_coefficients(const RgGrid &Rg) const {
        std::vector<double> coefficients(Rg.size() + 1);
        for (int k = 0; k <= Rg.size(); ++k) {
            coefficients[k] = rigidity_coefficient(Rg.face(k));
        }
        return coefficients;
    }

} // namespace cr_diffusion
