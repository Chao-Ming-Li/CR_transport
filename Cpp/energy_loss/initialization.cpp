#include "initialization.hpp"
#include <cmath>

void initialize_CR_source(Field3D &q, const SpatialGrid &grid, const RgGrid &Rg) {

    for (int k = 0; k < q.nrigidity(); ++k) {
        const double dN_dRg = 1.0  / Rg.center(k) / Rg.center(k); // dn/dRg ~ Rg^-2, a simple power-law spectrum.
        for (int i = 0; i < grid.nR(); ++i) {
            const double R1 = grid.radial_centroid(i) - 10.0;
            for (int j = 0; j < grid.nz(); ++j) {
                const double Z1 = grid.z().center(j) - 10.0;
                q(k, i, j) = std::exp(-R1 * R1 / 9.0 - Z1 * Z1 / 9.0) * dN_dRg;
            }
        }
    }
}

void initialize_wind_velocity(Field2D &vR, Field2D &vZ) {

    // Constant +z wind: convert 300 km/s to kpc/yr. The lower face
    // has zero flux to match the reflecting z=0 boundary.
    constexpr double v_z = 300.0 * 1e5 * 31557600.0 / 3.0856775814913673e21;
    for (int i = 1; i <= vR.nR(); ++i)
        for (int j = 0; j < vR.nz(); ++j)
            vR(i, j) = 0.0;
    for (int i = 0; i < vZ.nR(); ++i)
        for (int j = 1; j <= vZ.nz(); ++j)
            vZ(i, j) = j == 0 ? 0.0 : v_z;
}

void initialize_gas(Field2D &ndis_H, const double nH, const SpatialGrid &grid) {
    // Uniform cell-centered gas density, independent of position.
    for (int i = 0; i < grid.nR(); ++i) {
        for (int j = 0; j < grid.nz(); ++j) {
            ndis_H(i, j) = nH;
        }
    }
}
