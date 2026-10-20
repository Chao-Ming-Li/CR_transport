#include "field.hpp"
#include "initialization.hpp"
#include <cmath>
#include <stdexcept>

namespace {
    void validate_grid_size(const Field2D &field, const Grid2D &grid,
                            Field2D::Location location = Field2D::Location::Centroid) {
        if (!field.matches(grid, location)) {
            throw std::invalid_argument("Field dimensions and location must match the grid");
        }
    }
} // namespace

void initialize_CR_source(Field2D &ndis, const Grid2D &grid) {
    validate_grid_size(ndis, grid);
    // Gaussian centered at (R, z) = (10, 10) kpc, as in the original model.
    // Instantaneous injection at t=0; sigma_R = sigma_z = 3/sqrt(2) kpc.
    // Density amplitude is one, independent of the evolution timestep.
    // Ghost cells are left unchanged.
    for (int i = 0; i < grid.nR(); ++i) {
        const double R1 = grid.radial_centroid(i) - 10.0;
        for (int j = 0; j < grid.nz(); ++j) {
            const double Z1 = grid.z().center(j) - 10.0;
            ndis(i, j) = std::exp(-R1 * R1 / 9.0 - Z1 * Z1 / 9.0);
        }
    }
}

void initialize_wind_velocity(Field2D &vR, Field2D &vZ, const Grid2D &grid) {
    validate_grid_size(vR, grid, Field2D::Location::RadialFace);
    validate_grid_size(vZ, grid, Field2D::Location::VerticalFace);
    // Constant +z wind: convert 300 km/s to kpc/yr. The lower face
    // has zero flux to match the reflecting z=0 boundary.
    constexpr double v_z = 300.0 * 1e5 * 31557600.0 / 3.0856775814913673e21;
    for (int i = 0; i < vR.nR(); ++i)
        for (int j = 0; j < vR.nz(); ++j)
            vR(i, j) = 0.0;
    for (int i = 0; i < vZ.nR(); ++i)
        for (int j = 0; j < vZ.nz(); ++j)
            vZ(i, j) = j == 0 ? 0.0 : v_z;
}

void initialize_gas(Field2D &ndis_H, const double nH, const Grid2D &grid) {
    validate_grid_size(ndis_H, grid);
    // Uniform cell-centered gas density, independent of position.
    for (int i = 0; i < grid.nR(); ++i) {
        for (int j = 0; j < grid.nz(); ++j) {
            ndis_H(i, j) = nH;
        }
    }
}
