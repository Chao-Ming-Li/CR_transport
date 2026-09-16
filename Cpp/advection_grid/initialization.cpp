#include "field.hpp"
#include "initialization.hpp"
#include <cmath>
#include <stdexcept>

namespace {
void validate_grid_size(const Field2D& field, const Grid2D& grid,
                         Field2D::Location location = Field2D::Location::Centroid)
{
    if (!field.matches(grid, location)) {
        throw std::invalid_argument("Field dimensions and location must match the grid");
    }
}
} // namespace

void initialize_CR_source(Field2D& ndis, const Grid2D& grid)
{
    validate_grid_size(ndis, grid);
    // Gaussian centered at (R, z) = (30, 30) kpc, as in the original model.
    // Ghost cells are left unchanged.
    for (int i = 0; i < grid.nR(); ++i) {
        const double R1 = grid.radial_centroid(i) - 30.0;
        for (int j = 0; j < grid.nz(); ++j) {
            const double Z1 = grid.z().center(j) - 30.0;
            ndis(i, j) = std::exp(-R1 * R1 / 9.0 - Z1 * Z1 / 9.0) * DT;
        }
    }
}

void initialize_wind_velocity(Field2D& vR, Field2D& vZ, const Grid2D& grid)
{
    validate_grid_size(vR, grid, Field2D::Location::RadialFace);
    validate_grid_size(vZ, grid, Field2D::Location::VerticalFace);
    // Radial faces: (R.face(i), z.center(j)).
    for (int i = 0; i < vR.nR(); ++i)
        for (int j = 0; j < vR.nz(); ++j)
            vR(i, j) = 0.0;

    // Vertical faces: (R.center(i), z.face(j)); velocity in kpc / yr.
    for (int i = 0; i < vZ.nR(); ++i)
        for (int j = 0; j < vZ.nz(); ++j)
            vZ(i, j) = j == 0 ? 0.0 : 0.1;

}

void initialize_gas(Field2D& ndis_H, const double nH, const Grid2D& grid)
{
    validate_grid_size(ndis_H, grid);
    // Uniform cell-centered gas density, independent of position.
    for (int i = 0; i < grid.nR(); ++i) {
        for (int j = 0; j < grid.nz(); ++j) {
            ndis_H(i, j) = nH;
        }
    }
}
