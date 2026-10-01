#include <iostream>
#include "CR_transport_solver.hpp" 
#include "field.hpp"
#include "initialization.hpp"
#include <numbers>

#include <cmath>
#include <array>
#include <omp.h>

int main() {
    const Grid2D grid(
        AxisGrid::geometric_capped(0.0, 200.0, 0.05, 1.0, 1.03),
        AxisGrid::geometric_capped(0.0, 200.0, 0.01, 1.0, 1.03));

    Field2D ndis_C(grid);
    Field2D ndis_H(grid);
    Field2D vR(grid, Field2D::Location::RadialFace);
    Field2D vZ(grid, Field2D::Location::VerticalFace);

    initialize_gas(ndis_H, 0.001, grid);
    initialize_wind_velocity(vR, vZ, grid); 
    initialize_CR_source(ndis_C, grid);
    solve_advection_equation(ndis_C, ndis_H, vR, vZ, grid);

    return 0;
}
