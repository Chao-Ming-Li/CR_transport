#pragma once  // 现代写法，防止头文件被重复包含（等价于传统的 #ifndef/#define/#endif）

#include "field.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <span>



// Distribution initializers implemented in initialization.cpp.
// Field dimensions must match the supplied grid.


void apply_boundary_conditions(Field2D& ndis, const Grid2D& grid);
void write_array_to_bin(const std::string& filename, Field2D& arr, const std::size_t size);


// Conservative finite-volume advection of nonnegative cell averages.
// ndis/temp are distinct Center fields; vR is RadialFace; vZ is VerticalFace.
// Velocities are stationary over dT and are not modified. dT must be finite >= 0.
// Uses limited quadratic reconstruction and SSPRK2, with automatic substeps
// at outgoing CFL <= 0.45 (summed over active directions). Lower boundaries
// reflect; upper boundaries allow outflow with vacuum inflow.
// Historical TVD names are retained, but the smooth-curvature exception
// preserves resolved extrema and does not guarantee strict TVD for all data.
// Boundary cells use constant reconstruction. temp is stage scratch storage.
void advection_TVD(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT, const Grid2D& grid);

void advection_TVD_R(Field2D& ndis, Field2D& temp, Field2D& vR, double dT, const Grid2D& grid);
void advection_TVD_Z(Field2D& ndis, Field2D& temp, Field2D& vZ, double dT, const Grid2D& grid);


// inline double minmod(double a,double b);

// inline double MC_limiter(double sL,double sC,double sR);

// void advection_PLM(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT, const Grid2D& grid);

void solve_advection_equation(Field2D& ndis_C, Field2D& ndis_H, Field2D& vR, Field2D& vZ, const Grid2D& grid);