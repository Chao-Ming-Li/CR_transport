#pragma once // 现代写法，防止头文件被重复包含（等价于传统的 #ifndef/#define/#endif）

#include "field.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <span>

// Distribution initializers implemented in initialization.cpp.
// Field dimensions must match the supplied grid.

void apply_boundary_conditions(Field2D &ndis, const Grid2D &grid);
void write_array_to_bin(const std::string &filename, Field2D &arr, const std::size_t size);

void advection_upwind(Field2D &ndis, Field2D &temp, const Field2D &vR, const Field2D &vZ,
                      const double dT, const Grid2D &grid);
void advection_upwind_1(Field2D &ndis, Field2D &temp, const Field2D &vR, const Field2D &vZ,
                        const double dT, const Grid2D &grid);
void advection_TVD(Field2D &ndis, Field2D &temp, const Field2D &vR, const Field2D &vZ,
                   const double dT, const Grid2D &grid);

void advection_TVD_R(Field2D &ndis, Field2D &temp, const Field2D &vR, const double dT,
                     const Grid2D &grid);
void advection_TVD_Z(Field2D &ndis, Field2D &temp, const Field2D &vZ, const double dT,
                     const Grid2D &grid);
void advance_TVD_SSPRK2(Field2D &ndis, Field2D &temp, const Field2D &vR, const Field2D &vZ,
                        const double dT, const Grid2D &grid);

// MC-limited linear reconstruction; R then Z forward-Euler sweeps.
// Uses cached radial centroids and physical annular volumes. ndis/temp must
// be distinct Centroid fields. Caller supplies a stable dT and zero velocities
// at reflecting lower faces. This implementation is in CR_advection.cpp only.
void advection_PLM(Field2D &ndis, Field2D &temp, const Field2D &vR, const Field2D &vZ,
                   const double dT, const Grid2D &grid);

void solve_advection_equation(Field2D &ndis_C, Field2D &ndis_H, const Field2D &vR,
                              const Field2D &vZ, const Grid2D &grid);
