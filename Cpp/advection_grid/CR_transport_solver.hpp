#pragma once  
#include "field.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <span>



// Distribution initializers implemented in initialization.cpp.
// Field dimensions must match the supplied grid.

void apply_boundary_conditions(Field2D& ndis, const Grid2D& grid);
void write_array_to_bin(const std::string& filename, Field2D& arr, const std::size_t size);
void solve_advection_equation(Field2D& ndis_C, Field2D& ndis_H, const Field2D& vR, const Field2D& vZ, const Grid2D& grid);
