#pragma once  // 现代写法，防止头文件被重复包含（等价于传统的 #ifndef/#define/#endif）

#include "field.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <span>

// initialize constants for solver
#define DT 1.0 // Time step in years
#define NT 500 // Number of time steps


void initialize_CR_source(Field2D& ndis, const Grid2D& grid);
// vR must be RadialFace and vZ must be VerticalFace.
void initialize_wind_velocity(Field2D& vR, Field2D& vZ, const Grid2D& grid);
void initialize_gas(Field2D& ndis_H, const double nH, const Grid2D& grid);