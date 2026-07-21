#include <iostream>
#include "CR_transport.hpp"  
#include <cmath>
#include <array>
#include <omp.h>

int main() {
    // omp_set_num_threads(10);

    initialize_grids_log();

    Field2D ndis_C(NR, NZ);
    Field2D ndis_B(NR, NZ);
    Field2D ndis_H(NR, NZ);
    Field2D vR(NR, NZ);
    Field2D vZ(NR, NZ);
    constexpr std::array<double, 2> E = {1, 5};
    initialize_H(ndis_H, 0.001);
    initialize_wind_velocity(vR, vZ); // Assuming ndis_B is used for wind velocity in this context

    double D = 3.3e-8; // diffusion coefficient in kpc^2/yr
    initialize_disk_source(ndis_C);
    initialize_to_0(ndis_B);
    solve_diffusion_equation(ndis_C, ndis_B, ndis_H, D, 1.0); // Rg = 1.0 GV
    write_array_to_bin("ndis_C_1GeV_2nd_order.bin", ndis_C, NR * NZ);

    return 0;
}