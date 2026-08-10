#include <iostream>
#include "CR_advection.hpp" 
#include <cmath>
#include <array>
#include <omp.h>

int main() {
    // omp_set_num_threads(10);

    initialize_grids_log_inv();

    // Field2D ndis_C;
    // Field2D ndis_B;
    // Field2D ndis_H;
    // Field2D vR;
    // Field2D vZ;
    // constexpr std::array<double, 2> E = {1, 5};
    // initialize_H(ndis_H, 0.001);
    // initialize_wind_velocity(vR, vZ); 
    // write_array_to_bin("vR.bin", vR, (NR + 4) * (NZ + 4));
    // write_array_to_bin("vZ.bin", vZ, (NR + 4) * (NZ + 4));

    // double Rg = 1.0; // rigidity GV
    // double D = 3.3e-8; // diffusion coefficient in kpc^2/yr
    // initialize_disk_source(ndis_C);

    // solve_advection_equation(ndis_C, ndis_B, ndis_H, vR, vZ);
    // write_array_to_bin("ndis_C_linear.bin", ndis_C, (NR + 4) * (NZ + 4));
    // write_array_to_bin("ndis_B_1GeV_vbreeze1000.bin", ndis_B, (NR + 4) * (NZ + 4));    

    // solve_transport_equation(ndis_C, ndis_B, vR, vZ, ndis_H, D, Rg);
    // write_array_to_bin("ndis_C_1GeV_D1e28_vbreeze1000.bin", ndis_C, NR * NZ);
    // write_array_to_bin("ndis_B_1GeV_D1e28_vbreeze1000.bin", ndis_B, NR * NZ);
    // write_array_to_bin("ndis_C_dt10000_opt2.bin", ndis_C, (NR + 4) * (NZ + 4));
    // write_array_to_bin("ndis_B_dt10000_opt2.bin", ndis_B, (NR + 4) * (NZ + 4));

    return 0;
}