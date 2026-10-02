#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "diffusion_func.h"

int main(void) {

    // Set number of threads to the number of CPUs
    int ncore = omp_get_num_procs(); 
    omp_set_num_threads(ncore);

    // Dynamically allocate memory for ndis and temp
    double *ndis_O16 = malloc(NR * NZ * sizeof(double));
    double *ndis_N15 = malloc(NR * NZ * sizeof(double));
    double *ndis_N14 = malloc(NR * NZ * sizeof(double));
    double *ndis_C13 = malloc(NR * NZ * sizeof(double));
    double *ndis_C12 = malloc(NR * NZ * sizeof(double));
    double *ndis_B11 = malloc(NR * NZ * sizeof(double));
    double *ndis_B10 = malloc(NR * NZ * sizeof(double));
    double *ndis_Be10 = malloc(NR * NZ * sizeof(double));
    double *ndis_Be9 = malloc(NR * NZ * sizeof(double));
    double *ndis_H = malloc(NR * NZ * sizeof(double));
    double *D = malloc(NR * NZ * sizeof(double));

    double Rg[6] = {10., 16., 25., 40., 63., 100.}; // Rigidity
    double delta[6] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; // diffusion index
    double D0[8] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4}; // log10 of diffusion coefficient
    double dT, ss_time = 0.0, gas_type = -1.0, Dd, Dh;
    char grid_type[10] = "log";
    initialize_grids_log();
    initialize_H(ndis_H, gas_type);

    for (int i = 2; i < 6; i++) {
        for (int d = 0; d < 4; d++){
            Dd = pow(10.0, 28.0 + D0[d]) * pow(Rg[0], delta[i]) / 3.086e21 / 3.086e21 * 3.15576e7; // in kpc^2 / yr;
            dT = 1000.0;
            // Initialize particle density and apply boundary conditions, C and O are primary CRs, N and B are secondary CRs
            initialize_to_disk(ndis_O16, 1.0, dT); // O16 are all primary
            initialize_to_disk(ndis_N14, 0.1, dT); 
            initialize_to_disk(ndis_C12, 1.0, dT); 
            initialize_to_0(ndis_N15);
            initialize_to_0(ndis_C13);
            initialize_to_0(ndis_B11);
            initialize_to_0(ndis_B10);
            initialize_to_0(ndis_Be10);
            initialize_to_0(ndis_Be9);
            // initialize_D(D, Dd, Dh);
            solve_diffusion_equation_CN(ndis_O16, ndis_N15, ndis_N14, ndis_C13, ndis_C12, ndis_B11, ndis_B10, ndis_Be10, ndis_Be9, ndis_H, Dd, dT, Rg[0], &ss_time);
            // output as bin file
            char fn0[100], fn1[100], fn2[100], fn3[100], fn4[100], fn5[100], fn6[100], fn7[100], fn8[100], namebody[80];
            snprintf(namebody, sizeof(namebody), "%dGV_R%d_Z%d_Htot_t1e9_Rinj%d_D%.1f_delta%.1f.bin",(int)Rg[0], (int)round(R[NR-1]), (int)round(Z[NZ-1]),(int)R_inj, D0[d], delta[i]);
            snprintf(fn0, sizeof(fn0), "ndis_O16_%s", namebody);
            snprintf(fn1, sizeof(fn1), "ndis_N15_%s", namebody);
            snprintf(fn2, sizeof(fn2), "ndis_N14_%s", namebody);
            snprintf(fn3, sizeof(fn3), "ndis_C13_%s", namebody);
            snprintf(fn4, sizeof(fn4), "ndis_C12_%s", namebody);
            snprintf(fn5, sizeof(fn5), "ndis_B11_%s", namebody);
            snprintf(fn6, sizeof(fn6), "ndis_B10_%s", namebody);
            snprintf(fn7, sizeof(fn7), "ndis_Be10_%s", namebody);
            snprintf(fn8, sizeof(fn8), "ndis_Be9_%s", namebody);
            write_array_to_bin(fn0, ndis_O16, NR*NZ);
            write_array_to_bin(fn1, ndis_N15, NR*NZ);
            write_array_to_bin(fn2, ndis_N14, NR*NZ);           
            write_array_to_bin(fn3, ndis_C13, NR*NZ);           
            write_array_to_bin(fn4, ndis_C12, NR*NZ);           
            write_array_to_bin(fn5, ndis_B11, NR*NZ);           
            write_array_to_bin(fn6, ndis_B10, NR*NZ);           
            write_array_to_bin(fn7, ndis_Be10, NR*NZ);           
            write_array_to_bin(fn8, ndis_Be9, NR*NZ); 
        } 
    }

    // Free the dynamically allocated memory
    free(ndis_O16);
    free(ndis_N15);
    free(ndis_N14);
    free(ndis_C13);
    free(ndis_C12);
    free(ndis_B11);
    free(ndis_B10);
    free(ndis_Be10);
    free(ndis_Be9);
    free(ndis_H);
    free(D);
    return 0;
}



