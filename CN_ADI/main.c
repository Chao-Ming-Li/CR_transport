#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion_func.h"
int main(void) {

    double Rg[6] = {10., 16., 25., 40., 63., 100.}; // Rigidity
    double delta[6] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; // diffusion index
    double D0[10] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}; // log10 of diffusion coefficient
    double Rd[10] = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0}; // disk radius
    double ZH[6] = {3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // halo height
    // double delta[6] = {0.4, 0.45, 0.5, 0.55, 0.6}; // diffusion index
    // double D0[10] = { 0.7,0.75, 0.8,0.85, 0.9}; // log10 of diffusion coefficient
    // double Rd[10] = {14.0, 15.0, 16.0,17.0, 18.0}; // disk radius
    // double ZH[6] = { 6.0, 6.5, 7.0, 7.5, 8.0}; // halo height
    double BC_ratio[6] = {0.2857, 0.2543, 0.219, 0.1834, 0.1589, 0.1247}; // B/C ratio at different rigidities
    double BC_ratio_err[6] = {0.0060531, 0.00529528, 0.00472969, 0.00430116, 0.00429535, 0.00410488};
    double BeC_ratio[6] = {0.0978, 0.088, 0.0769, 0.068, 0.0584, 0.0471}; // Be/C ratio at different rigidities
    double BeC_ratio_err[6] = {0.00233452, 0.00208806, 0.00189737, 0.00187883, 0.00180278, 0.0019105};
    initialize_grids_log();
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
    double *ndis_Be7 = malloc(NR * NZ * sizeof(double));
    double *ndis_H = malloc(NR * NZ * sizeof(double));
    double BC_ratio_M[6], BeC_ratio_M[6], Be10_Be9_ratio_M[6];
    double M_flux[60];

    double gas_type = -3.0, Dd;
    initialize_H(ndis_H, gas_type);

    char fn[100];
    snprintf(fn, sizeof(fn), "D%.2f_delta%.2f_ef1.3_gasNS.bin",D0[0], delta[0]);

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 10; j++) {
            double chi2 = 0.0;
            for (int k = 0; k < 6; k++) {
                Dd = pow(10.0, 28.0 + D0[0]) * pow(Rg[k], delta[0]) / 3.086e21 / 3.086e21 * 3.15576e7; // in kpc^2 / yr;
                // Initialize particle density and apply boundary conditions, C and O are primary CRs, N and B are secondary CRs
                initialize_to_disk(ndis_O16, 1.0, Rd[j]); // O16 are all primary
                initialize_to_disk(ndis_N14, 0.1, Rd[j]); 
                initialize_to_disk(ndis_C12, 1.0, Rd[j]); 
                initialize_to_0(ndis_N15);
                initialize_to_0(ndis_C13);
                initialize_to_0(ndis_B11);
                initialize_to_0(ndis_B10);
                initialize_to_0(ndis_Be10);
                initialize_to_0(ndis_Be9);
                initialize_to_0(ndis_Be7);
                // printf("Thread %d processing Rigidity: %f GV with Dd: %e\n", omp_get_thread_num(), Rg[k], Dd);        
                solve_diffusion_equation_CN(ndis_O16, ndis_N15, ndis_N14, ndis_C13, ndis_C12, ndis_B11, ndis_B10, ndis_Be10, ndis_Be9, ndis_Be7, ndis_H, Dd, Rg[k], Rd[j], ZH[i]);
                double *arrays[10] = {ndis_O16, ndis_N15, ndis_N14, ndis_C13, ndis_C12, ndis_B11, ndis_B10, ndis_Be10, ndis_Be9, ndis_Be7};
                for (int i = 0; i < 10; i++) {
                    M_flux[i * 6 + k] = arrays[i][80 * NZ]; // Store the flux at R=8 kpc and Z=0 for each species and rigidity
                }
                // BC_ratio_M[k] = (ndis_B11[80 * NZ] + ndis_B10[80 * NZ]) / (ndis_C13[80 * NZ] + ndis_C12[80 * NZ]);
                // BeC_ratio_M[k] = (ndis_Be10[80 * NZ] + ndis_Be9[80 * NZ] + ndis_Be7[80 * NZ]) / (ndis_C13[80 * NZ] + ndis_C12[80 * NZ]);
                // Be10_Be9_ratio_M[k] = ndis_Be10[80 * NZ] / ndis_Be9[80 * NZ];
                // double BC_diff_M = BC_ratio_M[k] - BC_ratio[k];
                // double BeC_diff_M = BeC_ratio_M[k] - BeC_ratio[k];
                // chi2 += BC_diff_M * BC_diff_M / (BC_ratio_err[k] * BC_ratio_err[k]) + BeC_diff_M * BeC_diff_M / (BeC_ratio_err[k] * BeC_ratio_err[k]);
                // char fn0[100], fn1[100], fn2[100], fn3[100], fn4[100], fn5[100], fn6[100], fn7[100], fn8[100], fn9[100], namebody[80];
                // snprintf(namebody, sizeof(namebody), "%dGV_R%d_Z%d_HMC_t1e8_Rinj%d_ZH%d_dt5e4_D%.2f_delta%.2f.bin", (int)Rg[k],(int)round(R[NR-1]), (int)round(Z[NZ-1]),(int)Rd[j],(int)ZH[i], D0[6], delta[2]);
                // snprintf(fn0, sizeof(fn0), "ndis_O16_%s", namebody);
                // snprintf(fn1, sizeof(fn1), "ndis_N15_%s", namebody);
                // snprintf(fn2, sizeof(fn2), "ndis_N14_%s", namebody);
                // snprintf(fn3, sizeof(fn3), "ndis_C13_%s", namebody);
                // snprintf(fn4, sizeof(fn4), "ndis_C12_%s", namebody);
                // snprintf(fn5, sizeof(fn5), "ndis_B11_%s", namebody);
                // snprintf(fn6, sizeof(fn6), "ndis_B10_%s", namebody);
                // snprintf(fn7, sizeof(fn7), "ndis_Be10_%s", namebody);
                // snprintf(fn8, sizeof(fn8), "ndis_Be9_%s", namebody);
                // snprintf(fn9, sizeof(fn9), "ndis_Be7_%s", namebody);
                // write_array_to_bin(fn0, ndis_O16, NR*NZ);
                // write_array_to_bin(fn1, ndis_N15, NR*NZ);
                // write_array_to_bin(fn2, ndis_N14, NR*NZ);           
                // write_array_to_bin(fn3, ndis_C13, NR*NZ);           
                // write_array_to_bin(fn4, ndis_C12, NR*NZ);           
                // write_array_to_bin(fn5, ndis_B11, NR*NZ);           
                // write_array_to_bin(fn6, ndis_B10, NR*NZ);           
                // write_array_to_bin(fn7, ndis_Be10, NR*NZ);           
                // write_array_to_bin(fn8, ndis_Be9, NR*NZ); 
                // write_array_to_bin(fn9, ndis_Be7, NR*NZ);
            }
            // write_array_to_bin(fn, BC_ratio_M, 6);
            // write_array_to_bin(fn, BeC_ratio_M, 6);
            // write_array_to_bin(fn, Be10_Be9_ratio_M, 6);
            // write_array_to_bin(fn, &chi2, 1);
        write_array_to_bin(fn, M_flux, 60);
        }
    }
    free(ndis_O16);
    free(ndis_N15);
    free(ndis_N14);
    free(ndis_C13);
    free(ndis_C12);
    free(ndis_B11);
    free(ndis_B10);
    free(ndis_Be10);
    free(ndis_Be9);
    free(ndis_Be7);
    free(ndis_H);
    return 0;
}



