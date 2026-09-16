#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diff_adv_func.h"
int main(void) {

    initialize_grids_log();
    // Dynamically allocate memory for ndis and temp
    double *ndis_C12 = malloc(NR * NZ * sizeof(double));
    double *ndis_B11 = malloc(NR * NZ * sizeof(double));
    double *ndis_H = malloc(NR * NZ * sizeof(double));
    double *vR = malloc(NR * NZ * sizeof(double));
    double *vZ = malloc(NR * NZ * sizeof(double));

    double gas_type = 0.001, D0, delta, Rg, Dd;
    D0 = 28.0, delta = 1./3., Rg = 10.0;
    Dd = pow(10.0, D0) * pow(Rg, delta) / 3.086e21 / 3.086e21 * 3.15576e7; // in kpc^2 / yr;
    initialize_to_disk(ndis_C12, 1.0, R_disk); 
    initialize_to_0(ndis_B11);
    initialize_H(ndis_H, gas_type);
    initialize_velocity(vR, vZ);

    // solve_diff_adv_equation_CN(ndis_C12, ndis_B11, ndis_H, Dd, Rg);         
    solve_diff_adv_equation_CN_opt(ndis_C12, ndis_B11, ndis_H, Dd, Rg);         

    free(ndis_C12);
    free(ndis_B11);
    free(ndis_H);
    return 0;
}



