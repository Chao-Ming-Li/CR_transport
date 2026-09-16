#ifndef FUNC_H
#define FUNC_H

#include <stdio.h>

#define NR 361
#define NZ 171
#define NT 200001
#define dT 50000.0 
extern const int NR1, NZ1;
extern const double dR1, dZ1, dR2, dZ2, R_disk, Z_disk, R_inj;
extern double R[NR], Z[NZ], dR[NR], dZ[NZ];

typedef struct {
    int NR_BH, NZ_BH;
    double dR_BH, dZ_BH, dT_BH, Dd_BH, Dh_BH, ss_time_BH, ndis_H_BH, R_disk_BH;
    char description_BH[128];
} BinHeader;

// Cranck-Nicolson scheme matrix coefficients
typedef struct {

    double alpha_R[NR];
    double alpha_Z[NZ];

    double CN_Ra[NR];
    double CN_Rb[NR];
    double CN_Rc[NR];

    double CN_Za[NZ];
    double CN_Zb[NZ];
    double CN_Zc[NZ];

} CNMatrix;

typedef struct {

    double CN_Rcp[NR];
    double CN_Rdp[NR];

    double CN_Zcp[NZ];
    double CN_Zdp[NZ];

    double *ndis_temp;

} CNWorkspace;


void initialize_grids_linear(void);
void initialize_grids_log(void);
void initialize_to_disk(double *ndis, const double norm, const double R_disk);
void initialize_to_src(double *ndis);
void initialize_to_0(double *ndis);
void initialize_from_file(const char *filename, double *ndis);
void initialize_H(double *ndis_H, const double nH);
void apply_boundary_conditions(double *ndis);
void initialize_velocity(double *vR, double *vZ);
void CN_scheme(double *ndis, double *src, double *ndis_H, double *dT_tau, double *temp, double beta, double sigma, double tau_decay, double D);
void solve_diff_adv_equation_CN(double *ndis_C12, double *ndis_B11, double *ndis_H, double D, double Rg);

void write_array2D_to_txt(const char *filename, const double *ndis);
void write_array1D_to_txt(const char *filename, const double *N);
void write_array_to_bin(const char *filename, const double *arr, const int size);

void CN_scheme(double *ndis, double *src, double *ndis_H, double *dT_tau, double *temp, double beta, double sigma, double tau_decay, double D);
void solve_diff_adv_equation_CN(double *ndis_C12, double *ndis_B11, double *ndis_H, double D, double Rg);

void CN_scheme_opt(double *ndis, double *src, double *dT_tau, CNMatrix *matrix, CNWorkspace *workspace);

void solve_diff_adv_equation_CN_opt(double *ndis_C12, double *ndis_B11, double *ndis_H, double D, double Rg);


#endif
