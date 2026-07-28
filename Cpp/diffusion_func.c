#include "CR_transport.hpp"  
#include <omp.h>
#include <fstream>
#include <iostream>
#include <span>

static double c = 3e10, yr2s = 3.1536e7, ratio = 40.46285, tau_Be10 = 1.387e6, tau_stable = 1e100, rest_energy_p = 0.9315, enhancement_factor = 1.3; // tau in seconds, tau_stable for stable isotopes, rest_energy_p in GeV, enhancement_factor for cross-section to account for heavier ISM nuclei

static double sigma_C = 273e-27, sigma_B = 265e-27, sigma_C2B = 56e-27;

double R[NR], Z[NZ], dR[NR], dZ[NZ];

void initialize_grids_log(void){
    double Rc = 10.0;
    double Zc = 1.0;
    
    for (int i = 0; i < NR; i++) {
        if (i * 0.1 <= Rc)
        {
            R[i] = i * 0.1;
            dR[i] = 0.1;
        } 
        else {
            R[i] = pow(10.0, ((double)i - 100.0) / 200.0 + 1.0 ); // 400 steps for every order of magnitude, i.e. R[i] = 10**1.0000, R[i+1] = 10**1.0025, R[i+2] = 10**1.005, R3 = 10**1.0075....... 
            dR[i] = R[i] - R[i - 1]; 
        }
        // printf("R is %.5f and dR is %.5f : ", R[i], dR[i]);
    }
    for (int k = 0; k < NZ; k++) {
        if (k * 0.02 <= Zc){
            Z[k] = k * 0.02;
            dZ[k] = 0.02;
        } 
        else{
            Z[k] = pow(10.0, ((double)k - 50.0) / 50.0 - 1.0 ); //  100 steps for every order of magnitude: i.e. R[i] = 10**1.00, R[i+1] = 10**1.01, R[i+2] = 10**1.02, R[i+3] = 10**1.03....... 
            dZ[k] = Z[k] - Z[k - 1];
        }
        // printf("Z is %.5f and dZ is %.5f : ", Z[k], dZ[k]);
    }
}

void initialize_disk_source(Field2D& ndis) {
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis(i, j) = (i <= 80 && j <= 10) ? 1.0 : 0.0;
        }
    }
}

void initialize_to_0(Field2D& ndis) {
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis(i, j) = 0.0;
        }
    }
}

void apply_boundary_conditions(Field2D& ndis) {
    for (int i = 0; i < NR; i++) {
        ndis(i, NZ - 1) = 0.0;
        ndis(i, 0) = ndis(i, 1);
    }
    for (int j = 0; j < NZ; j++) {
        ndis(NR - 1, j) = 0.0;
        ndis(0, j) = ndis(1, j);
    }
}

void initialize_H(Field2D& ndis_H, const double nH) {
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis_H(i, j) = nH;
        }
    }
}


inline void thomas_solve( 
    std::span<const double> a, 
    std::span<const double> b, 
    std::span<const double> c, 
    std::span<const double> d, 
    std::span<double> x, 
    std::span<double> cp, 
    std::span<double> dp)
{
    std::size_t n = b.size();
    if (n == 0) return;

    // Modify the first coefficients
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // Forward sweep
    for (std::size_t i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * cp[i - 1];
        
        if (denom == 0.0) {
            throw std::runtime_error("Zero pivot in Thomas algorithm at index " + std::to_string(i));
        }
        
        double inv_denom = 1.0 / denom;
        cp[i] = (i < n - 1) ? c[i] * inv_denom : 0.0;
        dp[i] = (d[i] - a[i] * dp[i - 1]) * inv_denom;
    }
    // Backward substitution
    x[n - 1] = dp[n - 1];
    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
}

void CN_scheme_1st_order(Field2D& ndis, 
               const Field2D& src, 
               const Field2D& ndis_H, 
               Field2D& dT_tau, 
               Field2D& temp, 
               double beta, double sigma, double tau_decay, double D) 
{
    // 1. define ADI coefficients
    std::vector<double> alpha_R(NR), alpha_Z(NZ);
    for (int i = 0; i < NR - 1; i++) alpha_R[i] = D * dT / (dR[i] + dR[i + 1]);
    for (int i = 0; i < NZ - 1; i++) alpha_Z[i] = D * dT / (dZ[i] + dZ[i + 1]);

    // 2. fragmentation and decay term
    for (int j = 0; j < NR; j++) {
        for (int k = 0; k < NZ; k++) {
            dT_tau(j, k) = dT * yr2s * (ndis_H(j, k) * sigma * beta * c + 1.0 / tau_decay);
        }
    }

    // 3. 准备托马斯算法的工作缓冲区（分配在循环外部，避免重复开辟）
    std::vector<double> CN_a(std::max(NR, NZ)), CN_b(std::max(NR, NZ)), CN_c(std::max(NR, NZ));
    std::vector<double> CN_d(std::max(NR, NZ)), CN_x(std::max(NR, NZ)), CN_cp(std::max(NR, NZ)), CN_dp(std::max(NR, NZ));

    // 为托马斯算法建立安全的 std::span 视图，后面直接把视图切片丢进算法中
    std::span<double> s_a(CN_a), s_b(CN_b), s_c(CN_c), s_d(CN_d), s_x(CN_x), s_cp(CN_cp), s_dp(CN_dp);

    // ==========================================
    // Step 1: Implicit R, Explicit Z
    // ==========================================
    for (int k = 1; k < NZ - 1; k++) {     
        for (int j = 1; j < NR - 1; j++) {
            CN_d[j] = alpha_Z[k] / dZ[k] * ndis(j, k - 1) 
                    + (1.0 - alpha_Z[k] * (dZ[k + 1] + dZ[k]) / (dZ[k + 1] * dZ[k])) * ndis(j, k) 
                    + alpha_Z[k] / dZ[k + 1] * ndis(j, k + 1) 
                    - ndis(j, k) * dT_tau(j, k) / 4.0;
        }
        
        for (int j = 1; j < NR - 1; j++){     
            CN_a[j] = - alpha_R[j] * (1.0 / dR[j] - 0.5 / R[j]);
            CN_b[j] = 1.0 + alpha_R[j] * (dR[j + 1] + dR[j]) / (dR[j + 1] * dR[j]) + dT_tau(j, k) / 4.0;
            CN_c[j] = - alpha_R[j] * (1.0 / dR[j + 1] + 0.5 / R[j]);
        }
        
        CN_a[0] = 0.0; CN_b[0] = 1.0; CN_c[0] = -1.0; CN_d[0] = 0.0;

        // 使用 subspan 截取当前维度大小 (NR - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(0, NR - 1), s_b.subspan(0, NR - 1), s_c.subspan(0, NR - 1), s_d.subspan(0, NR - 1), s_x.subspan(0, NR - 1), s_cp.subspan(0, NR - 1), s_dp.subspan(0, NR - 1));
        
        CN_x[NR - 1] = 0.0; // enforce absorbing boundary at Rmax
        
        for (int j = 0; j < NR; j++) temp(j, k) = CN_x[j];
    }
    apply_boundary_conditions(temp);

    // ==========================================
    // Step 2: Implicit Z, Explicit R
    // ==========================================
    for (int j = 1; j < NR - 1; j++){
        for (int k = 1; k < NZ - 1; k++){
            CN_d[k] = alpha_R[j] * (1.0 / dR[j] - 0.5 / R[j]) * temp(j - 1, k) + (1.0 - alpha_R[j] * (dR[j + 1] + dR[j]) / (dR[j + 1] * dR[j])) * temp(j, k) + alpha_R[j] * (1.0 / dR[j + 1] + 0.5 / R[j]) * temp(j + 1, k) - temp(j, k) * dT_tau(j, k) / 4.0;
        }

        for (int k = 1; k < NZ - 1; k++){
            CN_a[k] = - alpha_Z[k] / dZ[k];
            CN_b[k] = 1.0 + alpha_Z[k] * (dZ[k + 1] + dZ[k]) / (dZ[k + 1] * dZ[k]) + dT_tau(j, k) / 4.0;
            CN_c[k] = - alpha_Z[k] / dZ[k + 1];
        }
        
        CN_a[0] = 0.0; CN_b[0] = 1.0; CN_c[0] = -1.0; CN_d[0] = 0.0;

        // 使用 subspan 截取当前维度大小 (NZ - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(0, NZ - 1), s_b.subspan(0, NZ - 1), s_c.subspan(0, NZ - 1), s_d.subspan(0, NZ - 1), s_x.subspan(0, NZ - 1), s_cp.subspan(0, NZ - 1), s_dp.subspan(0, NZ - 1));
        
        CN_x[NZ - 1] = 0.0; // enforce absorbing boundary at Zmax
        
        for (int k = 0; k < NZ; k++) ndis(j, k) = CN_x[k] + src(j, k);
    }
    apply_boundary_conditions(ndis);
}

void CN_scheme_2nd_order(Field2D& ndis, 
               const Field2D& src, 
               const Field2D& ndis_H, 
               Field2D& dT_tau, 
               Field2D& temp, 
               double beta, double sigma, double tau_decay, double D) 
{
    // 1. define ADI coefficients
    std::vector<double> alpha_R(NR), alpha_Z(NZ);
    for (int i = 0; i < NR - 1; i++) alpha_R[i] = D * dT / (dR[i] + dR[i + 1]);
    for (int i = 0; i < NZ - 1; i++) alpha_Z[i] = D * dT / (dZ[i] + dZ[i + 1]);

    // 2. fragmentation and decay term
    for (int j = 0; j < NR; j++) {
        for (int k = 0; k < NZ; k++) {
            dT_tau(j, k) = dT * yr2s * (ndis_H(j, k) * sigma * beta * c + 1.0 / tau_decay);
        }
    }

    // 3. 准备托马斯算法的工作缓冲区（分配在循环外部，避免重复开辟）
    std::vector<double> CN_a(std::max(NR, NZ)), CN_b(std::max(NR, NZ)), CN_c(std::max(NR, NZ));
    std::vector<double> CN_d(std::max(NR, NZ)), CN_x(std::max(NR, NZ)), CN_cp(std::max(NR, NZ)), CN_dp(std::max(NR, NZ));

    // 为托马斯算法建立安全的 std::span 视图，后面直接把视图切片丢进算法中
    std::span<double> s_a(CN_a), s_b(CN_b), s_c(CN_c), s_d(CN_d), s_x(CN_x), s_cp(CN_cp), s_dp(CN_dp);

    // ==========================================
    // Step 1: Implicit R, Explicit Z
    // ==========================================
    for (int k = 1; k < NZ - 1; k++) {     
        for (int j = 1; j < NR - 1; j++) {
            CN_d[j] = alpha_Z[k] / dZ[k] * ndis(j, k - 1) 
                    + (1.0 - alpha_Z[k] * (dZ[k + 1] + dZ[k]) / (dZ[k + 1] * dZ[k])) * ndis(j, k) 
                    + alpha_Z[k] / dZ[k + 1] * ndis(j, k + 1) 
                    - ndis(j, k) * dT_tau(j, k) / 4.0;
        }
        
        for (int j = 1; j < NR - 1; j++){     
            CN_a[j] = - alpha_R[j] / dR[j] * (1.0 - 0.5 * dR[j + 1] / R[j]);
            CN_b[j] = 1.0 + alpha_R[j] * (dR[j + 1] + dR[j]) / (dR[j + 1] * dR[j]) * (1.0 - (dR[j+1] - dR[j]) / 2.0 / R[j]) + dT_tau(j, k) / 4.0;
            CN_c[j] = - alpha_R[j] / dR[j + 1] * (1.0 + 0.5 * dR[j] / R[j]);
        }
        
        CN_a[0] = 0.0; CN_b[0] = 1.0; CN_c[0] = -1.0; CN_d[0] = 0.0;

        // 使用 subspan 截取当前维度大小 (NR - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(0, NR - 1), s_b.subspan(0, NR - 1), s_c.subspan(0, NR - 1), s_d.subspan(0, NR - 1), s_x.subspan(0, NR - 1), s_cp.subspan(0, NR - 1), s_dp.subspan(0, NR - 1));
        
        CN_x[NR - 1] = 0.0; // enforce absorbing boundary at Rmax
        
        for (int j = 0; j < NR; j++) temp(j, k) = CN_x[j];
    }
    apply_boundary_conditions(temp);

    // ==========================================
    // Step 2: Implicit Z, Explicit R
    // ==========================================
    for (int j = 1; j < NR - 1; j++){
        for (int k = 1; k < NZ - 1; k++){
            CN_d[k] = alpha_R[j] / dR[j] * (1.0 - 0.5 * dR[j + 1]/ R[j]) * temp(j - 1, k) + (1.0 - alpha_R[j] * (dR[j + 1] + dR[j]) / (dR[j + 1] * dR[j]) * (1.0 - (dR[j+1] - dR[j]) / 2.0 / R[j])) * temp(j, k) + alpha_R[j] / dR[j + 1] * (1.0 + 0.5 * dR[j] / R[j]) * temp(j + 1, k) - temp(j, k) * dT_tau(j, k) / 4.0;
        }

        for (int k = 1; k < NZ - 1; k++){
            CN_a[k] = - alpha_Z[k] / dZ[k];
            CN_b[k] = 1.0 + alpha_Z[k] * (dZ[k + 1] + dZ[k]) / (dZ[k + 1] * dZ[k]) + dT_tau(j, k) / 4.0;
            CN_c[k] = - alpha_Z[k] / dZ[k + 1];
        }
        
        CN_a[0] = 0.0; CN_b[0] = 1.0; CN_c[0] = -1.0; CN_d[0] = 0.0;

        // 使用 subspan 截取当前维度大小 (NZ - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(0, NZ - 1), s_b.subspan(0, NZ - 1), s_c.subspan(0, NZ - 1), s_d.subspan(0, NZ - 1), s_x.subspan(0, NZ - 1), s_cp.subspan(0, NZ - 1), s_dp.subspan(0, NZ - 1));
        
        CN_x[NZ - 1] = 0.0; // enforce absorbing boundary at Zmax
        
        for (int k = 0; k < NZ; k++) ndis(j, k) = CN_x[k] + src(j, k);
    }
    apply_boundary_conditions(ndis);
}



void solve_diffusion_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, double D, double Rg) {
    Field2D src_C(NR, NZ);
    Field2D src_B(NR, NZ);
    Field2D dT_tau(NR, NZ);
    Field2D temp(NR, NZ);
    initialize_disk_source(src_C);
    initialize_to_0(ndis_B);
    initialize_to_0(dT_tau);
    initialize_to_0(temp);
    double beta = Rg / sqrt(Rg * Rg + rest_energy_p * rest_energy_p); // beta = v/c

    for (int t = 0; t < NT; t++) {
        CN_scheme_2nd_order(ndis_C, src_C, ndis_H, dT_tau, temp, beta, sigma_C, tau_stable, D);
        for (int j = 0; j < NR; j++) {
            for (int k = 0; k < NZ; k++) {
                src_B(j, k) = ndis_C(j, k) * sigma_C2B * beta * ndis_H(j, k) * c * yr2s * dT;
            }
        }
        CN_scheme_2nd_order(ndis_B, src_B, ndis_H, dT_tau, temp, beta, sigma_B, tau_stable, D);
    }
}

