#include "CR_diffusion.hpp"  
#include "CNWorkspace.hpp"  

#include <omp.h>
#include <fstream>
#include <iostream>
#include <span>

static double c = 3e10, yr2s = 3.1536e7, ratio = 40.46285, tau_Be10 = 1.387e6, tau_stable = 1e100, rest_energy_p = 0.9315, enhancement_factor = 1.3; // tau in seconds, tau_stable for stable isotopes, rest_energy_p in GeV, enhancement_factor for cross-section to account for heavier ISM nuclei

static double sigma_C = 273e-27, sigma_B = 265e-27, sigma_C2B = 56e-27;

Array1D Rc(NR), Zc(NZ), Rf(NR + 1), Zf(NZ + 1), dR(NR), dZ(NZ); // Rf, Zf are face values of cells, while Rc, Zc are center values of cells

void initialize_grids_log(void){
    double Rs = 10.0;
    double Zs = 1.0;
    
    for (int i = -2; i < NR + 3; i++) {
        if (i * 0.1 <= Rs) { Rf(i) = i * 0.1;} 
        else { Rf(i) = pow(10.0, ((double)i - 100.0) / 200.0 + 1.0 ); }
    }
    for (int j = -2; j < NZ + 3; j++) {
        if (j * 0.02 <= Zs){ Zf(j) = j * 0.02; } 
        else{ Zf(j) = pow(10.0, ((double)j - 50.0) / 50.0 ); }
    }

    for (int i = -2; i < NR + 2; i++){
        dR(i) = Rf(i + 1) - Rf(i);
        Rc(i) = (Rf(i) + Rf(i+1)) / 2.0;
    }
    for (int j = -2; j < NZ + 2; j++){
        dZ(j) = Zf(j + 1) - Zf(j);
        Zc(j) = (Zf(j) + Zf(j+1)) / 2.0;
    }
}

void initialize_disk_source(Field2D& ndis) {
    // boundary values are initialized to 0 automatically
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis(i, j) = (Rc(i) <= 8.0 && Zc(j) <= 0.1) ? 1.0 * DT : 0.0;
        }
    }
}

void initialize_CN_workspace(CNWorkspace& CN, double D,double dT){

    for(int j = 0; j < NR - 1; j++){ CN.alpha_R(j) = D * dT / ( dR(j) + dR(j+1)); }

    for(int k = 0; k < NZ - 1; k++){ CN.alpha_Z(k) = D * dT / ( dZ(k) + dZ(k+1)); }

    for (int j = 0; j < NR - 1; j++){     
        CN.Ra(j) = - CN.alpha_R(j) / dR(j) * (1.0 - 0.5 * dR(j + 1) / Rc(j));
        CN.Rb(j) = 1.0 + CN.alpha_R(j) * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) * (1.0 - (dR(j + 1) - dR(j)) / 2.0 / Rc(j));
        CN.Rc(j) = - CN.alpha_R(j) / dR(j + 1) * (1.0 + 0.5 * dR(j) / Rc(j));
    }
    
    CN.Ra(-1) = 0.0; CN.Rb(-1) = 1.0; CN.Rc(-1) = -1.0; CN.Rd(-1) = 0.0;

    for (int k = 0; k < NZ - 1; k++){
        CN.Za(k) = - CN.alpha_Z(k) / dZ(k);
        CN.Zb(k) = 1.0 + CN.alpha_Z(k) * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k));
        CN.Zc(k) = - CN.alpha_Z(k) / dZ(k + 1);
    }

    CN.Za(-1) = 0.0; CN.Zb(-1) = 1.0; CN.Zc(-1) = -1.0; CN.Zd(-1) = 0.0;
}

void apply_boundary_conditions(Field2D& ndis) {
    for (int i = -2; i < NR + 2; i++) {
        ndis(i, NZ) = 0.0;
        ndis(i, NZ + 1) = 0.0;
        ndis(i, -1) = ndis(i, 0);
        ndis(i, -2) = ndis(i, 1);
    }
    for (int j = -2; j < NZ + 2; j++) {
        ndis(NR, j) = 0.0;
        ndis(NR + 1, j) = 0.0;
        ndis(-1, j) = ndis(0, j);
        ndis(-2, j) = ndis(1, j);
    }
}

void initialize_H(Field2D& ndis_H, const double nH) {
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis_H(i, j) = nH;
        }
    }
}

void write_array_to_bin(const std::string& filename, Field2D& arr, const std::size_t size) {
    // std::ios::binary 表示以二进制模式打开
    // std::ios::app    表示追加模式（对应原本的 "ab"）
    std::ofstream file(filename, std::ios::binary | std::ios::app);
    
    if (!file) {
        // 使用 std::cerr 代替 perror
        std::cerr << "File opening failed: " << filename << std::endl;
        return;
    }
    
    // reinterpret_cast 将 double 指针转换为 char 指针，符合 write 接口
    file.write(reinterpret_cast<const char*>(arr.data.data()), size * sizeof(double));
    
    // file 在离开作用域时会自动关闭，不需要显式调用 file.close()
}

void solve_diffusion_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, double D, double Rg, double dT) {

    CNWorkspace CN;
    initialize_CN_workspace(CN, D, dT); 

    Field2D src_C;
    Field2D src_B;
    Field2D dT_tau_C;
    Field2D dT_tau_B;

    initialize_disk_source(src_C);

    double beta = Rg / sqrt(Rg * Rg + rest_energy_p * rest_energy_p); // beta = v/c

    // 2. fragmentation and decay term
    for (int j = 0; j < NR; j++) {
        for (int k = 0; k < NZ; k++) {
            dT_tau_C(j, k) = dT * yr2s * (ndis_H(j, k) * sigma_C * beta * c + 1.0 / tau_stable);
            dT_tau_B(j, k) = dT * yr2s * (ndis_H(j, k) * sigma_B * beta * c + 1.0 / tau_stable);
        }
    }

    for (int t = 0; t < NT; t++) {

        CN.CN_scheme_2nd_order(ndis_C, src_C, ndis_H, dT_tau_C);

        for (int j = 0; j < NR; j++) {
            for (int k = 0; k < NZ; k++) {
                src_B(j, k) = ndis_C(j, k) * sigma_C2B * beta * ndis_H(j, k) * c * yr2s * DT;
            }
        }
        CN.CN_scheme_2nd_order(ndis_B, src_B, ndis_H, dT_tau_B);
    }
}
