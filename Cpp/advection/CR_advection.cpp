#include "CR_advection.hpp"  
#include <omp.h>
#include <fstream>
#include <iostream>
#include <span>
#include <cmath>   
static double c = 3e10, yr2s = 3.1536e7, ratio = 40.46285, tau_Be10 = 1.387e6, tau_stable = 1e100, rest_energy_p = 0.9315, enhancement_factor = 1.3; // tau in seconds, tau_stable for stable isotopes, rest_energy_p in GeV, enhancement_factor for cross-section to account for heavier ISM nuclei

static double sigma_C = 273e-27, sigma_B = 265e-27, sigma_C2B = 56e-27;

Array1D Rc(NR), Zc(NZ), Rf(NR + 1), Zf(NZ + 1), dR(NR), dZ(NZ); // Rf, Zf are face values of cells, while Rc, Zc are center values of cells


void initialize_grids_linear(void){
    
    for (int i = -2; i < NR + 3; i++) {
        Rf(i) = i * 0.1;
        Zf(i) = i * 0.1;
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
            ndis(i, j) = (Rc(i) <= 10.0 && Zc(j) <= 10.0) ? 1.0 * DT : 0.0;
            // double r2 = Rc(i) * Rc(i) + Zc(j) * Zc(j);
            // ndis(i,j) = 1.0 * DT / r2;
        }
    }
}

void initialize_wind_velocity(Field2D& vR, Field2D& vZ) {
    // TVD only need velocity at cell faces, we use 0 - NR to store NR+1 face values
    for (int i = 0; i <= NR; i++) {
        for (int j = 0; j <= NZ; j++) {
            double r = pow(Rf(i) * Rf(i) + Zf(j) * Zf(j), 0.5) / 1.0 + 1e-5; // normalize to 1 kpc, avoid 0
            // double v = 1000.0 * exp(- 1.0 / r) / r / r; // peak at r = 0.5 kpc, around 500 km/s, decreases as 1/r^2 at large r, and goes to 0 at r = 0
            // vR(i, j) = v * (Rf(i) / r); // velocity field in R direction, changes from -1 to 1 across the grid
            // vZ(i, j) = v * (Zf(j) / r); // velocity field in Z direction, changes from -1 to 1 across the grid
            vR(i, j) = 3e-7;  // kpc / yr
            // vZ(i ,j) = 3e-7;  // kpc / yr
            // vR(i, j) = 0.0;
            vZ(i, j) = 0.0;
        }
    }
    // for (int i = 0; i <= NR; i++){ vZ(i, 0) = 0.0; }
    for (int j = 0; j <= NZ; j++){ vR(0, j) = 0.0; }
}

void apply_boundary_conditions(Field2D& ndis) {
    for (int i = -2; i < NR + 2; i++) {
        ndis(i, NZ) = 0.0;
        ndis(i, NZ + 1) = 0.0;
        ndis(i, -1) = ndis(i, 0);
        ndis(i, -2) = ndis(i, 1);
        // ndis(i, -1) = 0.0;
        // ndis(i, -2) = 0.0;
    }
    for (int j = -2; j < NZ + 2; j++) {
        ndis(NR, j) = 0.0;
        ndis(NR + 1, j) = 0.0;
        ndis(-1, j) = ndis(0, j);
        ndis(-2, j) = ndis(1, j);
        // ndis(-1, j) = 0.0;
        // ndis(-2, j) = 0.0; 
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
    // std::ios::app    追加模式; std::ios::trunc    覆写模式
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

inline void Van_leer_limiter(double& phi, double ratio) {
    if (ratio <= 0.0) {
        phi = 0.0;
    } else {
        phi = (ratio + fabs(ratio)) / (1.0 + fabs(ratio));
    }
}

void advection_TVD_solution(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT) {
    // 2D TVD scheme for advection in cylindrical coordinates

    // ----- R direction -----
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            double n_ip_R, flux_ip_R = 0.0; // i+1/2
            double n_im_R, flux_im_R = 0.0; // i-1/2
            double phi_R; // TVD Limiter
            // ----- i+1/2 -----

            if (vR(i+1, j) >= 0.0) {
                double ratio_ip_R = (ndis(i,j) - ndis(i-1,j)) / dR(i) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * dR(i + 1); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i,j) + 0.5 * phi_R * (ndis(i+1,j) - ndis(i,j)); //flux in cylindrical coordinates is R * vR * ndis
            }
            else {
                double ratio_ip_R = (ndis(i+2,j) - ndis(i+1,j)) / dR(i + 2) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * dR(i + 1); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i+1,j) - 0.5 * phi_R * ( ndis(i+1,j) - ndis(i,j));
            }
            flux_ip_R = Rf(i + 1) * vR(i+1, j) * n_ip_R;
            
            // ----- i-1/2 -----
            if (vR(i,j) >= 0.0) {
                double ratio_im_R = (ndis(i-1,j) - ndis(i-2,j)) / dR(i - 1) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * dR(i); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i-1,j) + 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j));
            }
            else {
                double ratio_im_R = (ndis(i+1,j) - ndis(i,j)) / dR(i + 1) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * dR(i); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i,j) - 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j));
            }
            flux_im_R = Rf(i) * vR(i,j) * n_im_R;

            temp(i,j) = ndis(i,j) - dT * (flux_ip_R - flux_im_R) / dR(i) / Rc(i);
        }
    }

    apply_boundary_conditions(temp);

    // ----- z direction -----
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            double n_ip_z, flux_ip_z = 0.0; // i+1/2
            double n_im_z, flux_im_z = 0.0; // i-1/2
            double phi_z; // TVD Limiter
            
            // ----- i+1/2 -----
            if (vZ(i, j+1) >= 0.0) {
                double ratio_ip_z = (temp(i,j) - temp(i,j-1)) / dZ(j)/ (temp(i,j+1) - temp(i,j) + 1e-12) * dZ(j + 1);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j) + 0.5 * phi_z * (temp(i,j+1) -  temp(i,j));
            } else {
                double ratio_ip_z = (temp(i,j+2) - temp(i,j+1)) / dZ(j + 2) / (temp(i,j+1) - temp(i,j) + 1e-12) * dZ(j + 1);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j+1) - 0.5 * phi_z * (temp(i,j+1) - temp(i,j));
            }
            flux_ip_z = vZ(i, j+1) * n_ip_z;

            // ----- i-1/2 -----
            if (vZ(i, j) >= 0.0) {
                double ratio_im_z = (temp(i,j-1) - temp(i,j-2)) / dZ(j - 1) / (temp(i,j) - temp(i,j-1) + 1e-12) * dZ(j);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j-1) + 0.5 * phi_z * (temp(i,j) - temp(i,j-1));
            } else {
                double ratio_im_z = (temp(i,j+1) - temp(i,j)) / dZ(j + 1) / (temp(i,j) - temp(i,j-1) + 1e-12) * dZ(j);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j) - 0.5 * phi_z * (temp(i,j) - temp(i,j-1));
            }
            flux_im_z = vZ(i, j) * n_im_z;
            ndis(i,j) = temp(i,j) - dT * (flux_ip_z - flux_im_z) / dZ(j);
        }
    }
    apply_boundary_conditions(ndis);
}

void advection_TVD_R(Field2D& ndis, Field2D& temp, Field2D& vR, double dT) {
    // ----- R direction -----
    for (int i = 2; i < NR - 2; i++) {
        for (int j = 2; j < NZ - 2; j++) {
            double n_ip_R, flux_ip_R = 0.0; // i+1/2
            double n_im_R, flux_im_R = 0.0; // i-1/2
            double phi_R; // TVD Limiter
            // ----- i+1/2 -----
            double vR_ip = (vR(i,j) + vR(i+1,j)) / 2.0; // vR(i,j) + (vR(i+1,j) - vR(i,j)) * (R_ip - R[i]) / (R[i+1] - R[i]) = (vR(i,j) + vR(i+1,j)) / 2.0;

            if (vR_ip >= 0.0) {
                double ratio_ip_R = (ndis(i,j) - ndis(i-1,j)) / dR(i) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * dR(i + 1); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i,j) + 0.5 * phi_R * (ndis(i+1,j) - ndis(i,j)); //flux in cylindrical coordinates is R * vR * ndis
            }
            else {
                double ratio_ip_R = (ndis(i+2,j) - ndis(i+1,j)) / dR(i + 2) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * dR(i + 1); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i+1,j) - 0.5 * phi_R * ( ndis(i+1,j) - ndis(i,j));
            }
            flux_ip_R = Rf(i + 1) * vR_ip * n_ip_R;
            
            // ----- i-1/2 -----
            double vR_im = (vR(i-1,j) + vR(i,j)) / 2.0;

            if (vR_im >= 0.0) {
                double ratio_im_R = (ndis(i-1,j) - ndis(i-2,j)) / dR(i - 1) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * dR(i); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i-1,j) + 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j));
            }
            else {
                double ratio_im_R = (ndis(i+1,j) - ndis(i,j)) / dR(i + 1) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * dR(i); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i,j) - 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j));
            }
            flux_im_R = Rf(i) * vR_im * n_im_R;

            temp(i,j) = ndis(i,j) - dT * (flux_ip_R - flux_im_R) / dR(i) / Rc(i);
        }
    }
    apply_boundary_conditions(temp);
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis(i,j) = temp(i,j);
        }
    }
}

void advection_TVD_Z(Field2D& ndis, Field2D& temp, Field2D& vZ, double dT) {
    // ----- z direction -----
    for (int i = 2; i < NR - 2; i++) {
        for (int j = 2; j < NZ - 2; j++) {
            double n_ip_z, flux_ip_z = 0.0; // i+1/2
            double n_im_z, flux_im_z = 0.0; // i-1/2
            double phi_z; // TVD Limiter
            
            // ----- i+1/2 -----
            double vZ_ip = (vZ(i,j) + vZ(i,j+1)) / 2.0;

            if (vZ_ip >= 0.0) {
                double ratio_ip_z = (ndis(i,j) - ndis(i,j-1)) / dZ(j)/ (ndis(i,j+1) - ndis(i,j) + 1e-12) * dZ(j + 1);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = ndis(i,j) + 0.5 * phi_z * (ndis(i,j+1) -  ndis(i,j));
            } else {
                double ratio_ip_z = (ndis(i,j+2) - ndis(i,j+1)) / dZ(j + 2) / (ndis(i,j+1) - ndis(i,j) + 1e-12) * dZ(j + 1);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = ndis(i,j+1) - 0.5 * phi_z * (ndis(i,j+1) - ndis(i,j));
            }
            flux_ip_z = vZ_ip * n_ip_z;

            // ----- i-1/2 -----
            double vZ_im = (vZ(i,j-1) + vZ(i,j)) / 2.0;
            if (vZ_im >= 0.0) {
                double ratio_im_z = (ndis(i,j-1) - ndis(i,j-2)) / dZ(j - 1) / (ndis(i,j) - ndis(i,j-1) + 1e-12) * dZ(j);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = ndis(i,j-1) + 0.5 * phi_z * (ndis(i,j) - ndis(i,j-1));
            } else {
                double ratio_im_z = (ndis(i,j+1) - ndis(i,j)) / dZ(j + 1) / (ndis(i,j) - ndis(i,j-1) + 1e-12) * dZ(j);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = ndis(i,j) - 0.5 * phi_z * (ndis(i,j) - ndis(i,j-1));
            }
            flux_im_z = vZ_im * n_im_z;
            temp(i,j) = ndis(i,j) - dT * (flux_ip_z - flux_im_z) / dZ(j);
        }
    }
    apply_boundary_conditions(temp);
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis(i,j) = temp(i,j);
        }
    }
}

void solve_advection_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, Field2D& vR, Field2D& vZ){
    Field2D src_C;
    Field2D temp;
    initialize_disk_source(src_C);
    for (int t = 0; t < NT; t++){
        advection_TVD_solution(ndis_C, temp, vR, vZ, DT);
        // for (int j = 0; j < NR; j++){
        //     for (int k = 0; k < NZ; k++){
        //         // ndis_B(j, k) += ndis_C(j, k) * sigma_C2B * ndis_H(j, k) * c * yr2s * DT;
        //         ndis_C(j, k) += src_C(j, k);
        //     }
        // }
        // advection_TVD_solution(ndis_B, temp, vR, vZ, DT);
        if (t % 1000 == 0){
            write_array_to_bin("ndis_C_log_delta_inj_dt100_vR.bin", ndis_C, (NR + 4)*(NZ + 4));
        }
    }
}
