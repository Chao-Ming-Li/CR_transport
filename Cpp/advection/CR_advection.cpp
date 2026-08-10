#include "CR_advection.hpp"  
#include <omp.h>
#include <fstream>
#include <iostream>
#include <span>
#include <cmath>   
static double c = 3e10, yr2s = 3.1536e7, ratio = 40.46285, tau_Be10 = 1.387e6, tau_stable = 1e100, rest_energy_p = 0.9315, enhancement_factor = 1.3; // tau in seconds, tau_stable for stable isotopes, rest_energy_p in GeV, enhancement_factor for cross-section to account for heavier ISM nuclei

static double sigma_C = 273e-27, sigma_B = 265e-27, sigma_C2B = 56e-27;

Array1D Rc(NR), Zc(NZ), Rf(NR + 1), Zf(NZ + 1), dRf(NR), dZf(NZ), dRc(NR - 1), dZc(NZ - 1); // Rf, Zf are face values of cells, while Rc, Zc are center values of cells



// void initialize_grids_linear(void){
    
//     for (int i = -2; i < NR + 3; i++) {
//         Rf(i) = i * 0.1;
//         Zf(i) = i * 0.1;
//     }

//     for (int i = -2; i < NR + 2; i++){
//         dR(i) = Rf(i + 1) - Rf(i);
//         Rc(i) = (Rf(i) + Rf(i+1)) / 2.0;
//     }
//     for (int j = -2; j < NZ + 2; j++){
//         dZ(j) = Zf(j + 1) - Zf(j);
//         Zc(j) = (Zf(j) + Zf(j+1)) / 2.0;
//     }
// }

void initialize_grids_log(void){
    double Rs = 10.0;
    double Zs = 1.0;
    
    for (int i = -2; i < NR + 3; i++) {
        if (i * 0.1 <= Rs) { Rf(i) = i * 0.1;} 
        else { Rf(i) = pow(10.0, ((double)i - 100.0) / 200.0 + 1.0 ); }
    }
    for (int j = -2; j < NZ + 3; j++) {
        if (j * 0.02 <= Zs){ Zf(j) = j * 0.02; } 
        else{ Zf(j) = pow(10.0, ((double)j - 50.0) / 100.0 ); }
    }

    for (int i = -2; i < NR + 2; i++){
        dRf(i) = Rf(i + 1) - Rf(i);
        Rc(i) = (Rf(i) + Rf(i+1)) / 2.0;
    }
    for (int j = -2; j < NZ + 2; j++){
        dZf(j) = Zf(j + 1) - Zf(j);
        Zc(j) = (Zf(j) + Zf(j+1)) / 2.0;
    }

    for (int i = -2; i < NR + 1; i++){ dRc(i) = Rc(i + 1) - Rc(i); }
    for (int j = -2; j < NZ + 1; j++){ dZc(j) = Zc(j + 1) - Zc(j); }

}

void initialize_grids_log_inv(void){
    double Rs = 10.0;  // boundary between linear and log grid in R
    double Zs = 5.0;  // boundary between linear and log grid in Z
    
    Rf(0) = 0.0;
    for (int i = 0; i < NR + 2; i++) {
        if (Rf(i) < Rs) { dRf(i) = 0.01 * pow(1.1, i); } 
        else { dRf(i) = 1.0; }
        Rf(i + 1) = Rf(i) + dRf(i);
        Rc(i) = (Rf(i) + Rf(i + 1)) / 2.0;
    }

    Rf(-1) = -Rf(1);
    Rf(-2) = -Rf(2);
    dRf(-1) = dRf(0);
    dRf(-2) = dRf(1);
    Rc(-1) = -Rc(0);
    Rc(-2) = -Rc(1);


    Zf(0) = 0.0;
    for (int i = 0; i < NZ + 2; i++) {
        if (Zf(i) < Zs) { dZf(i) = 0.01 * pow(1.2, i); } 
        else { dZf(i) = 1.0; }
        Zf(i + 1) = Zf(i) + dZf(i);
        Zc(i) = (Zf(i) + Zf(i + 1)) / 2.0;
    }

    Zf(-1) = -Zf(1);
    Zf(-2) = -Zf(2);
    dZf(-1) = dZf(0);
    dZf(-2) = dZf(1);
    Zc(-1) = -Zc(0);
    Zc(-2) = -Zc(1);

    for (int i = -2; i < NR + 1; i++){ dRc(i) = Rc(i + 1) - Rc(i); }
    for (int j = -2; j < NZ + 1; j++){ dZc(j) = Zc(j + 1) - Zc(j); }

    for (int i = -2; i < NZ + 3; i++){ std::cout<<Zf(i)<<" "; }
    for (int i = -2; i < NZ + 2; i++){ std::cout<<dZf(i)<<" "; }
    for (int i = -2; i < NZ + 2; i++){ std::cout<<Zc(i)<<" "; }


}


void initialize_disk_source(Field2D& ndis) {
    // boundary values are initialized to 0 automatically
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            // ndis(i, j) = (Rc(i) <= 10.0 && Zc(j) <= 1.0) ? 1.0 * DT : 0.0;
            // double r2 = Rc(i) * Rc(i) + Zc(j) * Zc(j);
            // ndis(i,j) = 1.0 * DT / r2;
            double R1 = fabs(Rc(i) - 30.0);
            double Z1 = fabs(Zc(j) - 30.0);
            ndis(i, j) = 1000.0 * exp(- R1 * R1 / 9.0 - Z1 * Z1 / 9.0) * DT; // exponential disk source, scale length 5 kpc, scale height 1 kpc

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
            // vR(i, j) = 3e-7;  // kpc / yr
            vZ(i ,j) = 3e-7;  // kpc / yr
            vR(i, j) = 0.0;
            // vZ(i, j) = 0.0;
        }
    }
    for (int i = 0; i <= NR; i++){ vZ(i, 0) = 0.0; }
    // for (int j = 0; j <= NZ; j++){ vR(0, j) = 0.0; }
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

void advection_TVD(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT) {
    // 2D TVD scheme for advection in cylindrical coordinates

    // ----- R direction -----
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            double n_ip_R, flux_ip_R = 0.0; // i+1/2
            double n_im_R, flux_im_R = 0.0; // i-1/2
            double phi_R; // TVD Limiter
            // ----- i+1/2 -----

            if (vR(i+1, j) >= 0.0) {
                double ratio_ip_R = (ndis(i,j) - ndis(i-1,j)) / dRc(i - 1) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * dRc(i); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i,j) + 0.5 * phi_R * (ndis(i+1,j) - ndis(i,j)) * dRf(i) / dRc(i); //flux in cylindrical coordinates is R * vR * ndis
            }
            else {
                double ratio_ip_R = (ndis(i+2,j) - ndis(i+1,j)) / dRc(i + 1) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * dRc(i); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i+1,j) - 0.5 * phi_R * ( ndis(i+1,j) - ndis(i,j)) * dRf(i + 1) / dRc(i);
            }
            flux_ip_R = Rf(i + 1) * vR(i+1, j) * n_ip_R;
            
            // ----- i-1/2 -----
            if (vR(i,j) >= 0.0) {
                double ratio_im_R = (ndis(i-1,j) - ndis(i-2,j)) / dRc(i - 2) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * dRc(i - 1); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i-1,j) + 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j)) / dRc(i - 1) * dRf(i - 1);
            }
            else {
                double ratio_im_R = (ndis(i+1,j) - ndis(i,j)) / dRc(i) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * dRc(i - 1); // for non-uniform grids, r is upstream slope / downstream slope, instead of the ratio of differences
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i,j) - 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j)) / dRc(i - 1) * dRf(i);
            }
            flux_im_R = Rf(i) * vR(i,j) * n_im_R;

            temp(i,j) = ndis(i,j) - dT * (flux_ip_R - flux_im_R) / dRf(i) / Rc(i);
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
                double ratio_ip_z = (temp(i,j) - temp(i,j-1)) / dZc(j - 1)/ (temp(i,j+1) - temp(i,j) + 1e-12) * dZc(j);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j) + 0.5 * phi_z * (temp(i,j+1) -  temp(i,j)) / dZc(j) * dZf(j);
            } else {
                double ratio_ip_z = (temp(i,j+2) - temp(i,j+1)) / dZc(j + 1) / (temp(i, j+1) - temp(i,j) + 1e-12) * dZc(j);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j+1) - 0.5 * phi_z * (temp(i,j+1) - temp(i,j)) / dZc(j) * dZf(j + 1);
            }
            flux_ip_z = vZ(i, j+1) * n_ip_z;

            // ----- i-1/2 -----
            if (vZ(i, j) >= 0.0) {
                double ratio_im_z = (temp(i,j-1) - temp(i,j-2)) / dZc(j - 2) / (temp(i,j) - temp(i,j-1) + 1e-12) * dZc(j - 1);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j-1) + 0.5 * phi_z * (temp(i,j) - temp(i,j-1)) / dZc(j - 1) * dZf(j - 1);
            } else {
                double ratio_im_z = (temp(i,j+1) - temp(i,j)) / dZc(j) / (temp(i,j) - temp(i,j-1) + 1e-12) * dZc(j - 1);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j) - 0.5 * phi_z * (temp(i,j) - temp(i,j-1)) / dZc(j - 1) * dZf(j);
            }
            flux_im_z = vZ(i, j) * n_im_z;
            ndis(i,j) = temp(i,j) - dT * (flux_ip_z - flux_im_z) / dZf(j);
        }
    }
    apply_boundary_conditions(ndis);
}

inline double minmod(double a,double b)
{
    if(a*b<=0.0)
        return 0.0;
    return fabs(a)<fabs(b)?a:b;
}

inline double MC_limiter(double sL,double sC,double sR)
{
    return minmod(2.0 * minmod(sL, sR), sC);
}


void advection_PLM(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT) {

    // 2D PLM scheme for advection in cylindrical coordinates
    Field2D flux_R, flux_Z; // slopes for PLM

    // ----- R direction -----
    for (int i = 0; i < NR + 1; i++) {
        for (int j = 0; j < NZ + 1; j++) {
            if (vR(i, j) >= 0.0){
                double sF_L = (ndis(i - 1, j) - ndis(i - 2, j)) / dRc(i - 2);
                double sF_R = (ndis(i, j) - ndis(i - 1, j)) / dRc(i - 1);
                double sC = (ndis(i, j) - ndis(i - 2, j)) / (dRc(i - 1) + dRc(i - 2));
                double slope = MC_limiter(sF_L, sC, sF_R);
                double ndis_face = ndis(i - 1, j) + slope * (Rf(i) - Rc(i - 1)); // reconstruct left state at face i
                flux_R(i, j) = Rf(i) * vR(i, j) * ndis_face;
            }
            else{
                double sF_L = (ndis(i, j) - ndis(i - 1, j)) / dRc(i - 1);
                double sF_R = (ndis(i + 1, j) - ndis(i, j)) / dRc(i);
                double sC = (ndis(i + 1, j) - ndis(i - 1, j)) / (dRc(i - 1) + dRc(i));
                double slope = MC_limiter(sF_L, sC, sF_R);
                double ndis_face = ndis(i, j) - slope * (Rc(i) - Rf(i)); // reconstruct right state at face i
                flux_R(i, j) = Rf(i) * vR(i, j) * ndis_face;
            }
        }
    }
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            temp(i, j) = ndis(i, j) - dT * (flux_R(i + 1, j) - flux_R(i, j)) / dRf(i) / Rc(i);
        }
    }
    apply_boundary_conditions(temp);

    // ----- Z direction -----
    for (int i = 0; i < NR + 1; i++) {
        for (int j = 0; j < NZ + 1; j++) {
            if (vZ(i, j) >= 0.0){
                double sF_L = (temp(i, j - 1) - temp(i, j - 2)) / dZc(j - 2);
                double sF_R = (temp(i, j) - temp(i, j - 1)) / dZc(j - 1);
                double sC = (temp(i, j) - temp(i, j - 2)) / (dZc(j - 1) + dZc(j - 2));
                double slope = MC_limiter(sF_L, sC, sF_R);
                double temp_face = temp(i, j - 1) + slope * (Zf(j) - Zc(j - 1)); // reconstruct left state at face j
                flux_Z(i, j) = vZ(i, j) * temp_face;
            }
            else{
                double sF_L = (temp(i, j) - temp(i, j - 1)) / dZc(j - 1);
                double sF_R = (temp(i, j + 1) - temp(i, j)) / dZc(j);
                double sC = (temp(i, j + 1) - temp(i, j - 1)) / (dZc(j - 1) + dZc(j));
                double slope = MC_limiter(sF_L, sC, sF_R);
                double temp_face = temp(i, j) - slope * (Zc(j) - Zf(j)); // reconstruct right state at face j
                flux_Z(i, j) = vZ(i, j) * temp_face;
            }

        }
    }
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            ndis(i, j) = temp(i, j) - dT * (flux_Z(i, j + 1) - flux_Z(i, j)) / dZf(j);
        }
    }
    apply_boundary_conditions(ndis);

}


void solve_advection_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, Field2D& vR, Field2D& vZ){
    Field2D src_C;
    Field2D temp;
    initialize_disk_source(src_C);
    for (int t = 0; t < NT; t++){
        advection_PLM(ndis_C, temp, vR, vZ, DT);
        // for (int j = 0; j < NR; j++){
        //     for (int k = 0; k < NZ; k++){
        //         // ndis_B(j, k) += ndis_C(j, k) * sigma_C2B * ndis_H(j, k) * c * yr2s * DT;
        //         ndis_C(j, k) += src_C(j, k);
        //     }
        // }
        // advection_PLM(ndis_B, temp, vR, vZ, DT);
        if (t % 100 == 0){
            write_array_to_bin("ndis_C_log_gauss_inj_dt10000_vZ_dZ2_PLM.bin", ndis_C, (NR + 4)*(NZ + 4));
        }
    }
}

