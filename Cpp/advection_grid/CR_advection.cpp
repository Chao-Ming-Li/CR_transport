#include "field.hpp"
#include "CR_advection.hpp"  
#include "initialization.hpp"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <span>
#include <cmath>   

// void initialize_grids_log(void){
//     double Rs = 10.0;
//     double Zs = 1.0;
//     for (int i = -2; i < NR + 3; i++) {
//         if (i * 0.1 <= Rs) { Rf(i) = i * 0.1;} 
//         else { Rf(i) = pow(10.0, ((double)i - 100.0) / 200.0 + 1.0 ); }
//     }
//     for (int j = -2; j < NZ + 3; j++) {
//         if (j * 0.02 <= Zs){ Zf(j) = j * 0.02; } 
//         else{ Zf(j) = pow(10.0, ((double)j - 50.0) / 100.0 ); }
//     }
//     for (int i = -2; i < NR + 2; i++){
//         dRf(i) = Rf(i + 1) - Rf(i);
//         Rc(i) = (Rf(i) + Rf(i+1)) / 2.0;
//     }
//     for (int j = -2; j < NZ + 2; j++){
//         dZf(j) = Zf(j + 1) - Zf(j);
//         Zc(j) = (Zf(j) + Zf(j+1)) / 2.0;
//     }
//     for (int i = -2; i < NR + 1; i++){ dRc(i) = Rc(i + 1) - Rc(i); }
//     for (int j = -2; j < NZ + 1; j++){ dZc(j) = Zc(j + 1) - Zc(j); }
// }

// void initialize_grids_log_inv(void){
//     double Rs = 10.0;  // boundary between linear and log grid in R
//     double Zs = 5.0;  // boundary between linear and log grid in Z 
//     Rf(0) = 0.0;
//     for (int i = 0; i < NR + 2; i++) {
//         if (Rf(i) < Rs) { dRf(i) = 0.01 * pow(1.1, i); } 
//         else { dRf(i) = 1.0; }
//         Rf(i + 1) = Rf(i) + dRf(i);
//         Rc(i) = (Rf(i) + Rf(i + 1)) / 2.0;
//     }
//     Rf(-1) = -Rf(1);
//     Rf(-2) = -Rf(2);
//     dRf(-1) = dRf(0);
//     dRf(-2) = dRf(1);
//     Rc(-1) = -Rc(0);
//     Rc(-2) = -Rc(1);
//     Zf(0) = 0.0;
//     for (int i = 0; i < NZ + 2; i++) {
//         if (Zf(i) < Zs) { dZf(i) = 0.01 * pow(1.2, i); } 
//         else { dZf(i) = 1.0; }
//         Zf(i + 1) = Zf(i) + dZf(i);
//         Zc(i) = (Zf(i) + Zf(i + 1)) / 2.0;
//     }
//     Zf(-1) = -Zf(1);
//     Zf(-2) = -Zf(2);
//     dZf(-1) = dZf(0);
//     dZf(-2) = dZf(1);
//     Zc(-1) = -Zc(0);
//     Zc(-2) = -Zc(1);
//     for (int i = -2; i < NR + 1; i++){ dRc(i) = Rc(i + 1) - Rc(i); }
//     for (int j = -2; j < NZ + 1; j++){ dZc(j) = Zc(j + 1) - Zc(j); }
//     for (int i = -2; i < NZ + 3; i++){ std::cout<<Zf(i)<<" "; }
//     for (int i = -2; i < NZ + 2; i++){ std::cout<<dZf(i)<<" "; }
//     for (int i = -2; i < NZ + 2; i++){ std::cout<<Zc(i)<<" "; }
// }


namespace {
void validate_legacy_field(const Field2D& field, const Grid2D& grid,
                         Field2D::Location location = Field2D::Location::Center)
{
    if (!field.matches(grid, location)) {
        throw std::invalid_argument("Field dimensions and location must match the grid");
    }
}
} // namespace

void apply_boundary_conditions(Field2D& ndis, const Grid2D& grid) {
    validate_legacy_field(ndis, grid);
    for (int i = -Field2D::NG; i < ndis.nR() + Field2D::NG; i++) {
        ndis(i, ndis.nz()) = 0.0;     // ghost cell in upper absorbing boundary in z
        ndis(i, ndis.nz() + 1) = 0.0; // ghost cell in upper absorbing boundary in z
        ndis(i, -1) = ndis(i, 0);     // ghost cell in lower reflecting boundary in z
        ndis(i, -2) = ndis(i, 1);     // ghost cell in lower reflecting boundary in z
    }
    for (int j = -Field2D::NG; j < ndis.nz() + Field2D::NG; j++) {
        ndis(ndis.nR(), j) = 0.0;     // ghost cell in outer absorbing boundary in R
        ndis(ndis.nR() + 1, j) = 0.0; // ghost cell in outer absorbing boundary in R
        ndis(-1, j) = ndis(0, j);     // ghost cell in inner reflecting boundary in R
        ndis(-2, j) = ndis(1, j);     // ghost cell in inner reflecting boundary in R
    }
}

void write_array_to_bin(const std::string& filename, Field2D& arr, const std::size_t size) {
    if (size > arr.data.size()) {
        throw std::invalid_argument("Requested output exceeds field storage");
    }
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

void advection_TVD(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT, const Grid2D& grid) {
    validate_legacy_field(ndis, grid);
    validate_legacy_field(temp, grid);
    validate_legacy_field(vR, grid, Field2D::Location::RadialFace);
    validate_legacy_field(vZ, grid, Field2D::Location::VerticalFace);

    // Ensure ghost cells are set before computing fluxes
    apply_boundary_conditions(ndis, grid); 

    //Ensure zero flux at the inner-R/lower-z boundaries
    for (int i = -Field2D::NG; i < ndis.nR() + Field2D::NG; i++) {
        vZ(i, 0) = 0.0; // lower-z boundary
    }
    for (int j = -Field2D::NG; j < ndis.nz() + Field2D::NG; j++) {
        vR(0, j) = 0.0; // inner-R boundary
    }

    // CFL condition check, start from 1 because velocity is defined at faces, vR(0,j) = 0, vZ(i,0) = 0, so the first cell is not used for CFL check  
    for (int i = 1; i < grid.nR(); i++) {
        for (int j = 1; j < grid.nz(); j++) {
            double dt_R = grid.R().center_distance(i) / (fabs(vR(i, j)) + 1e-12);
            double dt_Z = grid.z().center_distance(j) / (fabs(vZ(i, j)) + 1e-12);
            double dt_min = std::min(dt_R, dt_Z);
            if (dT > dt_min) {
                throw std::runtime_error("CFL condition violated: dT exceeds minimum cell crossing time");
            }
        }
    }

    // 2D TVD scheme for advection in cylindrical coordinates
    // ----- R direction -----
    for (int i = 0; i < grid.nR(); i++) {
        for (int j = 0; j < grid.nz(); j++) {
            double n_ip_R, flux_ip_R = 0.0; // i+1/2
            double n_im_R, flux_im_R = 0.0; // i-1/2
            double phi_R; // TVD Limiter
            // ----- i+1/2 -----

            if (vR(i+1, j) >= 0.0) {
                double ratio_ip_R = (ndis(i,j) - ndis(i-1,j)) / grid.R().center_distance(i - 1) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * grid.R().center_distance(i); 
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i,j) + 0.5 * phi_R * (ndis(i+1,j) - ndis(i,j)) * grid.R().width(i) / grid.R().center_distance(i); 
            }
            else {
                double ratio_ip_R = (ndis(i+2,j) - ndis(i+1,j)) / grid.R().center_distance(i + 1) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * grid.R().center_distance(i); 
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i+1,j) - 0.5 * phi_R * ( ndis(i+1,j) - ndis(i,j)) * grid.R().width(i + 1) / grid.R().center_distance(i);
            }
            flux_ip_R = grid.R().face(i + 1) * vR(i+1, j) * n_ip_R;
            
            // ----- i-1/2 -----
            if (vR(i,j) >= 0.0) {
                double ratio_im_R = (ndis(i-1,j) - ndis(i-2,j)) / grid.R().center_distance(i - 2) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * grid.R().center_distance(i - 1); 
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i-1,j) + 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j)) / grid.R().center_distance(i - 1) * grid.R().width(i - 1);
            }
            else {
                double ratio_im_R = (ndis(i+1,j) - ndis(i,j)) / grid.R().center_distance(i) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * grid.R().center_distance(i - 1); 
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i,j) - 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j)) / grid.R().center_distance(i - 1) * grid.R().width(i);
            }
            flux_im_R = grid.R().face(i) * vR(i,j) * n_im_R;

            temp(i,j) = ndis(i,j) + dT * (flux_im_R - flux_ip_R) / grid.R().width(i) / grid.R().center(i);
            if (temp(i,j) < 0.0) {
                throw std::runtime_error("Negative density encountered after R advection step");
            }
        }
    }

    apply_boundary_conditions(temp, grid);

    // ----- z direction -----
    for (int i = 0; i < grid.nR(); i++) {
        for (int j = 0; j < grid.nz(); j++) {
            double n_ip_z, flux_ip_z = 0.0; // i+1/2
            double n_im_z, flux_im_z = 0.0; // i-1/2
            double phi_z; // TVD Limiter
            
            // ----- i+1/2 -----
            if (vZ(i, j+1) >= 0.0) {
                double ratio_ip_z = (temp(i,j) - temp(i,j-1)) / grid.z().center_distance(j - 1)/ (temp(i,j+1) - temp(i,j) + 1e-12) * grid.z().center_distance(j);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j) + 0.5 * phi_z * (temp(i,j+1) -  temp(i,j)) / grid.z().center_distance(j) * grid.z().width(j);
            } else {
                double ratio_ip_z = (temp(i,j+2) - temp(i,j+1)) / grid.z().center_distance(j + 1) / (temp(i, j+1) - temp(i,j) + 1e-12) * grid.z().center_distance(j);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j+1) - 0.5 * phi_z * (temp(i,j+1) - temp(i,j)) / grid.z().center_distance(j) * grid.z().width(j + 1);
            }
            flux_ip_z = vZ(i, j+1) * n_ip_z;

            // ----- i-1/2 -----
            if (vZ(i, j) >= 0.0) {
                double ratio_im_z = (temp(i,j-1) - temp(i,j-2)) / grid.z().center_distance(j - 2) / (temp(i,j) - temp(i,j-1) + 1e-12) * grid.z().center_distance(j - 1);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j-1) + 0.5 * phi_z * (temp(i,j) - temp(i,j-1)) / grid.z().center_distance(j - 1) * grid.z().width(j - 1);
            } else {
                double ratio_im_z = (temp(i,j+1) - temp(i,j)) / grid.z().center_distance(j) / (temp(i,j) - temp(i,j-1) + 1e-12) * grid.z().center_distance(j - 1);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j) - 0.5 * phi_z * (temp(i,j) - temp(i,j-1)) / grid.z().center_distance(j - 1) * grid.z().width(j);
            }
            flux_im_z = vZ(i, j) * n_im_z;
            ndis(i,j) = temp(i,j) + dT * (flux_im_z - flux_ip_z) / grid.z().width(j);
            if (ndis(i,j) < 0.0) {
                throw std::runtime_error("Negative density encountered after Z advection step");
            }
        }
    }
    apply_boundary_conditions(ndis, grid);
}


void advection_TVD_R(Field2D& ndis, Field2D& temp, Field2D& vR, double dT, const Grid2D& grid) {

    // 2D TVD scheme for advection in cylindrical coordinates
    // ----- R direction -----
    for (int i = 0; i < grid.nR(); i++) {
        for (int j = 0; j < grid.nz(); j++) {
            double n_ip_R, flux_ip_R = 0.0; // i+1/2
            double n_im_R, flux_im_R = 0.0; // i-1/2
            double phi_R; // TVD Limiter
            // ----- i+1/2 -----

            if (vR(i+1, j) >= 0.0) {
                double ratio_ip_R = (ndis(i,j) - ndis(i-1,j)) / grid.R().center_distance(i - 1) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * grid.R().center_distance(i); 
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i,j) + 0.5 * phi_R * (ndis(i+1,j) - ndis(i,j)) * grid.R().width(i) / grid.R().center_distance(i); 
            }
            else {
                double ratio_ip_R = (ndis(i+2,j) - ndis(i+1,j)) / grid.R().center_distance(i + 1) / (ndis(i+1,j) - ndis(i,j) + 1e-12) * grid.R().center_distance(i); 
                Van_leer_limiter(phi_R, ratio_ip_R);
                n_ip_R = ndis(i+1,j) - 0.5 * phi_R * ( ndis(i+1,j) - ndis(i,j)) * grid.R().width(i + 1) / grid.R().center_distance(i);
            }
            flux_ip_R = grid.R().face(i + 1) * vR(i+1, j) * n_ip_R;
            
            // ----- i-1/2 -----
            if (vR(i,j) >= 0.0) {
                double ratio_im_R = (ndis(i-1,j) - ndis(i-2,j)) / grid.R().center_distance(i - 2) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * grid.R().center_distance(i - 1); 
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i-1,j) + 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j)) / grid.R().center_distance(i - 1) * grid.R().width(i - 1);
            }
            else {
                double ratio_im_R = (ndis(i+1,j) - ndis(i,j)) / grid.R().center_distance(i) / (ndis(i,j) - ndis(i-1,j) + 1e-12) * grid.R().center_distance(i - 1); 
                Van_leer_limiter(phi_R, ratio_im_R);
                n_im_R = ndis(i,j) - 0.5 * phi_R * (ndis(i,j) - ndis(i-1,j)) / grid.R().center_distance(i - 1) * grid.R().width(i);
            }
            flux_im_R = grid.R().face(i) * vR(i,j) * n_im_R;

            temp(i,j) = ndis(i,j) + dT * (flux_im_R - flux_ip_R) / grid.R().width(i) / grid.R().center(i);
            if (temp(i,j) < 0.0) {
                throw std::runtime_error("Negative density encountered after R advection step");
            }
        }
    }
    for (int i = 0; i < grid.nR(); i++) {
        for (int j = 0; j < grid.nz(); j++) {
            ndis(i,j) = temp(i,j);
        }
    }
    apply_boundary_conditions(ndis, grid);
}


void advection_TVD_Z(Field2D& ndis, Field2D& temp, Field2D& vZ, double dT, const Grid2D& grid) {

    // 2D TVD scheme for advection in cylindrical coordinates

    // ----- z direction -----
    for (int i = 0; i < grid.nR(); i++) {
        for (int j = 0; j < grid.nz(); j++) {
            double n_ip_z, flux_ip_z = 0.0; // i+1/2
            double n_im_z, flux_im_z = 0.0; // i-1/2
            double phi_z; // TVD Limiter
            
            // ----- i+1/2 -----
            if (vZ(i, j+1) >= 0.0) {
                double ratio_ip_z = (ndis(i,j) - ndis(i,j-1)) / grid.z().center_distance(j - 1)/ (ndis(i,j+1) - ndis(i,j) + 1e-12) * grid.z().center_distance(j);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = ndis(i,j) + 0.5 * phi_z * (ndis(i,j+1) -  ndis(i,j)) / grid.z().center_distance(j) * grid.z().width(j);
            } else {
                double ratio_ip_z = (ndis(i,j+2) - ndis(i,j+1)) / grid.z().center_distance(j + 1) / (ndis(i, j+1) - ndis(i,j) + 1e-12) * grid.z().center_distance(j);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = ndis(i,j+1) - 0.5 * phi_z * (ndis(i,j+1) - ndis(i,j)) / grid.z().center_distance(j) * grid.z().width(j + 1);
            }
            flux_ip_z = vZ(i, j+1) * n_ip_z;

            // ----- i-1/2 -----
            if (vZ(i, j) >= 0.0) {
                double ratio_im_z = (ndis(i,j-1) - ndis(i,j-2)) / grid.z().center_distance(j - 2) / (ndis(i,j) - ndis(i,j-1) + 1e-12) * grid.z().center_distance(j - 1);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = ndis(i,j-1) + 0.5 * phi_z * (ndis(i,j) - ndis(i,j-1)) / grid.z().center_distance(j - 1) * grid.z().width(j - 1);
            } else {
                double ratio_im_z = (ndis(i,j+1) - ndis(i,j)) / grid.z().center_distance(j) / (ndis(i,j) - ndis(i,j-1) + 1e-12) * grid.z().center_distance(j - 1);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = ndis(i,j) - 0.5 * phi_z * (ndis(i,j) - ndis(i,j-1)) / grid.z().center_distance(j - 1) * grid.z().width(j);
            }
            flux_im_z = vZ(i, j) * n_im_z;
            temp(i,j) = ndis(i,j) + dT * (flux_im_z - flux_ip_z) / grid.z().width(j);
            if (temp(i,j) < 0.0) {
                throw std::runtime_error("Negative density encountered after Z advection step");
            }
        }
    }
    for (int i = 0; i < grid.nR(); i++) {
        for (int j = 0; j < grid.nz(); j++) {
            ndis(i,j) = temp(i,j);
        }
    }
    apply_boundary_conditions(ndis, grid);
}


// inline double minmod(double a,double b)
// {
//     if(a*b<=0.0)
//         return 0.0;
//     return fabs(a)<fabs(b)?a:b;
// }

// inline double MC_limiter(double sL,double sC,double sR)
// {
//     return minmod(2.0 * minmod(sL, sR), sC);
// }

// void advection_PLM(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT, const Grid2D& grid) {
//     validate_legacy_field(ndis, grid);
//     validate_legacy_field(temp, grid);
//     validate_legacy_field(vR, grid, Field2D::Location::RadialFace);
//     validate_legacy_field(vZ, grid, Field2D::Location::VerticalFace);
//     // 2D PLM scheme for advection in cylindrical coordinates
//     Field2D flux_R(ndis.nR(), ndis.nz()), flux_Z(ndis.nR(), ndis.nz()); // slopes for PLM
//     // ----- R direction -----
//     for (int i = 0; i < grid.nR() + 1; i++) {
//         for (int j = 0; j < grid.nz() + 1; j++) {
//             if (vR(i, j) >= 0.0){
//                 double sF_L = (ndis(i - 1, j) - ndis(i - 2, j)) / grid.R().center_distance(i - 2);
//                 double sF_R = (ndis(i, j) - ndis(i - 1, j)) / grid.R().center_distance(i - 1);
//                 double sC = (ndis(i, j) - ndis(i - 2, j)) / (grid.R().center_distance(i - 1) + grid.R().center_distance(i - 2));
//                 double slope = MC_limiter(sF_L, sC, sF_R);
//                 double ndis_face = ndis(i - 1, j) + slope * (grid.R().face(i) - grid.R().center(i - 1)); // reconstruct left state at face i
//                 flux_R(i, j) = grid.R().face(i) * vR(i, j) * ndis_face;
//             }
//             else{
//                 double sF_L = (ndis(i, j) - ndis(i - 1, j)) / grid.R().center_distance(i - 1);
//                 double sF_R = (ndis(i + 1, j) - ndis(i, j)) / grid.R().center_distance(i);
//                 double sC = (ndis(i + 1, j) - ndis(i - 1, j)) / (grid.R().center_distance(i - 1) + grid.R().center_distance(i));
//                 double slope = MC_limiter(sF_L, sC, sF_R);
//                 double ndis_face = ndis(i, j) - slope * (grid.R().center(i) - grid.R().face(i)); // reconstruct right state at face i
//                 flux_R(i, j) = grid.R().face(i) * vR(i, j) * ndis_face;
//             }
//         }
//     }
//     for (int i = 0; i < grid.nR(); i++) {
//         for (int j = 0; j < grid.nz(); j++) {
//             temp(i, j) = ndis(i, j) - dT * (flux_R(i + 1, j) - flux_R(i, j)) / grid.R().width(i) / grid.R().center(i);
//         }
//     }
//     apply_boundary_conditions(temp, grid);
//     // ----- Z direction -----
//     for (int i = 0; i < grid.nR() + 1; i++) {
//         for (int j = 0; j < grid.nz() + 1; j++) {
//             if (vZ(i, j) >= 0.0){
//                 double sF_L = (temp(i, j - 1) - temp(i, j - 2)) / grid.z().center_distance(j - 2);
//                 double sF_R = (temp(i, j) - temp(i, j - 1)) / grid.z().center_distance(j - 1);
//                 double sC = (temp(i, j) - temp(i, j - 2)) / (grid.z().center_distance(j - 1) + grid.z().center_distance(j - 2));
//                 double slope = MC_limiter(sF_L, sC, sF_R);
//                 double temp_face = temp(i, j - 1) + slope * (grid.z().face(j) - grid.z().center(j - 1)); // reconstruct left state at face j
//                 flux_Z(i, j) = vZ(i, j) * temp_face;
//             }
//             else{
//                 double sF_L = (temp(i, j) - temp(i, j - 1)) / grid.z().center_distance(j - 1);
//                 double sF_R = (temp(i, j + 1) - temp(i, j)) / grid.z().center_distance(j);
//                 double sC = (temp(i, j + 1) - temp(i, j - 1)) / (grid.z().center_distance(j - 1) + grid.z().center_distance(j));
//                 double slope = MC_limiter(sF_L, sC, sF_R);
//                 double temp_face = temp(i, j) - slope * (grid.z().center(j) - grid.z().face(j)); // reconstruct right state at face j
//                 flux_Z(i, j) = vZ(i, j) * temp_face;
//             }
//         }
//     }
//     for (int i = 0; i < grid.nR(); i++) {
//         for (int j = 0; j < grid.nz(); j++) {
//             ndis(i, j) = temp(i, j) - dT * (flux_Z(i, j + 1) - flux_Z(i, j)) / grid.z().width(j);
//         }
//     }
//     apply_boundary_conditions(ndis, grid);
// }


void solve_advection_equation(Field2D& ndis_C, Field2D& ndis_H, Field2D& vR, Field2D& vZ, const Grid2D& grid) {
    validate_legacy_field(ndis_C, grid);
    validate_legacy_field(ndis_H, grid);
    validate_legacy_field(vR, grid, Field2D::Location::RadialFace);
    validate_legacy_field(vZ, grid, Field2D::Location::VerticalFace);

    Field2D src_C(grid);
    Field2D temp(grid);
    initialize_CR_source(src_C, grid);
    apply_boundary_conditions(ndis_C, grid);

    // for (int t = 0; t < NT; t++){
    //     advection_TVD(ndis_C, temp, vR, vZ, DT, grid);
    //     write_array_to_bin("ndis_C.bin", ndis_C, ndis_C.data.size());
    // }

    // Strang splitting with alternate directional updates
    if (NT > 0) {
        advection_TVD_R(ndis_C, temp, vR, 0.5 * DT, grid);

        for (int step = 0; step < NT; ++step) {
            advection_TVD_Z(ndis_C, temp, vZ, DT, grid);

            const double dt_R = (step == NT - 1) ? 0.5 * DT : DT;
            advection_TVD_R(ndis_C, temp, vR, dt_R, grid);
        }
    }
    
    write_array_to_bin("ndis_C_t_ssad_1.bin", ndis_C, ndis_C.data.size());

}

