#include "field.hpp"
#include "initialization.hpp"
#include "CR_transport_solver.hpp"
#include "CR_advection.hpp"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <span>
#include <cmath>   
#include <algorithm>
#include <limits>

void apply_boundary_conditions(Field2D& ndis, const Grid2D& grid) {
    if (!ndis.matches(grid, Field2D::Location::Centroid))
        throw std::invalid_argument("Boundary density must match the grid");
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

void solve_advection_equation(Field2D& ndis_C, Field2D& ndis_H, const Field2D& vR, const Field2D& vZ, const Grid2D& grid, double D) {

    // Gas density is reserved for future interaction/loss terms.
    (void)ndis_H;
    cr_advection::Options options;
    options.reconstruction = cr_advection::Reconstruction::PLM;
    options.limiter = cr_advection::Limiter::MC;
    options.integrator = cr_advection::TimeIntegrator::SSPRK2;
    options.splitting = cr_advection::Splitting::Unsplit;

    cr_advection::Solver advection(grid, options);
    cr_diffusion::Options diffusion_options;
    cr_diffusion::Solver diffusion(grid, D, DT, diffusion_options);

    apply_boundary_conditions(ndis_C, grid);
    // Sequential advection/diffusion is first-order operator splitting.
    for (int t = 0; t < NT; ++t) {
        advection.advance(ndis_C, vR, vZ, DT);
        diffusion.advance(ndis_C);
        if (t % 10 == 0)
            write_array_to_bin(D == 0.0 ? "ndis_C_transport_dR50_dz10_r1.03_dt1e4_MC_fix.bin"
                                      : "ndis_C_advection_diffusion.bin", ndis_C, ndis_C.data.size());
    }
}
