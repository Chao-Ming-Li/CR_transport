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

void apply_boundary_conditions(Field2D &ndis, const Grid2D &grid) {
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

void write_array_to_bin(const std::string &filename, Field2D &arr, const std::size_t size) {
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
    file.write(reinterpret_cast<const char *>(arr.data.data()), size * sizeof(double));

    // file 在离开作用域时会自动关闭，不需要显式调用 file.close()
}

void solve_advection_equation(Field2D &ndis_C, Field2D &ndis_H, const Field2D &vR,
                              const Field2D &vZ, const Grid2D &grid, double D) {

    double max_advection_rate = 0.0;
    double max_diffusion_rate = 0.0;
    for (int i = 0; i < grid.nR(); ++i) {
        const double radial_volume = grid.R().center(i) * grid.R().width(i);
        const double diffusion_R =
            D *
            ((i > 0 ? grid.R().face(i) / grid.radial_centroid_distance(i - 1) : 0.0) +
             grid.R().face(i + 1) / grid.radial_centroid_distance(i)) /
            radial_volume;
        for (int j = 0; j < grid.nz(); ++j) {
            // Lower boundaries reflect: their velocities do not contribute flux.
            const double vr_lower = i > 0 ? vR(i, j) : 0.0;
            const double vr_upper = vR(i + 1, j);
            const double vz_lower = j > 0 ? vZ(i, j) : 0.0;
            const double vz_upper = vZ(i, j + 1);

            const double advection_rate =
                (grid.R().face(i) * std::max(-vr_lower, 0.0) +
                 grid.R().face(i + 1) * std::max(vr_upper, 0.0)) /
                    radial_volume +
                (std::max(-vz_lower, 0.0) + std::max(vz_upper, 0.0)) / grid.z().width(j);
            max_advection_rate = std::max(max_advection_rate, advection_rate);
            const double diffusion_z = D *
                                       ((j > 0 ? 1.0 / grid.z().center_distance(j - 1) : 0.0) +
                                        1.0 / grid.z().center_distance(j)) /
                                       grid.z().width(j);
            max_diffusion_rate = std::max(max_diffusion_rate, std::max(diffusion_R, diffusion_z));
        }
    }

    // MC-limited PLM face values can reach twice the donor cell average.
    // This outgoing-flux bound also accounts for cylindrical/stretched cells.
    constexpr double advection_cfl_limit = 0.5;
    if (DT * max_advection_rate > advection_cfl_limit)
        throw std::invalid_argument("Advection CFL condition violated: reduce DT to at most " +
                                    std::to_string(advection_cfl_limit / max_advection_rate) +
                                    " years");

    // Each ADI explicit half-step is nonnegative when (DT/2)*rate <= 1.
    // Exceeding this bound affects positivity, not implicit diffusion stability.
    if (0.5 * DT * max_diffusion_rate > 1.0)
        std::cerr << "Warning: diffusion positivity condition exceeded; DT should be at most "
                  << 2.0 / max_diffusion_rate << " years to guarantee nonnegative ADI steps.\n";

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
                                        : "ndis_C_advection_diffusion.bin",
                               ndis_C, ndis_C.data.size());
    }
}
