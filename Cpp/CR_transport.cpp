#include "CR_transport.hpp"  
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

void initialize_wind_velocity(Field2D& vR, Field2D& vZ) {
    // TVD only need velocity at cell faces, we use 0 - NR to store NR+1 face values
    for (int i = 0; i <= NR + 1; i++) {
        for (int j = 0; j <= NZ + 1; j++) {
            double r = pow(Rf(i) * Rf(i) + Zf(j) * Zf(j), 0.5) / 1.0 + 1e-5; // normalize to 1 kpc, avoid 0
            double v = 1000.0 * exp(- 1.0 / r) / r / r; // peak at r = 0.5 kpc, around 500 km/s, decreases as 1/r^2 at large r, and goes to 0 at r = 0
            vR(i, j) = v * (Rf(i) / r); // velocity field in R direction, changes from -1 to 1 across the grid
            vZ(i, j) = v * (Zf(j) / r); // velocity field in Z direction, changes from -1 to 1 across the grid
        }
    }
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
        
        if (std::abs(denom) < 1e-15) {
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
               double beta, double sigma, double tau_decay, double D, double dT) 
{
    // 1. define ADI coefficients
    std::vector<double> alpha_R(NR), alpha_Z(NZ);
    for (int i = 0; i < NR - 1; i++) alpha_R[i] = D * dT / (dR(i) + dR(i + 1));
    for (int i = 0; i < NZ - 1; i++) alpha_Z[i] = D * dT / (dZ(i) + dZ(i + 1));

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
            CN_d[j] = alpha_Z[k] / dZ(k) * ndis(j, k - 1) 
                    + (1.0 - alpha_Z[k] * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k))) * ndis(j, k) 
                    + alpha_Z[k] / dZ(k + 1) * ndis(j, k + 1) 
                    - ndis(j, k) * dT_tau(j, k) / 4.0;
        }
        
        for (int j = 1; j < NR - 1; j++){     
            CN_a[j] = - alpha_R[j] * (1.0 / dR(j) - 0.5 / Rc(j));
            CN_b[j] = 1.0 + alpha_R[j] * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) + dT_tau(j, k) / 4.0;
            CN_c[j] = - alpha_R[j] * (1.0 / dR(j + 1) + 0.5 / Rc(j));
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
            CN_d[k] = alpha_R[j] * (1.0 / dR(j) - 0.5 / Rc(j)) * temp(j - 1, k) + (1.0 - alpha_R[j] * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j))) * temp(j, k) + alpha_R[j] * (1.0 / dR(j + 1) + 0.5 / Rc(j)) * temp(j + 1, k) - temp(j, k) * dT_tau(j, k) / 4.0;
        }

        for (int k = 1; k < NZ - 1; k++){
            CN_a[k] = - alpha_Z[k] / dZ(k);
            CN_b[k] = 1.0 + alpha_Z[k] * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k)) + dT_tau(j, k) / 4.0;
            CN_c[k] = - alpha_Z[k] / dZ(k + 1);
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
               double beta, double sigma, double tau_decay, double D, double dT) 
{
    // 1. define ADI coefficients
    std::vector<double> alpha_R(NR), alpha_Z(NZ);
    for (int i = 0; i < NR - 1; i++) alpha_R[i] = D * dT / (dR(i) + dR(i + 1));
    for (int i = 0; i < NZ - 1; i++) alpha_Z[i] = D * dT / (dZ(i) + dZ(i + 1));

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
            CN_d[j] = alpha_Z[k] / dZ(k) * ndis(j, k - 1) 
                    + (1.0 - alpha_Z[k] * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k))) * ndis(j, k) 
                    + alpha_Z[k] / dZ(k + 1) * ndis(j, k + 1) 
                    - ndis(j, k) * dT_tau(j, k) / 4.0;
        }
        
        for (int j = 1; j < NR - 1; j++){     
            CN_a[j] = - alpha_R[j] / dR(j) * (1.0 - 0.5 * dR(j + 1) / Rc(j));
            CN_b[j] = 1.0 + alpha_R[j] * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) * (1.0 - (dR(j + 1) - dR(j)) / 2.0 / Rc(j)) + dT_tau(j, k) / 4.0;
            CN_c[j] = - alpha_R[j] / dR(j + 1) * (1.0 + 0.5 * dR(j) / Rc(j));
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
            CN_d[k] = alpha_R[j] / dR(j) * (1.0 - 0.5 * dR(j + 1)/ Rc(j)) * temp(j - 1, k) + (1.0 - alpha_R[j] * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) * (1.0 - (dR(j + 1) - dR(j)) / 2.0 / Rc(j))) * temp(j, k) + alpha_R[j] / dR(j + 1) * (1.0 + 0.5 * dR(j) / Rc(j)) * temp(j + 1, k) - temp(j, k) * dT_tau(j, k) / 4.0;
        }

        for (int k = 1; k < NZ - 1; k++){
            CN_a[k] = - alpha_Z[k] / dZ(k);
            CN_b[k] = 1.0 + alpha_Z[k] * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k)) + dT_tau(j, k) / 4.0;
            CN_c[k] = - alpha_Z[k] / dZ(k + 1);
        }
        
        CN_a[0] = 0.0; CN_b[0] = 1.0; CN_c[0] = -1.0; CN_d[0] = 0.0;

        // 使用 subspan 截取当前维度大小 (NZ - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(0, NZ - 1), s_b.subspan(0, NZ - 1), s_c.subspan(0, NZ - 1), s_d.subspan(0, NZ - 1), s_x.subspan(0, NZ - 1), s_cp.subspan(0, NZ - 1), s_dp.subspan(0, NZ - 1));
        
        CN_x[NZ - 1] = 0.0; // enforce absorbing boundary at Zmax
        
        for (int k = 0; k < NZ; k++) ndis(j, k) = CN_x[k] + src(j, k);
    }
    apply_boundary_conditions(ndis);
}


void CN_scheme_2nd_order_gc(Field2D& ndis, 
               const Field2D& src, 
               const Field2D& ndis_H, 
               Field2D& dT_tau, 
               Field2D& temp, 
               Array1D alpha_R,
               Array1D alpha_Z,
               double beta, double sigma, double tau_decay, double D, double dT) 
{

    // 3. 准备托马斯算法的工作缓冲区（分配在循环外部，避免重复开辟）
    // std::vector<double> CN_a(std::max(NR, NZ)), CN_b(std::max(NR, NZ)), CN_c(std::max(NR, NZ));
    // std::vector<double> CN_d(std::max(NR, NZ)), CN_x(std::max(NR, NZ)), CN_cp(std::max(NR, NZ)), CN_dp(std::max(NR, NZ));
    Array1D CN_a(std::max(NR, NZ)), CN_b(std::max(NR, NZ)), CN_c(std::max(NR, NZ));
    Array1D CN_d(std::max(NR, NZ)), CN_x(std::max(NR, NZ)), CN_cp(std::max(NR, NZ)), CN_dp(std::max(NR, NZ));

    // 为托马斯算法建立安全的 std::span 视图，后面直接把视图切片丢进算法中
    std::span<double> s_a(CN_a.arr), s_b(CN_b.arr), s_c(CN_c.arr), s_d(CN_d.arr), s_x(CN_x.arr), s_cp(CN_cp.arr), s_dp(CN_dp.arr);

    // ==========================================
    // Step 1: Implicit R, Explicit Z
    // ==========================================
    for (int k = 0; k < NZ - 1; k++) {     
        for (int j = 0; j < NR - 1; j++) {
            CN_d(j) = alpha_Z(k) / dZ(k) * ndis(j, k - 1) 
                    + (1.0 - alpha_Z(k) * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k))) * ndis(j, k) 
                    + alpha_Z(k) / dZ(k + 1) * ndis(j, k + 1) 
                    - ndis(j, k) * dT_tau(j, k) / 4.0;
        }
        
        for (int j = 0; j < NR - 1; j++){     
            CN_a(j) = - alpha_R(j) / dR(j) * (1.0 - 0.5 * dR(j + 1) / Rc(j));
            CN_b(j) = 1.0 + alpha_R(j) * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) * (1.0 - (dR(j + 1) - dR(j)) / 2.0 / Rc(j)) + dT_tau(j, k) / 4.0;
            CN_c(j) = - alpha_R(j) / dR(j + 1) * (1.0 + 0.5 * dR(j) / Rc(j));
        }
        
        CN_a(-1) = 0.0; CN_b(-1) = 1.0; CN_c(-1) = -1.0; CN_d(-1) = 0.0;

        // 使用 subspan 截取当前维度大小 (NR - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(1, NR), s_b.subspan(1, NR), s_c.subspan(1, NR), s_d.subspan(1, NR), s_x.subspan(1, NR), s_cp.subspan(1, NR), s_dp.subspan(1, NR));
        
        CN_x(NR - 1) = 0.0; // enforce absorbing boundary at Rmax
        
        for (int j = 0; j < NR; j++) temp(j, k) = CN_x(j);
    }
    apply_boundary_conditions(temp);

    // ==========================================
    // Step 2: Implicit Z, Explicit R
    // ==========================================
    for (int j = 0; j < NR - 1; j++){
        for (int k = 0; k < NZ - 1; k++){
            CN_d(k) = alpha_R(j) / dR(j) * (1.0 - 0.5 * dR(j + 1)/ Rc(j)) * temp(j - 1, k) + (1.0 - alpha_R(j) * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) * (1.0 - (dR(j + 1) - dR(j)) / 2.0 / Rc(j))) * temp(j, k) + alpha_R(j) / dR(j + 1) * (1.0 + 0.5 * dR(j) / Rc(j)) * temp(j + 1, k) - temp(j, k) * dT_tau(j, k) / 4.0;
        }

        for (int k = 0; k < NZ - 1; k++){
            CN_a(k) = - alpha_Z(k) / dZ(k);
            CN_b(k) = 1.0 + alpha_Z(k) * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k)) + dT_tau(j, k) / 4.0;
            CN_c(k) = - alpha_Z(k) / dZ(k + 1);
        }
        
        CN_a(-1) = 0.0; CN_b(-1) = 1.0; CN_c(-1) = -1.0; CN_d(-1) = 0.0;

        // 使用 subspan 截取当前维度大小 (NZ - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(1, NZ), s_b.subspan(1, NZ), s_c.subspan(1, NZ), s_d.subspan(1, NZ), s_x.subspan(1, NZ), s_cp.subspan(1, NZ), s_dp.subspan(1, NZ));
        
        CN_x(NZ - 1) = 0.0; // enforce absorbing boundary at Zmax
        
        for (int k = 0; k < NZ; k++) ndis(j, k) = CN_x(k) + src(j, k);
    }
    apply_boundary_conditions(ndis);
}



void solve_diffusion_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, double D, double Rg, double dT) {
    Field2D src_C;
    Field2D src_B;
    Field2D dT_tau;
    Field2D temp;
    initialize_disk_source(src_C);

    double beta = Rg / sqrt(Rg * Rg + rest_energy_p * rest_energy_p); // beta = v/c

    // 1. define ADI coefficients
    Array1D alpha_R(NR), alpha_Z(NZ);
    for (int i = 0; i < NR - 1; i++) alpha_R(i) = D * dT / (dR(i) + dR(i + 1));
    for (int i = 0; i < NZ - 1; i++) alpha_Z(i) = D * dT / (dZ(i) + dZ(i + 1));

    // 2. fragmentation and decay term
    for (int j = 0; j < NR; j++) {
        for (int k = 0; k < NZ; k++) {
            dT_tau(j, k) = dT * yr2s * (ndis_H(j, k) * sigma_C * beta * c + 1.0 / tau_stable);
        }
    }

    for (int t = 0; t < NT; t++) {
        CN_scheme_2nd_order_gc(ndis_C, src_C, ndis_H, dT_tau, temp, alpha_R, alpha_Z, beta, sigma_C, tau_stable, D, DT);
        for (int j = 0; j < NR; j++) {
            for (int k = 0; k < NZ; k++) {
                src_B(j, k) = ndis_C(j, k) * sigma_C2B * beta * ndis_H(j, k) * c * yr2s * DT;
            }
        }
        CN_scheme_2nd_order_gc(ndis_B, src_B, ndis_H, dT_tau, temp, alpha_R, alpha_Z, beta, sigma_B, tau_stable, D, DT);
    }
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

    // ----- z direction -----
    for (int i = 2; i < NR - 2; i++) {
        for (int j = 2; j < NZ - 2; j++) {
            double n_ip_z, flux_ip_z = 0.0; // i+1/2
            double n_im_z, flux_im_z = 0.0; // i-1/2
            double phi_z; // TVD Limiter
            
            // ----- i+1/2 -----
            double vZ_ip = (vZ(i,j) + vZ(i,j+1)) / 2.0;

            if (vZ_ip >= 0.0) {
                double ratio_ip_z = (temp(i,j) - temp(i,j-1)) / dZ(j)/ (temp(i,j+1) - temp(i,j) + 1e-12) * dZ(j + 1);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j) + 0.5 * phi_z * (temp(i,j+1) -  temp(i,j));
            } else {
                double ratio_ip_z = (temp(i,j+2) - temp(i,j+1)) / dZ(j + 2) / (temp(i,j+1) - temp(i,j) + 1e-12) * dZ(j + 1);
                Van_leer_limiter(phi_z, ratio_ip_z);
                n_ip_z = temp(i,j+1) - 0.5 * phi_z * (temp(i,j+1) - temp(i,j));
            }
            flux_ip_z = vZ_ip * n_ip_z;

            // ----- i-1/2 -----
            double vZ_im = (vZ(i,j-1) + vZ(i,j)) / 2.0;
            if (vZ_im >= 0.0) {
                double ratio_im_z = (temp(i,j-1) - temp(i,j-2)) / dZ(j - 1) / (temp(i,j) - temp(i,j-1) + 1e-12) * dZ(j);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j-1) + 0.5 * phi_z * (temp(i,j) - temp(i,j-1));
            } else {
                double ratio_im_z = (temp(i,j+1) - temp(i,j)) / dZ(j + 1) / (temp(i,j) - temp(i,j-1) + 1e-12) * dZ(j);
                Van_leer_limiter(phi_z, ratio_im_z);
                n_im_z = temp(i,j) - 0.5 * phi_z * (temp(i,j) - temp(i,j-1));
            }
            flux_im_z = vZ_im * n_im_z;
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
    Field2D src_B;
    Field2D temp;
    initialize_disk_source(src_C);
    for (int t = 0; t < NT; t++){
        advection_TVD_solution(ndis_C, temp, vR, vZ, DT);
        for (int j = 0; j < NR; j++){
            for (int k = 0; k < NZ; k++){
                src_B(j, k) = ndis_C(j, k) * sigma_C2B * ndis_H(j, k) * c * yr2s * DT;
            }
        }
        advection_TVD_solution(ndis_B, temp, vR, vZ, DT);
    }
}

void solve_transport_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& vR, Field2D& vZ, Field2D& ndis_H, double D, double Rg) {
    Field2D src_C;
    Field2D src_B;
    Field2D dT_tau;
    Field2D temp;
    initialize_disk_source(src_C);

    double beta = Rg / sqrt(Rg * Rg + rest_energy_p * rest_energy_p); // beta = v/c

    // first half step of Strang splitting
    advection_TVD_R(ndis_C, temp, vR, DT / 2.0);
    advection_TVD_R(ndis_B, temp, vR, DT / 2.0);
    advection_TVD_Z(ndis_C, temp, vZ, DT / 2.0);
    advection_TVD_Z(ndis_B, temp, vZ, DT / 2.0);

    // full step of Srang splitting, adjecent advection half steps are combined into one full step of advection; but L_R and L_Z do not commute, so the order of accuracy is not O(DT^2); swap the direction of advection every other time step can help to reduce the splitting error
    for (int t = 0; t < NT; t++) {
        CN_scheme_2nd_order(ndis_C, src_C, ndis_H, dT_tau, temp, beta, sigma_C, tau_stable, D, DT); // Crank-Nicolson + Peaceman–Rachford ADI is already second order accurate in time, and the spatial discretization is also symmetric in R, Z directions and second order accurate

        for (int j = 0; j < NR; j++) {
            for (int k = 0; k < NZ; k++) {
                src_B(j, k) = ndis_C(j, k) * sigma_C2B * beta * ndis_H(j, k) * c * yr2s * DT;
            }
        }

        CN_scheme_2nd_order(ndis_B, src_B, ndis_H, dT_tau, temp, beta, sigma_B, tau_stable, D, DT);

        if (t < NT - 1) {
            if (t % 2){
                advection_TVD_R(ndis_C, temp, vR, DT);
                advection_TVD_R(ndis_B, temp, vR, DT);
                advection_TVD_Z(ndis_C, temp, vZ, DT);
                advection_TVD_Z(ndis_B, temp, vZ, DT);
            }
            else {
                advection_TVD_Z(ndis_C, temp, vZ, DT);
                advection_TVD_Z(ndis_B, temp, vZ, DT);
                advection_TVD_R(ndis_C, temp, vR, DT);
                advection_TVD_R(ndis_B, temp, vR, DT);
            }
        }
    }
    advection_TVD_Z(ndis_C, temp, vZ, DT / 2.0);
    advection_TVD_Z(ndis_B, temp, vZ, DT / 2.0);
    advection_TVD_R(ndis_C, temp, vR, DT / 2.0);
    advection_TVD_R(ndis_B, temp, vR, DT / 2.0);

}

