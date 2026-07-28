#include <vector>
#include <fstream>
#include <iostream>
#include <span>
// 常量声明放头文件，方便main.cpp和diffusion.cpp共享
constexpr int NR = 300;
constexpr int NZ = 150;
constexpr int NT = 1;
constexpr double DT = 1000.0;

Array1D Rc(NR), Zc(NZ), Rf(NR + 1), Zf(NZ + 1), dR(NR), dZ(NZ); // Rf, Zf are face values of cells, while Rc, Zc are center values of cells

class Field2D
{
public:
    static constexpr int NG = 2;
    std::vector<double> data;

    // 直接使用全局常量 NR 和 NZ
    Field2D() 
        : data((NR + 2 * NG) * (NZ + 2 * NG)) {}

    inline double& operator()(int i, int k)
    {
        return data[(i + NG) * (NZ + 2 * NG) + k + NG]; // 这里明确使用的是全局常量
    }
    inline const double& operator()(int i, int k) const
    {
        return data[(i + NG) * (NZ + 2 * NG) + k + NG]; // 这里明确使用的是全局常量
    }
};

class Array1D
{
public:
    static constexpr int NG = 2; 
    std::vector<double> arr;

    explicit Array1D(int N)
        : arr(N + 2 * NG, 0.0) {}    

    inline double& operator()(int i) { return arr[i + NG]; }
    inline const double& operator()(int i) const { return arr[i + NG]; }   
    
    // Pointer to the start of memory (including ghost cells)
    double* data() noexcept { return arr.data(); }
    const double* data() const noexcept { return arr.data(); }

    // Total allocated size (including ghost cells)
    std::size_t size() const noexcept { return arr.size(); }

    // Optional helper: Pointer directly to the physical domain (skipping left ghost cells)
    double* interior_data() noexcept { return arr.data() + NG; }
};

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


void CN_scheme_2nd_order_gc(Field2D& ndis, 
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
            CN_d(j) = alpha_Z[k] / dZ(k) * ndis(j, k - 1) 
                    + (1.0 - alpha_Z[k] * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k))) * ndis(j, k) 
                    + alpha_Z[k] / dZ(k + 1) * ndis(j, k + 1) 
                    - ndis(j, k) * dT_tau(j, k) / 4.0;
        }
        
        for (int j = 0; j < NR - 1; j++){     
            CN_a(j) = - alpha_R[j] / dR(j) * (1.0 - 0.5 * dR(j + 1) / Rc(j));
            CN_b(j) = 1.0 + alpha_R[j] * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) * (1.0 - (dR(j + 1) - dR(j)) / 2.0 / Rc(j)) + dT_tau(j, k) / 4.0;
            CN_c(j) = - alpha_R[j] / dR(j + 1) * (1.0 + 0.5 * dR(j) / Rc(j));
        }
        
        CN_a(-1) = 0.0; CN_b(-1) = 1.0; CN_c(-1) = -1.0; CN_d(-1) = 0.0;

        // 使用 subspan 截取当前维度大小 (NR - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(-1, NR - 1), s_b.subspan(-1, NR - 1), s_c.subspan(-1, NR - 1), s_d.subspan(-1, NR - 1), s_x.subspan(-1, NR - 1), s_cp.subspan(-1, NR - 1), s_dp.subspan(-1, NR - 1));
        
        CN_x(NR - 1) = 0.0; // enforce absorbing boundary at Rmax
        
        for (int j = 0; j < NR; j++) temp(j, k) = CN_x(j);
    }
    apply_boundary_conditions(temp);

    // ==========================================
    // Step 2: Implicit Z, Explicit R
    // ==========================================
    for (int j = 0; j < NR - 1; j++){
        for (int k = 0; k < NZ - 1; k++){
            CN_d(k) = alpha_R[j] / dR(j) * (1.0 - 0.5 * dR(j + 1)/ Rc(j)) * temp(j - 1, k) + (1.0 - alpha_R[j] * (dR(j + 1) + dR(j)) / (dR(j + 1) * dR(j)) * (1.0 - (dR(j + 1) - dR(j)) / 2.0 / Rc(j))) * temp(j, k) + alpha_R[j] / dR(j + 1) * (1.0 + 0.5 * dR(j) / Rc(j)) * temp(j + 1, k) - temp(j, k) * dT_tau(j, k) / 4.0;
        }

        for (int k = 0; k < NZ - 1; k++){
            CN_a(k) = - alpha_Z[k] / dZ(k);
            CN_b(k) = 1.0 + alpha_Z[k] * (dZ(k + 1) + dZ(k)) / (dZ(k + 1) * dZ(k)) + dT_tau(j, k) / 4.0;
            CN_c(k) = - alpha_Z[k] / dZ(k + 1);
        }
        
        CN_a(-1) = 0.0; CN_b(-1) = 1.0; CN_c(-1) = -1.0; CN_d(-1) = 0.0;

        // 使用 subspan 截取当前维度大小 (NZ - 1) 的视窗丢给托马斯算法
        thomas_solve(s_a.subspan(-1, NZ - 1), s_b.subspan(-1, NZ - 1), s_c.subspan(-1, NZ - 1), s_d.subspan(-1, NZ - 1), s_x.subspan(-1, NZ - 1), s_cp.subspan(-1, NZ - 1), s_dp.subspan(-1, NZ - 1));
        
        CN_x(NZ - 1) = 0.0; // enforce absorbing boundary at Zmax
        
        for (int k = 0; k < NZ; k++) ndis(j, k) = CN_x(k) + src(j, k);
    }
    apply_boundary_conditions(ndis);
}


void solve_diffusion_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, double D, double Rg) {
    Field2D src_C;
    Field2D src_B;
    Field2D dT_tau;
    Field2D temp;
    initialize_disk_source(src_C);

    double beta = Rg / sqrt(Rg * Rg + rest_energy_p * rest_energy_p); // beta = v/c

    for (int t = 0; t < NT; t++) {
        CN_scheme_2nd_order_gc(ndis_C, src_C, ndis_H, dT_tau, temp, beta, sigma_C, tau_stable, D, DT);
        for (int j = 0; j < NR; j++) {
            for (int k = 0; k < NZ; k++) {
                src_B(j, k) = ndis_C(j, k) * sigma_C2B * beta * ndis_H(j, k) * c * yr2s * DT;
            }
        }
        CN_scheme_2nd_order_gc(ndis_B, src_B, ndis_H, dT_tau, temp, beta, sigma_B, tau_stable, D, DT);
    }
}
