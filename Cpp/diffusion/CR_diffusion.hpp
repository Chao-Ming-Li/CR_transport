#pragma once  // 现代写法，防止头文件被重复包含（等价于传统的 #ifndef/#define/#endif）

#include <vector>
#include <fstream>
#include <iostream>
#include <span>
// 常量声明放头文件，方便main.cpp和diffusion.cpp共享
constexpr int NR = 300;
constexpr int NZ = 150;
constexpr int NT = 10000;
constexpr double DT = 10000.0;

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

// 只写函数原型（声明），不写函数体
void initialize_grids_log(void);

void initialize_disk_source(Field2D& ndis);

void initialize_H(Field2D& ndis_H, const double nH);

void apply_boundary_conditions(Field2D& ndis);

void write_array_to_bin(const std::string& filename, Field2D& arr, const std::size_t size);

void solve_diffusion_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, double D, double Rg, double dT);
