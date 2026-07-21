#pragma once  // 现代写法，防止头文件被重复包含（等价于传统的 #ifndef/#define/#endif）

#include <vector>
#include <fstream>
#include <iostream>
#include <span>
// 常量声明放头文件，方便main.cpp和diffusion.cpp共享
constexpr int NR = 301;
constexpr int NZ = 151;
constexpr int NT = 1'000'00;
constexpr double DT = 1000.0;

class Field2D
{
public:
    int NR;
    int NZ;

    std::vector<double> data;

    Field2D(int nr, int nz)
        : NR(nr), NZ(nz), data(nr * nz) {}

    // 直接在类内部实现的函数，编译器会自动视为 inline，
    // 所以即使不显式写 inline 也是可以的。
    inline double& operator()(int i, int k)
    {
        return data[i * NZ + k];
    }

    inline const double& operator()(int i, int k) const
    {
        return data[i * NZ + k];
    }
};

// 只写函数原型（声明），不写函数体
void initialize_grids_log(void);
void initialize_disk_source(Field2D& ndis);
void initialize_to_0(Field2D& ndis);
void initialize_wind_velocity(Field2D& vR, Field2D& vZ);
void initialize_H(Field2D& ndis_H, const double nH);
void apply_boundary_conditions(Field2D& ndis);
void write_array_to_bin(const std::string& filename, Field2D& arr, const std::size_t size);
inline void thomas_solve( 
    std::span<const double> a, 
    std::span<const double> b, 
    std::span<const double> c, 
    std::span<const double> d, 
    std::span<double> x, 
    std::span<double> cp, 
    std::span<double> dp);

void CN_scheme_1st_order(Field2D& ndis, 
               const Field2D& src, 
               const Field2D& ndis_H, 
               Field2D& dT_tau, 
               Field2D& temp, 
               double beta, double sigma, double tau_decay, double D, double dT);

void CN_scheme_2nd_order(Field2D& ndis, 
               const Field2D& src, 
               const Field2D& ndis_H, 
               Field2D& dT_tau, 
               Field2D& temp, 
               double beta, double sigma, double tau_decay, double D, double dT);

void solve_diffusion_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& ndis_H, double D, double Rg);

void advection_TVD_R(Field2D& ndis, Field2D& temp, Field2D& vR, double dT);

void advection_TVD_Z(Field2D& ndis, Field2D& temp, Field2D& vZ, double dT);

inline void Van_leer_limiter(double& phi, double ratio);

void advection_TVD_solution(Field2D& ndis, Field2D& temp, Field2D& vR, Field2D& vZ, double dT);

void solve_transport_equation(Field2D& ndis_C, Field2D& ndis_B, Field2D& vR, Field2D& vZ, Field2D& ndis_H, double D, double Rg);