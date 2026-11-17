#pragma once

#include "cross_section.hpp"
#include "field.hpp"
#include "initialization.hpp"

#include <filesystem>
#include <string>
#include <vector>

// Runtime values are read from a complete configuration file; units are part of names.
struct ModelParameters {
    std::string mode;
    std::string output_directory;
    int steps;
    int nR;
    int nz;
    int nRg;
    double dt_yr;
    double R_max_kpc;
    double z_max_kpc;
    double Rg_min_GV;
    double Rg_max_GV;
    double D0_cm2_per_s;
    double delta;
    double vA_km_per_s;
    double nH_cm3;
    double helium_enhancement;
    double source_rate;
    double wind_km_per_s;
    double source_R_kpc;
    double source_z_kpc;
    double source_width_kpc;
    double source_index;
    std::vector<double> source_weights;
};

// Dimensionless dt times the maximum discrete outgoing/diagonal loss rate.
struct CFLNumbers {
    double advection = 0.0;
    double advection_limit = 0.5; // Conservative PLM + SSPRK2 bound on the main uniform grid.
    double spatial_diffusion_R = 0.0;
    double spatial_diffusion_z = 0.0;
    double spatial_diffusion = 0.0;
    double momentum_diffusion = 0.0;
    double adiabatic = 0.0;
    double momentum = 0.0; // Combined adiabatic and rigidity-diffusion diagonal rate.
    double nuclear_loss = 0.0;
};

CFLNumbers calculate_cfl_numbers(const ModelParameters &p, const SpatialGrid &grid, const RgGrid &Rg,
                                const Field2D &vR, const Field2D &vZ, const Field2D &gas,
                                const cr_diffusion::DiffusionModel &model, const std::vector<cross_section::Isotope> &isotopes);

ModelParameters read_model_parameters(const std::filesystem::path &path);
void validate_model_parameters(const ModelParameters &p);
std::vector<cross_section::Isotope> select_isotopes(const ModelParameters &p);
void write_model_parameters(const std::filesystem::path &path, const ModelParameters &p);
void write_run_parameters(const std::filesystem::path &path, const ModelParameters &p, const std::vector<cross_section::Isotope> &isotopes, double D0, const CFLNumbers &cfl);

// Raw native double values in [rigidity][R][z] order, z contiguous; physical cells only.
void write_density_distributions(const std::filesystem::path &directory, const std::vector<cross_section::Isotope> &isotopes, const std::vector<Field3D> &density);
