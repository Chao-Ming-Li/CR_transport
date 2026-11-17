#include "cross_section.hpp"
#include "configuration.hpp"
#include "initialization.hpp"
#include "spatial_solver.hpp"
#include "momentum_solver.hpp"
#include "nuclear_solver.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

int run_model(int argc, char **argv) {
    // 1. Read all run parameters from the configuration file and convert the physical units.
    if (argc > 2)
        throw std::runtime_error("Usage: cr_solver [model.cfg]");
    const ModelParameters settings = read_model_parameters(argc == 2 ? argv[1] : "model.cfg");
    const std::string &mode = settings.mode;
    const int steps = settings.steps;
    const std::filesystem::path output_directory = settings.output_directory;
    const double dt = settings.dt_yr;
    const double D0 = settings.D0_cm2_per_s * physical_constants::seconds_per_year / (physical_constants::cm_per_kpc * physical_constants::cm_per_kpc);
    const double delta = settings.delta, vA = settings.vA_km_per_s;
    const double nH = settings.nH_cm3, helium_enhancement = settings.helium_enhancement;
    const double source_rate = settings.source_rate;
    const int nR = settings.nR, nz = settings.nz, nRg = settings.nRg;
    const double R_max = settings.R_max_kpc, z_max = settings.z_max_kpc;
    const double Rg_min = settings.Rg_min_GV, Rg_max = settings.Rg_max_GV;
    const cr_diffusion::DiffusionModel model(D0, delta, vA);

    // 2. Construct grids
    const SpatialGrid grid(AxisGrid::linear(nR, 0.0, R_max), AxisGrid::linear(nz, 0.0, z_max));
    const RgGrid Rg = RgGrid::logarithmic(nRg + 1, Rg_min, Rg_max);
    
    // 3. Select isotopes and their source weights.
    const std::vector<cross_section::Isotope> isotopes = select_isotopes(settings);
    const std::vector<double> &source_weights = settings.source_weights;

    // 4. Initialize zero density, the continuous source, gas and face velocities.
    const int count = static_cast<int>(isotopes.size());
    std::vector<Field3D> density(count, Field3D(grid, Rg));
    Field3D source_profile(grid, Rg);
    initialize_CR_source(source_profile, grid, Rg, settings.source_R_kpc, settings.source_z_kpc, settings.source_width_kpc, settings.source_index);
    Field2D gas(grid);
    initialize_gas(gas, nH, grid);
    Field2D vR(grid, Field2D::Location::RadialFace), vZ(grid, Field2D::Location::VerticalFace);
    initialize_wind_velocity(vR, vZ, settings.wind_km_per_s);

    const CFLNumbers cfl = calculate_cfl_numbers(settings, grid, Rg, vR, vZ, gas, model, isotopes);

    // Cache the fixed grids, coefficients and timestep in each solver.
    cr_advection::Solver advection(grid, vR, vZ, dt);
    cr_diffusion::Solver diffusion(grid, Rg, model, dt);
    cr_momentum::Solver momentum(grid, Rg, vR, vZ, model, dt, cr_momentum::Boundary::Outflow, cr_momentum::Boundary::Outflow);
    cr_nuclear::Solver nuclear(grid, Rg, gas, isotopes, dt, helium_enhancement);

    // 5. Evolve: source -> diffusion -> advection -> momentum -> nuclear reactions.
    // Sequential splitting and nuclear backward Euler give first-order time accuracy.
    // Spatial ghosts are refreshed by spatial solvers.
    const bool evolve_advection = settings.wind_km_per_s != 0.0;
    const bool evolve_momentum = evolve_advection || vA != 0.0;
    const bool evolve_nuclear = mode != "proton";
    const int nrigidity = Rg.size();
#ifdef _OPENMP
    const int source_threads = std::min(omp_get_max_threads(), nrigidity);
#endif
    for (int step = 0; step < steps; ++step) {
        for (int s = 0; s < count; ++s) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(static) num_threads(source_threads) if(source_threads > 1)
#endif
            for (int k = 0; k < nrigidity; ++k) {
                for (int i = 0; i < grid.nR(); ++i) {
                    for (int j = 0; j < grid.nz(); ++j) {
                        density[s](k, i, j) += dt * source_rate * source_weights[s] * source_profile(k, i, j);
                    }
                }
            }
            diffusion.advance(density[s]);
            if (evolve_advection) {
                advection.advance(density[s]);
            }
            if (evolve_momentum) {
                momentum.advance(density[s]);
            }
        }
        if (evolve_nuclear) {
            nuclear.advance(density);
        }
    }

    // 6. Write physical density distributions and run parameters.
    write_density_distributions(output_directory, isotopes, density);
    write_model_parameters(output_directory / "resolved.cfg", settings);
    write_run_parameters(output_directory / "settings.txt", settings, isotopes, D0, cfl);
    std::cout << "Completed " << count << " isotopes, t=" << steps * dt << " yr; distributions: " << output_directory << '\n';
    return 0;
}

int main(int argc, char **argv) {
    try {
        return run_model(argc, argv);
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
