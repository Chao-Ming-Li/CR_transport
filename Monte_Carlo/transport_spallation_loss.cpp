#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include "cross_section.hpp"
#define PI 3.14159265358979323846
#define SEED 12345
#define N_particles 210000
#define T_end 1e8 // total simulation time in years
#define c_pc_yr 0.30660139 // speed of light in pc/yr
#define c_cm_s 2.99792458e10 // speed of light in cm/s
#define m_p 0.938272 // proton mass in GeV/c^2
#define yr2s 31557600.0 // seconds in a year
#define R_EARTH 8000.0 // Earth's galactocentric radius in pc

struct EarthArrival {
    double energy; // total nuclear energy in GeV
    double species_id;
    double time_yr; // elapsed propagation time at the recorded step endpoint
};

struct VelocityField {
    double u_x; // pc/yr
    double u_y; // pc/yr
    double u_z; // pc/yr
    double divergence; // du_x/dx + du_y/dy + du_z/dz, yr^-1
};

VelocityField velocity_field(double x_pc, double y_pc, double z_pc, double time_yr){
    const double a = 1e-8; // vertical wind gradient in yr^-1
    (void)x_pc;
    (void)y_pc;
    (void)time_yr;
    return {0.0, 0.0, a * z_pc, a};
}

double mean_free_path(double E, cross_section::Isotope species){
    const auto& nucleus = cross_section::nuclide(species);
    const double mass = nucleus.mass_number * m_p; // nuclear rest energy in GeV
    const double rigidity = std::sqrt((E - mass) * (E + mass)) / nucleus.charge; // GV
    // lambda = 2 pc at 10 GV; E is total nuclear energy in GeV.
    return 2.0 * std::pow(rigidity / 10.0, 1.0 / 3.0);
}

double distance_to_earth_ring(double x_pc, double y_pc, double z_pc, double lambda){
    const double dR_pc = std::hypot(x_pc, y_pc) - R_EARTH;
    if (std::abs(dR_pc) >= lambda || std::abs(z_pc) >= lambda)
        return lambda;
    return std::hypot(dR_pc, z_pc);
}

double apply_adiabatic_loss(double E, cross_section::Isotope species,
                            double divergence, double dt){
    // dp/dt = -div(u)*p/3, with total nuclear E [GeV] and p [GeV/c].
    // Approximate the nuclear rest energy as A*m_p (binding energy neglected).
    // Freeze the local divergence over this substep. This is exact for constant
    // divergence; spatial/time-dependent fields require timestep convergence.
    // Negative divergence describes compression and produces an energy gain.
    const double mass = cross_section::nuclide(species).mass_number * m_p;
    const double p2 = (E - mass) * (E + mass);
    return std::sqrt(mass * mass + p2 * std::exp(-2.0 * divergence * dt / 3.0));
}

// Galactocentric Cartesian positions in pc; hydrogen nuclei density in cm^-3.
double hydrogen_density(double x_pc, double y_pc, double z_pc){
    const double R = std::hypot(x_pc, y_pc) / 1000.0;
    const double Z = std::abs(z_pc) / 1000.0;
    if (R >= 7.8 && R <= 8.2 && Z < 0.2)
        return 0.05;

    const double z_H2 = (R > 2.5) ? 0.059 : 0.015;
    const double X_CO = 0.9 * std::exp(R / 7.1) + 0.5 * std::exp(-R / 3.0);
    const double n_H2 = (150.0 * std::exp(-std::pow((R - 0.05) / 0.2, 2))
                       + 0.49 * std::exp(-std::pow((R - 5.3) / 2.2, 2))
                       + 0.24 * std::exp(-std::pow((R - 8.1) / 1.5, 2)))
                       * std::exp(-std::log(2.0) * std::pow(Z / z_H2, 2)) / X_CO;
    double z_HI;
    if (R < 2.5) z_HI = 0.045;
    else if (R < 8.5) z_HI = 0.115;
    else z_HI = 0.115 * std::exp((R - 8.5) / 6.7);
    double n0;
    if (R < 3.0)
        n0 = 0.57 * std::exp((R - 3.0) / 0.4)
             + 8.0 * std::exp(-std::pow(R / 0.2, 2));
    else if (R < 13.0) n0 = 0.57;
    else n0 = 0.57 * std::exp(-(R - 13.0) / 4.0);
    const double n_HI = n0 * std::exp(-std::log(2.0) * std::pow(Z / z_HI, 2));
    const double n_HII = 0.025 * std::exp(-std::pow(R / 20.0, 2) - Z / 1.0)
                        + 0.2 * std::exp(-std::pow((R - 4.0) / 2.0, 2) - Z / 0.15);
    return 2.0 * n_H2 + n_HI + n_HII;
}

std::vector<cross_section::Isotope> initialize_species(){
    // O16:N14:C12 = 10:1:10, with exact particle counts.
    static_assert(N_particles % 21 == 0, "N_particles must be divisible by 21 for the exact species ratio");
    constexpr int n_O16 = N_particles * 10 / 21;
    constexpr int n_N14 = N_particles / 21;
    std::vector<cross_section::Isotope> species(N_particles,cross_section::Isotope::C12);
    std::fill_n(species.begin(), n_O16, cross_section::Isotope::O16);
    std::fill_n(species.begin() + n_O16, n_N14, cross_section::Isotope::N14);

    // Avoid correlating species with the ordered initial energies, randomly rearranges all elements in species
    std::mt19937_64 species_rng(SEED);
    std::shuffle(species.begin(), species.end(), species_rng);
    return species;
}

void initialize_injection_positions(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z){
    // Uniform source in a cylinder: radius 10 kpc and total thickness 100 pc.
    // Positions are stored in pc. sqrt(U) gives a uniform areal density.
    constexpr double source_radius_pc = 10.0e3;
    constexpr double source_half_height_pc = 100.0;
    std::mt19937_64 position_rng(SEED + 1);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (std::size_t i = 0; i < x.size(); ++i) {
        const double radius = source_radius_pc * std::sqrt(uniform(position_rng));
        const double azimuth = 2.0 * PI * uniform(position_rng);
        x[i] = radius * std::cos(azimuth);
        y[i] = radius * std::sin(azimuth);
        z[i] = (2.0 * uniform(position_rng) - 1.0) * source_half_height_pc;
    }
}

// Single-residual fragmentation model; Count denotes destruction into untracked nuclei.
// Density and speed are held constant over this transport step. Multiple reactions
// are allowed, so a daughter can fragment again during the remaining time.
// Preserve kinetic energy per nucleon in the straight-ahead approximation.
// All species use the same rest energy per nucleon m_p, as in the cooling update.
void spallate(cross_section::Isotope& species, double& E, double nH, double dt,
              std::mt19937_64& rng){
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    while (dt > 0.0 && species != cross_section::Isotope::Count) {
        const auto parent = cross_section::index(species);
        const double A_parent = cross_section::nuclide(species).mass_number;
        const double kinetic_energy_per_nucleon = E / A_parent - m_p;
        const double sigma_total = cross_section::sigma[parent][parent];
        const double rate = nH * c_cm_s * sigma_total * yr2s; // yr^-1; v approximately c
        const double probability = -std::expm1(-rate * dt);
        if (uniform(rng) >= probability) return;

        // Sample the event time conditional on an interaction within dt, so that the remaining time is still available for further reactions.
        const double event_time = -std::log1p(-uniform(rng) * probability) / rate;
        dt -= event_time;
        const double draw = uniform(rng) * sigma_total;
        double cumulative = 0.0;
        species = cross_section::Isotope::Count;
        for (std::size_t daughter = 0; daughter < cross_section::isotope_count; ++daughter) {
            if (daughter == parent) continue;
            cumulative += cross_section::sigma[parent][daughter];
            if (draw < cumulative) {
                species = static_cast<cross_section::Isotope>(daughter);
                const double A_daughter = cross_section::nuclide(species).mass_number;
                E = A_daughter * (kinetic_energy_per_nucleon + m_p);
                // Equivalent to E_daughter = E_parent*A_daughter/A_parent.
                // Other reaction products carry the remaining nuclear energy.
                break;
            }
        }
        // The unused part of sigma_total removes the particle from the network.
    }
}


bool output2bin(const std::vector<double>& x, const std::vector<double>& y,
                const std::vector<double>& z, const std::vector<double>& E,
                const std::vector<cross_section::Isotope>& species){
    // Five native-endian doubles per record (40 bytes): x, y, z [pc],
    // Total nuclear E [GeV], species ID. IDs follow Isotope in cross_section.hpp:
    // 0=O16, 1=N15, 2=N14, 3=C13, 4=C12, 5=B11, 6=B10,
    // 7=Be10, 8=Be9, 9=Be7; -1=destroyed into untracked fragments.
    // Keep removed particles in the output to preserve initial particle order.
    // Their E is the last tracked parent's total energy, not a fragment energy.
    // Layout is unchanged, but legacy files stored E in GeV/n.
    std::ofstream output("transport_spallation_loss.bin", std::ios::binary | std::ios::trunc);
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double species_id = species[i] == cross_section::Isotope::Count? -1.0 : static_cast<double>(species[i]);
        const double par[] = {x[i], y[i], z[i], E[i], species_id};
        output.write(reinterpret_cast<const char*>(par), sizeof(par));
    }
    output.close();
    return true;
}

bool output_earth_arrivals(const std::vector<std::vector<EarthArrival>>& arrivals){
    // Three native-endian doubles per record (24 bytes): total nuclear energy
    // [GeV], species ID, elapsed propagation time [yr].
    std::ofstream output("earth_arrivals.bin", std::ios::binary | std::ios::trunc);

    std::size_t count = 0;
    for (const auto& particle_arrivals : arrivals) {
        for (const EarthArrival& arrival : particle_arrivals) {
            const double record[] = {arrival.energy, arrival.species_id, arrival.time_yr};
            output.write(reinterpret_cast<const char*>(record), sizeof(record));
            ++count;
        }
    }
    output.close();

    std::cout << "Recorded Earth arrivals: " << count << '\n';
    return true;
}

int main(){
    const double E_min = 10.0; // minimum injected total energy per nucleon in GeV/n
    const double E_max = 1e6; // maximum injected total energy per nucleon in GeV/n
    std::vector<double> x(N_particles);
    std::vector<double> y(N_particles);
    std::vector<double> z(N_particles);
    std::vector<double> E(N_particles); // total nuclear energy in GeV
    std::vector<double> t(N_particles); // propagation time in years
    std::vector<cross_section::Isotope> species = initialize_species();
    std::vector<std::vector<EarthArrival>> earth_arrivals(N_particles);

    initialize_injection_positions(x, y, z);

    // Preserve the injected log grid in energy per nucleon, including endpoints.
    // Store total nuclear energy using each particle's initial mass number.
    const double log_E_min = std::log10(E_min);
    const double log_E_step = (std::log10(E_max) - log_E_min) / (N_particles - 1);
    for (int i = 0; i < N_particles; ++i) {
        const double A = cross_section::nuclide(species[i]).mass_number;
        E[i] = A * std::pow(10.0, log_E_min + i * log_E_step); // dN/dE \propto E^-1
    }

    // Each particle owns its random stream, independent of thread scheduling.
    #pragma omp parallel for default(none) shared(x, y, z, E, t, species, earth_arrivals) schedule(dynamic, 1)
    for (int i = 0; i < N_particles; ++i) {
        std::seed_seq particle_seed{SEED, i, 1};
        std::mt19937_64 particle_rng(particle_seed);
        std::uniform_real_distribution<double> particle_uniform(0.0, 1.0);
        while (t[i] < T_end && species[i] != cross_section::Isotope::Count) {
            // Explicit step: evaluate gas density at the starting position.
            // const double nH = hydrogen_density(x[i], y[i], z[i]);
            const VelocityField wind = velocity_field(x[i], y[i], z[i], t[i]);
            const double nH = 0.1; // constant hydrogen density in cm^-3
            double cos_theta_v = 2.0 * particle_uniform(particle_rng) - 1.0;
            double sin_theta_v = std::sqrt(1.0 - cos_theta_v * cos_theta_v);
            double phi_v = 2.0 * PI * particle_uniform(particle_rng);
            double lambda = mean_free_path(E[i], species[i]); // pc
            double dt = std::min(lambda / 2.0 / c_pc_yr, T_end - t[i]); // time step in years
            const double step_length = 2.0 * c_pc_yr * dt; // step length = lambda, until the last step, which may be shorter than lambda.
            const double dx = step_length * sin_theta_v * std::cos(phi_v) + wind.u_x * dt;
            const double dy = step_length * sin_theta_v * std::sin(phi_v) + wind.u_y * dt;
            const double dz = step_length * cos_theta_v + wind.u_z * dt;

            t[i] += dt; // update propagation time in years
            x[i] += dx;
            y[i] += dy;
            z[i] += dz;
            E[i] = apply_adiabatic_loss(E[i], species[i], wind.divergence, dt);
            spallate(species[i], E[i], nH, dt, particle_rng);

            // Check the endpoint against the axisymmetric solar-ring region.
            // Every step ending inside the region adds one energy/species record.
            if (species[i] != cross_section::Isotope::Count) {
                const double lambda_end = mean_free_path(E[i], species[i]);
                if (distance_to_earth_ring(x[i], y[i], z[i], lambda_end) < lambda_end) {
                    earth_arrivals[i].push_back({E[i], static_cast<double>(species[i]), t[i]});
                }
            }
        }   
    }

    if (!output2bin(x, y, z, E, species)) return 1;
    if (!output_earth_arrivals(earth_arrivals)) return 1;

    return 0;
}
