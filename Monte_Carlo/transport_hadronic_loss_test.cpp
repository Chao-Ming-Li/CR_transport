#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

#define PI 3.14159265358979323846
#define SEED 12345
#define N_particles 10000
#define T_end 1e8 // total simulation time in years
#define c_pc_yr 0.30660139 // speed of light in pc/yr
#define c_cm_s 2.99792458e10 // speed of light in cm/s
#define m_p 0.938272 // proton mass in GeV/c^2
#define yr2s 31557600.0 // seconds in a year
struct Proton {
    double momentum; // in GeV/c   
    double x; 
    double y;
    double z;
};

double sigma_pp(double T){
    constexpr double T_threshold = 0.2797;
    if (T <= T_threshold)
        return 0.0;
    double sigma = (30.7 - 0.96 * std::log(T / T_threshold) + 0.18 * std::log(T / T_threshold) * std::log(T / T_threshold)) * std::pow(1.0 - std::pow(T_threshold / T, 1.9),3.0) * 1e-27; // in cm^2
    return sigma;
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

int main(){
    const double v_z = 3e-2;  // speed in units of pc / yr, equivalent to 30000 km/s
    const double E_min = 10.0; // minimum energy in GeV
    const double E_max = 1e6; // maximum energy in GeV
    std::vector<double> x(N_particles);
    std::vector<double> y(N_particles);
    std::vector<double> z(N_particles);
    std::vector<double> E(N_particles); // total energy in unit of GeV
    std::vector<double> t(N_particles); // propagation time in years

    // Random engine
    std::mt19937_64 rng(SEED);  // fixed seed -> reproducible
    // Uniform random number in [0, 1)
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (int i = 0; i < N_particles; ++i) {
        const double u = uniform(rng);
        // E[i] = 1.0 / ((1.0 - u) / E_min + u / E_max); // initialize energy to power-law distribution with -2 slope 
        E[i] = E_min * std::pow(E_max / E_min, u); // initialize energy to power-law distribution with -1 slope 
        }

    // Each particle owns its random stream, independent of thread scheduling.
    // Keep initialization above serial to preserve the initial energy sample.
    #pragma omp parallel for default(none) shared(x, y, z, E, t, v_z) schedule(static)
    for (int i = 0; i < N_particles; ++i) {
        std::seed_seq particle_seed{SEED, i, 1};
        std::mt19937_64 particle_rng(particle_seed);
        std::uniform_real_distribution<double> particle_uniform(0.0, 1.0);
        while (t[i] < T_end) {
            // Explicit step: evaluate gas density at the starting position.
            // const double nH = hydrogen_density(x[i], y[i], z[i]);
            const double nH = 1.0; // constant hydrogen density in cm^-3
            double cos_theta_v = 2.0 * particle_uniform(particle_rng) - 1.0;
            double sin_theta_v = std::sqrt(1.0 - cos_theta_v * cos_theta_v);
            double phi_v = 2.0 * PI * particle_uniform(particle_rng);
            double lambda = 2.0 * std::pow(E[i], 1./3); // mean free path
            double dt = std::min(lambda / 2.0 / c_pc_yr, T_end - t[i]); // time step in years
            t[i] += dt; // update propagation time in years
            x[i] += lambda * sin_theta_v * std::cos(phi_v);
            y[i] += lambda * sin_theta_v * std::sin(phi_v);
            z[i] += lambda * cos_theta_v + v_z * dt;  // add vertical motion
            // Compute the energy loss due to hadronic interactions
            double T = E[i] - m_p; // kinetic energy in GeV
            double dT_dt = -0.5 * nH * sigma_pp(T) * c_cm_s * T * yr2s; // energy loss rate in GeV/yr
            E[i] += dT_dt * dt; // update energy
        }   
    }

    // Raw native-endian doubles: x0, y0, z0 (pc), E0 (GeV), x1, ...
    std::ofstream output("particles_loss_t1e8_nH1.0_test.bin", std::ios::binary | std::ios::trunc);
    if (!output) {
        std::cerr << "Cannot open particles_loss_t1e8_nH1.0_test.bin for writing\n";
        return 1;
    }
    for (std::size_t i = 0; i < N_particles; ++i) {
        const double par[] = {x[i], y[i], z[i], E[i]};
        output.write(reinterpret_cast<const char*>(par), sizeof(par));
    }
    output.close();

}
