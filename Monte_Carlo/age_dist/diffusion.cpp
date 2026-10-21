#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

#define PI 3.14159265358979323846
#define SEED 12345
#define N_particles 10000
#define T_end 1e9 // total simulation time in years
#define c_pc_yr 0.30660139 // speed of light in pc/yr
#define c_cm_s 2.99792458e10 // speed of light in cm/s
#define m_p 0.938272 // proton mass in GeV/c^2
#define yr2s 31557600.0 // seconds in a year
#define R_EARTH 8000.0 // solar-ring galactocentric radius in pc

struct EarthArrival {
    double energy; // total proton energy in GeV
    double time_yr; // elapsed propagation time at the step endpoint
};

// Same axisymmetric solar-ring selection as transport_spallation_loss.cpp.
// The tube radius is the mean free path at the recorded endpoint energy.
double distance_to_earth_ring(double x_pc, double y_pc, double z_pc, double lambda){
    const double dR_pc = std::hypot(x_pc, y_pc) - R_EARTH;
    if (std::abs(dR_pc) >= lambda || std::abs(z_pc) >= lambda)
        return lambda;
    return std::hypot(dR_pc, z_pc);
}


bool output_earth_arrivals(const std::vector<std::vector<EarthArrival>>& arrivals){
    // Two native-endian doubles per record (16 bytes): E [GeV], time [yr].
    // Each step ending inside the ring contributes a record, including repeat
    // visits. Records are grouped by particle, not sorted by arrival time.
    std::ofstream output("proton_earth_t1e9_D0.5.bin", std::ios::binary | std::ios::trunc);
    for (const auto& particle_arrivals : arrivals) {
        for (const EarthArrival& arrival : particle_arrivals) {
            const double record[] = {arrival.energy, arrival.time_yr};
            output.write(reinterpret_cast<const char*>(record), sizeof(record));
        }
    }
    output.close();
    return true;
}


double mean_free_path(double E){
    // Proton charge Z=1; E is total proton energy in GeV.
    const double rigidity = std::sqrt((E - m_p) * (E + m_p)); // GV
    // Same normalization as transport_spallation_loss.cpp: 2 pc at 10 GV.
    return 2.0 * std::pow(rigidity / 10.0, 1.0 / 2.0);
}


int main(){
    const double E_min = 10.0; // minimum energy in GeV
    const double E_max = 1e7; // maximum energy in GeV
    // const double nH = 0.1; // constant hydrogen density in cm^-3

    std::vector<double> x(N_particles);
    std::vector<double> y(N_particles);
    std::vector<double> z(N_particles);
    std::vector<double> E(N_particles); // total energy in unit of GeV
    std::vector<double> t(N_particles); // propagation time in years

    std::vector<std::vector<EarthArrival>> earth_arrivals(N_particles);

    // Random engine
    std::mt19937_64 rng(SEED);  // fixed seed -> reproducible
    // Uniform random number in [0, 1)
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    const double log_E_min = std::log10(E_min);
    const double log_E_step = (std::log10(E_max) - log_E_min) / (N_particles - 1);
    for (int i = 0; i < N_particles; ++i) {
        E[i] = std::pow(10.0, log_E_min + i * log_E_step); // dN/dE \propto E^-1
        // E[i] = 100000.0;
    }

    // Each particle owns its random stream, independent of thread scheduling.
    // Keep initialization above serial to preserve the initial energy sample.
    #pragma omp parallel for default(none) shared(x, y, z, E, t, earth_arrivals) schedule(dynamic, 1)
    for (int i = 0; i < N_particles; ++i) {
        std::seed_seq particle_seed{SEED, i, 1};
        std::mt19937_64 particle_rng(particle_seed);
        std::uniform_real_distribution<double> particle_uniform(0.0, 1.0);
        while (t[i] < T_end) {
            // Explicit step: evaluate gas density at the starting position.
            // const double nH = hydrogen_density(x[i], y[i], z[i]);
            double cos_theta_v = 2.0 * particle_uniform(particle_rng) - 1.0;
            double sin_theta_v = std::sqrt(1.0 - cos_theta_v * cos_theta_v);
            double phi_v = 2.0 * PI * particle_uniform(particle_rng);
            double lambda = mean_free_path(E[i]); // pc
            const double dt_full = lambda / 2.0 / c_pc_yr;
            double dt = std::min(dt_full, T_end - t[i]); // time step in years
            // Preserve <dr^2> = 6*D*dt when the final step is shortened.
            const double step_length = lambda * std::sqrt(dt / dt_full);
            t[i] += dt; // update propagation time in years
            x[i] += step_length * sin_theta_v * std::cos(phi_v);
            y[i] += step_length * sin_theta_v * std::sin(phi_v);
            z[i] += step_length * cos_theta_v;

            // Record the cooled energy and elapsed time at the step endpoint.
            if (distance_to_earth_ring(x[i], y[i], z[i], lambda) < lambda) {
                earth_arrivals[i].push_back({E[i], t[i]});
            }
        }   
    }

    if (!output_earth_arrivals(earth_arrivals)) return 1;
    return 0;
}
