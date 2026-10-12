#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>

#define PI 3.14159265358979323846

int main(){
    const std::size_t N_particles = 10000;
    const int N_steps = 10000;
    double l_scatter = 10.0;  // fixed diffusion step in pc, not a sampled mean free path
    double dt = 16.0;  // time step in years
    double v_z = 3e-2;  // speed in units of pc / yr, equivalent to 30000 km/s
    std::vector<double> x(N_particles);
    std::vector<double> y(N_particles);
    std::vector<double> z(N_particles);

    // Random engine
    std::mt19937_64 rng(12345);  // fixed seed -> reproducible

    // Uniform random number in [0, 1)
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (int t = 0; t < N_steps; t++){
        for (std::size_t i = 0; i < N_particles; ++i) {
            double cos_theta_v = 2.0 * uniform(rng) - 1.0;
            double sin_theta_v = std::sqrt(1.0 - cos_theta_v * cos_theta_v);
            double phi_v = 2.0 * PI * uniform(rng);
            x[i] += l_scatter * sin_theta_v * std::cos(phi_v);
            y[i] += l_scatter * sin_theta_v * std::sin(phi_v);
            z[i] += l_scatter * cos_theta_v + v_z * dt;  // add vertical motion
        }   
    }

    // Raw native-endian doubles in pc: x0, y0, z0, x1, y1, z1, ...
    std::ofstream output("particles.bin", std::ios::binary | std::ios::trunc);
    if (!output) {
        std::cerr << "Cannot open particles.bin for writing\n";
        return 1;
    }
    for (std::size_t i = 0; i < N_particles; ++i) {
        const double position[] = {x[i], y[i], z[i]};
        output.write(reinterpret_cast<const char*>(position), sizeof(position));
    }
    output.close();
    if (!output) {
        std::cerr << "Failed to write particles.bin\n";
        return 1;
    }
    std::cout << "Saved " << N_particles << " final positions to particles.bin"
              << " at t = " << N_steps * dt << " yr\n";
}
