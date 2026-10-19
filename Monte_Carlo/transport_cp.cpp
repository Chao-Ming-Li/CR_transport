#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>

// Same instantaneous axisymmetric Gaussian burst as advection_grid:
// n(R,z,0) = exp(-((R-10)^2 + (z-10)^2)/9), positions in kpc.
// No energy evolution, interactions, decay, or additional injection.
int main() {
    constexpr int N_particles = 10000;
    constexpr int N_steps = 10000;
    constexpr double T_end = 1e8; // yr
    constexpr double dt = T_end / N_steps;
    constexpr double kpc_cm = 3.0856775814913673e21;
    constexpr double year_s = 31557600.0;
    constexpr double D = 1e29 * year_s / (kpc_cm * kpc_cm); // kpc^2/yr
    constexpr double v_z = 300.0 * 1e5 * year_s / kpc_cm; // kpc/yr
    constexpr double outer_boundary = 200.0; // kpc, absorbing R and z
    const double sigma = 3.0 / std::sqrt(2.0);
    const double diffusion_step = std::sqrt(2.0 * D * dt);
    const double pi = std::acos(-1.0);
    std::vector<double> x(N_particles), y(N_particles), z(N_particles);
    std::vector<int> alive(N_particles, 1);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N_particles; ++i) {
        // Per-particle streams make results independent of thread scheduling.
        std::seed_seq seed{12345, i};
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        std::normal_distribution<double> source(10.0, sigma);
        std::normal_distribution<double> normal(0.0, 1.0);
        double R;
        // Probability in dR is proportional to R*n(R), not just n(R).
        do {
            R = source(rng);
        } while (R < 0.0 || R >= outer_boundary || uniform(rng) >= R / outer_boundary);
        const double phi = 2.0 * pi * uniform(rng);
        x[i] = R * std::cos(phi);
        y[i] = R * std::sin(phi);
        do {
            z[i] = source(rng);
        } while (z[i] < 0.0 || z[i] >= outer_boundary);

        for (int step = 0; step < N_steps; ++step) {
            // Independent Gaussian increments give <dx_a^2> = 2*D*dt.
            // Cartesian diffusion includes the cylindrical radial geometry.
            x[i] += diffusion_step * normal(rng);
            y[i] += diffusion_step * normal(rng);
            z[i] = std::abs(z[i] + v_z * dt + diffusion_step * normal(rng));
            // abs implements the reflecting z=0 boundary with finite-step
            // error; reduce dt for boundary convergence. R=0 is an axis,
            // requiring no explicit reflection in Cartesian coordinates.
            if (std::hypot(x[i], y[i]) >= outer_boundary || z[i] >= outer_boundary) {
                alive[i] = 0;
                break;
            }
        }
    }

    // Preserve particles.bin contract: native float64 (N_surviving,3), pc.
    // Absorbed particles are omitted. Each surviving particle has weight 1/N.
    std::ofstream output("particles.bin", std::ios::binary | std::ios::trunc);
    if (!output) {
        std::cerr << "Cannot open particles.bin\n";
        return 1;
    }
    int survivors = 0;
    for (int i = 0; i < N_particles; ++i) {
        if (!alive[i])
            continue;
        const double position[] = {1000.0 * x[i], 1000.0 * y[i], 1000.0 * z[i]};
        output.write(reinterpret_cast<const char *>(position), sizeof(position));
        ++survivors;
    }
    output.close();
    if (!output) {
        std::cerr << "Failed to write particles.bin\n";
        return 1;
    }
    std::cout << "Saved " << survivors << "/" << N_particles
              << " particles at t=" << T_end << " yr; D=" << D
              << " kpc^2/yr, vz=" << v_z << " kpc/yr\n";
}
