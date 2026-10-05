// Compare the actual production routines, without changing their integrators.
#include "CR_advection.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

constexpr double speed = 0.1, duration = 500.0;
double exact(double l, double r, double shift, bool hat) {
    if (hat)
        return std::max(0.0, std::min(r, 35.0 + shift) - std::max(l, 25.0 + shift)) / (r - l);
    // erfc avoids cancellation in the Gaussian tails.
    l = (l - 30.0 - shift) / 3.0;
    r = (r - 30.0 - shift) / 3.0;
    double integral = l >= 0   ? std::erfc(l) - std::erfc(r)
                      : r <= 0 ? std::erfc(-r) - std::erfc(-l)
                               : std::erf(r) - std::erf(l);
    return std::sqrt(std::acos(-1.0)) * integral / (2 * (r - l));
}
struct Metrics {
    double mass = 0, mean = 0, variance = 0, peak = 0, low = 1e300, tv = 0, error = 0;
};
Metrics measure(const Field2D &n, const Grid2D &g, double shift, bool hat) {
    Metrics m;
    double refmass = 0;
    for (int j = 0; j < g.nz(); ++j) {
        double q = n(0, j), z = g.z().center(j), dz = g.z().width(j);
        double ref = exact(g.z().face(j), g.z().face(j + 1), shift, hat);
        m.mass += q * dz;
        m.mean += q * z * dz;
        m.variance += q * (z * z + dz * dz / 12) * dz;
        m.peak = std::max(m.peak, q);
        m.low = std::min(m.low, q);
        m.error += std::abs(q - ref) * dz;
        refmass += ref * dz;
        m.tv += std::abs(q - (j ? n(0, j - 1) : 0.0));
    }
    m.tv += std::abs(n(0, g.nz() - 1));
    m.mean /= m.mass;
    m.variance = m.variance / m.mass - m.mean * m.mean;
    m.error /= refmass;
    return m;
}
int main() {
    std::ofstream profiles("ssprk2_profiles.csv"), metrics("ssprk2_metrics.csv");
    profiles << std::setprecision(17) << "profile,nz,dt,scheme,z,initial,exact,numerical\n";
    metrics << std::setprecision(12)
            << "profile,nz,dt,CFL,scheme,mass_rel_change,mean_error,variance_change,D_effective,"
               "peak,min,TV_ratio,relative_L1,status\n";
    for (bool hat : {false, true})
        for (int nz : {200, 400})
            for (double dt : {0.1, 1.0, 2.0, 5.0}) {
                if (speed * dt / (200.0 / nz) > 0.5)
                    continue;
                Grid2D g(AxisGrid::linear(2, 0, 1), AxisGrid::linear(nz, 0, 200));
                Field2D initial(g), vr(g, Field2D::Location::RadialFace),
                    vz(g, Field2D::Location::VerticalFace);
                for (int i = 0; i < g.nR(); ++i) {
                    for (int j = 0; j < nz; ++j)
                        initial(i, j) = exact(g.z().face(j), g.z().face(j + 1), 0, hat);
                    for (int j = 1; j <= nz; ++j)
                        vz(i, j) = speed;
                }
                auto before = measure(initial, g, 0, hat);
                for (std::string scheme : {"vanLeer_FE", "vanLeer_SSPRK2"}) {
                    Field2D n = initial, temp(g);
                    std::string status = "ok";
                    try {
                        int steps = std::lround(duration / dt);
                        for (int k = 0; k < steps; ++k) {
                            apply_boundary_conditions(n, g);
                            if (scheme == "vanLeer_FE")
                                advection_TVD(n, temp, vr, vz, dt, g);
                            else
                                advance_TVD_SSPRK2(n, temp, vr, vz, dt, g);
                        }
                    } catch (const std::exception &e) {
                        status = e.what();
                    }
                    metrics << (hat ? "top_hat" : "gaussian") << ',' << nz << ',' << dt << ','
                            << speed * dt / (200.0 / nz) << ',' << scheme << ',';
                    if (status != "ok") {
                        metrics << "nan,nan,nan,nan,nan,nan,nan,nan," << status << '\n';
                        continue;
                    }
                    auto after = measure(n, g, speed * duration, hat);
                    metrics << after.mass / before.mass - 1 << ','
                            << after.mean - before.mean - speed * duration << ','
                            << after.variance - before.variance << ','
                            << (after.variance - before.variance) / (2 * duration) << ','
                            << after.peak << ',' << after.low << ',' << after.tv / before.tv << ','
                            << after.error << ",ok\n";
                    for (int j = 0; j < nz; ++j)
                        profiles << (hat ? "top_hat" : "gaussian") << ',' << nz << ',' << dt << ','
                                 << scheme << ',' << g.z().center(j) << ',' << initial(0, j) << ','
                                 << exact(g.z().face(j), g.z().face(j + 1), speed * duration, hat)
                                 << ',' << n(0, j) << '\n';
                }
            }
    std::cout << "Wrote ssprk2_metrics.csv and ssprk2_profiles.csv\n";
}
