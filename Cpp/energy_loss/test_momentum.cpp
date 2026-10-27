#include "momentum_solver.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
    void check(bool ok) {
        if (!ok) throw std::runtime_error("Momentum test failed");
    }
    double number(const Field3D &q, const RgGrid &Rg) {
        double sum = 0.0;
        for (int k = 0; k < Rg.size(); ++k) sum += q(k, 0, 0) * Rg.width(k);
        return sum;
    }
    double mean(const Field3D &q, const RgGrid &Rg) {
        double sum = 0.0;
        for (int k = 0; k < Rg.size(); ++k)
            sum += q(k, 0, 0) * Rg.width(k) * Rg.center(k);
        return sum / number(q, Rg);
    }
}

int main() {
    const SpatialGrid grid(AxisGrid::linear(2, 0.0, 2.0), AxisGrid::linear(2, 0.0, 2.0));
    const auto Rg = RgGrid::logarithmic(121, 1.0, 1e6);
    Field2D vr(grid, Field2D::Location::RadialFace);
    Field2D vz(grid, Field2D::Location::VerticalFace);
    std::vector<double> D(Rg.size() + 1, 0.0);
    Field3D q(grid, Rg);
    for (int k = 0; k < Rg.size(); ++k)
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                q(k, i, j) = std::pow(Rg.center(k), -1.0);
    const auto initial = q;
    cr_momentum::Solver identity(grid, Rg, vr, vz, D, 0.1);
    identity.advance(q);
    for (std::size_t l = 0; l < q.data().size(); ++l)
        check(q.data()[l] == initial.data()[l]);

    // vR=a*R, vz=a*z => cylindrical divergence=3*a, including R=0.
    const double expansion = 0.1;
    for (int i = 0; i < vr.nR(); ++i)
        for (int j = 0; j < vr.nz(); ++j) vr(i,j) = expansion * grid.R().face(i);
    for (int i = 0; i < vz.nR(); ++i)
        for (int j = 0; j < vz.nz(); ++j) vz(i,j) = expansion * grid.z().face(j);
    cr_momentum::Solver cooling(grid, Rg, vr, vz, D, 0.01);
    const double N0 = number(q, Rg), mean0 = mean(q, Rg);
    q(0, -2, -2) = 7.0;
    for (int t = 0; t < 100; ++t) cooling.advance(q);
    check(std::abs(number(q, Rg) / N0 - 1.0) < 1e-12);
    check(std::abs(mean(q, Rg) / mean0 - std::exp(-expansion)) < 0.01);
    check(q(0, -2, -2) == 7.0);
    for (double value : q.data()) check(value >= 0.0);

    // Pure diffusion leaves q/Rg^2 constant and conserves number.
    for (double &v : vr.data) v = 0.0;
    for (double &v : vz.data) v = 0.0;
    for (int k = 0; k <= Rg.size(); ++k) D[k] = 0.1 * Rg.face(k) * Rg.face(k);
    for (int k = 0; k < Rg.size(); ++k)
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j) q(k,i,j) = Rg.center(k)*Rg.center(k);
    cr_momentum::Solver diffusion(grid, Rg, vr, vz, D, 0.01);
    diffusion.advance(q);
    for (int k = 0; k < Rg.size(); ++k)
        check(std::abs(q(k,0,0)/(Rg.center(k)*Rg.center(k))-1.0) < 1e-11);
    q = initial;
    diffusion.advance(q);
    check(std::abs(number(q, Rg) / N0 - 1.0) < 1e-11);
    check(mean(q, Rg) > mean0);
    for (double value : q.data()) check(value >= 0.0);

    // Strong cooling at the low boundary: outflow loses particles.
    for (int i = 0; i < vz.nR(); ++i)
        for (int j = 0; j < vz.nz(); ++j) vz(i,j) = 3.0 * grid.z().face(j);
    std::fill(D.begin(), D.end(), 0.0);
    q = initial;
    cr_momentum::Solver outflow(grid, Rg, vr, vz, D, 0.01,
                               cr_momentum::Boundary::Outflow);
    for (int t = 0; t < 1000; ++t) outflow.advance(q);
    check(number(q, Rg) < N0);
    for (double value : q.data()) check(value >= 0.0);
    // Two-cell diffusion has one decaying mode with an exact exponential.
    // Compare at fixed spatial/rigidity resolution to isolate time accuracy.
    for (double &v : vz.data) v = 0.0;
    const auto small_Rg = RgGrid::logarithmic(3, 1.0, 4.0);
    const std::vector<double> small_D(3, 0.1);
    const double x0 = small_Rg.center(0), x1 = small_Rg.center(1);
    const double w0 = small_Rg.width(0), w1 = small_Rg.width(1);
    const double conductance = std::pow(small_Rg.face(1), 2) * small_D[1] / (x1 - x0);
    const double lambda = conductance * (1.0 / (x0*x0*w0) + 1.0 / (x1*x1*w1));
    const auto error = [&](int steps) {
        Field3D density(grid, small_Rg);
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j) {
                density(0,i,j) = x0*x0 + 0.5;
                density(1,i,j) = x1*x1 - 0.5*w0/w1;
            }
        const double initial_number = number(density, small_Rg);
        cr_momentum::Solver solver(grid, small_Rg, vr, vz, small_D, 1.0/steps);
        for (int t = 0; t < steps; ++t) solver.advance(density);
        check(std::abs(number(density, small_Rg) / initial_number - 1.0) < 1e-13);
        return std::abs(density(0,0,0) - (x0*x0 + 0.5*std::exp(-lambda)));
    };
    const double coarse_error = error(4), fine_error = error(8);
    check(coarse_error / fine_error > 3.9 && coarse_error / fine_error < 4.1);

    // An isolated outflow cell follows the scalar CN decay factor, at both ends.
    for (const double rate : {-0.1, 0.1}) {
        for (int i = 0; i < vz.nR(); ++i)
            for (int j = 0; j < vz.nz(); ++j)
                vz(i,j) = -3.0 * rate * grid.z().face(j);
        Field3D density(grid, small_Rg);
        const int k = rate < 0.0 ? 1 : 0;
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j) density(k,i,j) = 1.0;
        const std::vector<double> zero_D(3, 0.0);
        const double dt = 0.1;
        cr_momentum::Solver solver(grid, small_Rg, vr, vz, zero_D, dt,
                                   cr_momentum::Boundary::Outflow,
                                   cr_momentum::Boundary::Outflow);
        const double face = small_Rg.face(rate < 0.0 ? k : k + 1);
        const double decay = std::abs(rate) * face / small_Rg.width(k);
        solver.advance(density);
        check(std::abs(density(k,0,0) - (1.0 - 0.5*dt*decay)/(1.0 + 0.5*dt*decay)) < 1e-13);
    }
    std::cout << "Momentum: identity, cooling, conservation, diffusion equilibrium, heating, outflow, small-step positivity and CN second-order convergence passed.\n";
}
