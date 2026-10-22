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
    const auto Rg = RgGrid::logarithmic(121, 1.0, 1e4);
    Field2D vr(grid, Field2D::Location::RadialFace);
    Field2D vz(grid, Field2D::Location::VerticalFace);
    std::vector<double> D(Rg.size() + 1, 0.0);
    Field3D q(grid, Rg);
    for (int k = 0; k < Rg.size(); ++k)
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                q(k, i, j) = std::exp(-std::pow(std::log(Rg.center(k) / 100.0), 2) / 0.08);
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
    cr_momentum::Solver diffusion(grid, Rg, vr, vz, D, 10.0);
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
    cr_momentum::Solver outflow(grid, Rg, vr, vz, D, 10.0,
                               cr_momentum::Boundary::Outflow);
    outflow.advance(q);
    check(number(q, Rg) < N0);
    for (double value : q.data()) check(value >= 0.0);
    std::cout << "Momentum: identity, cooling, conservation, diffusion equilibrium, heating, outflow and positivity passed.\n";
}
