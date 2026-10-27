#include "spatial_solver.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
    void check(bool ok) {
        if (!ok) throw std::runtime_error("Spatial solver test failed");
    }
}

int main() {
    const SpatialGrid grid(AxisGrid::geometric(9, 0.0, 3.0, 1.1),
                           AxisGrid::geometric(12, 0.0, 2.0, 1.12));
    const auto Rg = RgGrid::logarithmic(5, 1.0, 100.0);
    Field2D vr(grid, Field2D::Location::RadialFace);
    Field2D vz(grid, Field2D::Location::VerticalFace);
    // Closed face velocities isolate conservation, with both signs inside.
    for (int i = 1; i < grid.nR(); ++i)
        for (int j = 0; j < grid.nz(); ++j) vr(i,j) = 0.1*std::cos(i+j);
    for (int i = 0; i < grid.nR(); ++i)
        for (int j = 1; j < grid.nz(); ++j) vz(i,j) = 0.1*std::sin(i-j);

    for (auto reconstruction : {cr_advection::Reconstruction::Constant, cr_advection::Reconstruction::PLM})
        for (auto limiter : {cr_advection::Limiter::Minmod, cr_advection::Limiter::VanLeer, cr_advection::Limiter::MC})
            for (auto integrator : {cr_advection::TimeIntegrator::Euler, cr_advection::TimeIntegrator::SSPRK2}) {
                Field3D q(grid, Rg);
                std::vector<Field2D> reference;
                std::vector<double> initial_mass;
                for (int k = 0; k < Rg.size(); ++k) {
                    reference.emplace_back(grid);
                    double mass = 0.0;
                    for (int i = 0; i < grid.nR(); ++i)
                        for (int j = 0; j < grid.nz(); ++j) {
                            const double value = (k+1)*std::exp(-std::pow(grid.radial_centroid(i)-0.3*k, 2)
                                                               -std::pow(grid.z().center(j)-0.2*k, 2));
                            q(k,i,j) = reference[k](i,j) = value;
                            mass += value*grid.volume(i,j);
                        }
                    initial_mass.push_back(mass);
                }
                cr_advection::Options options;
                options.reconstruction = reconstruction;
                options.limiter = limiter;
                options.integrator = integrator;
                cr_advection::Solver advection(grid, vr, vz, 0.001, options);
                cr_diffusion::Solver diffusion(grid, 0.02, 0.001);
                cr_advection::Solver reference_advection(grid, vr, vz, 0.001, options);
                cr_diffusion::Solver reference_diffusion(grid, 0.02, 0.001);
                for (int step = 0; step < 20; ++step) {
                    advection.advance(q);
                    for (int k = 0; k < Rg.size(); ++k) {
                        reference_advection.advance(reference[k]);
                        auto values = q.spatial_plane(k);
                        check(std::equal(values.begin(), values.end(), reference[k].data.begin()));
                        if (step == 0) {
                            double mass = 0.0;
                            for (int i = 0; i < grid.nR(); ++i)
                                for (int j = 0; j < grid.nz(); ++j)
                                    mass += q(k,i,j)*grid.volume(i,j);
                            check(std::abs(mass/initial_mass[k]-1.0) < 1e-13);
                        }
                    }
                    diffusion.advance(q);
                    for (int k = 0; k < Rg.size(); ++k) {
                        reference_diffusion.advance(reference[k]);
                        auto values = q.spatial_plane(k);
                        check(std::equal(values.begin(), values.end(), reference[k].data.begin()));
                        check(q(k,-1,2) == q(k,0,2));
                        check(q(k,2,-1) == q(k,2,0));
                        check(q(k,grid.nR(),2) == 0.0);
                        check(q(k,2,grid.nz()) == 0.0);
                    }
                }
                const std::vector<double> before(q.data().begin(), q.data().end());
                cr_advection::Solver(grid, vr, vz, 0.0, options).advance(q);
                cr_diffusion::Solver(grid, 0.0, 0.001).advance(q);
                check(std::equal(q.data().begin(), q.data().end(), before.begin()));
            }
    std::cout << "Spatial: rigidity-plane equivalence, limiter/integrator combinations, conservation, boundaries and disabled operators passed.\n";
}
