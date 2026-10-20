#include <iostream>
#include "CR_transport_solver.hpp"
#include "field.hpp"
#include "initialization.hpp"
#include <numbers>

#include <cmath>
#include <array>
#include <omp.h>

int main() {
    // End both axes exactly at 200 kpc (geometric_capped may overshoot).
    const auto make_axis = [] {
        const auto capped = AxisGrid::geometric_capped(0.0, 200.0, 0.25, 1.0, 1.03);
        std::vector<double> faces;
        for (int i = 0; i < capped.size(); ++i)
            faces.push_back(capped.face(i));
        faces.push_back(200.0);
        return AxisGrid::from_faces(faces);
    };
    const Grid2D grid(make_axis(), make_axis());

    // Save physical cell geometry for readers; binary snapshots include ghosts.
    std::ofstream geometry("transport_grid.csv", std::ios::trunc);
    if (!geometry)
        throw std::runtime_error("Cannot open transport_grid.csv");
    geometry.precision(17);
    geometry << "i,R_centroid_kpc,z_center_kpc,dR_kpc,dz_kpc\n";
    for (int i = 0; i < grid.nR(); ++i)
        geometry << i << ',' << grid.radial_centroid(i) << ',' << grid.z().center(i)
                 << ',' << grid.R().width(i) << ',' << grid.z().width(i) << '\n';
    geometry.close();
    if (!geometry)
        throw std::runtime_error("Failed to write transport_grid.csv");

    Field2D ndis_C(grid);
    Field2D ndis_H(grid);
    Field2D vR(grid, Field2D::Location::RadialFace);
    Field2D vZ(grid, Field2D::Location::VerticalFace);

    constexpr double kpc_cm = 3.0856775814913673e21;
    constexpr double year_s = 31557600.0;
    const double D = 1e29 * year_s / (kpc_cm * kpc_cm); // kpc^2/yr
    initialize_gas(ndis_H, 0.001, grid);
    initialize_wind_velocity(vR, vZ, grid);
    initialize_CR_source(ndis_C, grid);
    solve_transport_equation(ndis_C, ndis_H, vR, vZ, D, grid);

    return 0;
}
