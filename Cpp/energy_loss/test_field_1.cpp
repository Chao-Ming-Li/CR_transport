#include "field.hpp"
#include <iostream>


int main() {
    int NG = 2;
    const Grid2D grid(AxisGrid::geometric_capped(0.0, 30.0, 0.25, 1.0, 1.03),
                      AxisGrid::linear(5, 0.0, 5.0));

    // N cells have N+1 physical faces, including the upper boundary.
    for (int i = -NG; i <= grid.nR() + NG; ++i) {
        std::cout << "R[" << i << "]: face = " << grid.R().face(i) << '\n';
    }
    for (int i = -NG; i <= grid.nz() + NG; ++i) {
        std::cout << "Z[" << i << "]: face = " << grid.z().face(i) << '\n';
    }

    const auto rigidity = RigidityGrid::logarithmic(61, 1.0, 1e6);
    for (int k = 0; k < rigidity.size(); ++k) {
        std::cout << "Rigidity[" << k << "]: face = " << rigidity.face(k)
                  << ", center = " << rigidity.center(k) << '\n';
    }
    std::cout << "Upper rigidity face = " << rigidity.face(rigidity.size()) << '\n';

}
