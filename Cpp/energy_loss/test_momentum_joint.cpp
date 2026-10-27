#include "momentum_solver.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

int main() {
    // Rigidity [GV], time [yr]; div(v)=3*a and Dpp=D on all faces.
    const int cells = 512;
    const double a = 0.01, D = 1.0, dt = 0.001;
    const double R0 = 20.0, sigma = 2.0;
    const SpatialGrid grid(AxisGrid::linear(2, 0.0, 2.0), AxisGrid::linear(2, 0.0, 2.0));
    const auto Rg = RgGrid::logarithmic(cells + 1, 0.1, 200.0);
    Field2D vr(grid, Field2D::Location::RadialFace);
    Field2D vz(grid, Field2D::Location::VerticalFace);
    std::vector<double> Dpp(Rg.size() + 1, D);
    Field3D q(grid, Rg);

    for (int k = 0; k < Rg.size(); ++k)
        for (int i = 0; i < grid.nR(); ++i)
            for (int j = 0; j < grid.nz(); ++j)
                q(k,i,j) = 0.5 * (std::erf((Rg.face(k+1)-R0)/(std::sqrt(2.0)*sigma))
                                - std::erf((Rg.face(k)-R0)/(std::sqrt(2.0)*sigma))) / Rg.width(k);
    const auto initial = q;
    for (int i = 0; i < vr.nR(); ++i)
        for (int j = 0; j < vr.nz(); ++j) vr(i,j) = a * grid.R().face(i);
    for (int i = 0; i < vz.nR(); ++i)
        for (int j = 0; j < vz.nz(); ++j) vz(i,j) = a * grid.z().face(j);

    cr_momentum::Solver joint(grid, Rg, vr, vz, Dpp, dt,
                             cr_momentum::Boundary::Outflow,
                             cr_momentum::Boundary::Outflow);
    int step = 0;
    for (int time : {1, 5, 10, 50, 100}) {
        const int target = static_cast<int>(std::round(time/dt));
        for (; step < target; ++step) joint.advance(q);
        // Five native float64 columns: lower face, upper face, center,
        // initial cell-average dn/dRg, evolved cell-average dn/dRg.
        std::ofstream output("momentum_joint_t" + std::to_string(time) + ".bin",
                             std::ios::binary | std::ios::trunc);
        for (int k = 0; k < Rg.size(); ++k) {
            const double record[] = {Rg.face(k), Rg.face(k+1), Rg.center(k),
                                     initial(k,0,0), q(k,0,0)};
            output.write(reinterpret_cast<const char *>(record), sizeof(record));
        }
    }
}
