#include "field.hpp"
#include "CR_advection.hpp"
#include "initialization.hpp"
#include <cassert>
#include <cmath>
#include <stdexcept>

int main()
{
    const Grid2D grid(AxisGrid::from_faces({28., 32., 38.}),
                      AxisGrid::from_faces({29., 31., 35., 39.}));
    Field2D source(grid), gas(grid), vR(grid, Field2D::Location::RadialFace), vZ(grid, Field2D::Location::VerticalFace);
    assert(source.location() == Field2D::Location::Center);
    assert(vR.nR() == 3 && vR.nz() == 3 && vR.data.size() == 7 * 7);
    assert(vZ.nR() == 2 && vZ.nz() == 4 && vZ.data.size() == 6 * 8);
    assert(vR.matches(grid, Field2D::Location::RadialFace));
    assert(!vR.matches(grid, Field2D::Location::Center));
    // Exercise all storage, including two ghosts beyond the last physical face.
    for (Field2D* field : {&vR, &vZ}) {
        int value = 0;
        for (int i = -Field2D::NG; i < field->nR() + Field2D::NG; ++i)
            for (int j = -Field2D::NG; j < field->nz() + Field2D::NG; ++j)
                (*field)(i, j) = ++value;
        for (std::size_t k = 0; k < field->data.size(); ++k)
            assert(field->data[k] == static_cast<double>(k + 1));
    }
    assert(source.nR() == 2 && source.nz() == 3);
    assert(source.data.size() == 6 * 7);
    int value = 0;
    for (int i = -Field2D::NG; i < source.nR() + Field2D::NG; ++i)
        for (int j = -Field2D::NG; j < source.nz() + Field2D::NG; ++j)
            source(i, j) = ++value;
    const Field2D& view = source;
    value = 0;
    for (int i = -Field2D::NG; i < source.nR() + Field2D::NG; ++i)
        for (int j = -Field2D::NG; j < source.nz() + Field2D::NG; ++j)
            assert(view(i, j) == ++value);
    initialize_CR_source(source, grid);
    assert(source(0, 0) == DT);
    assert(std::abs(source(1, 1) / (DT) - std::exp(-34. / 9.)) < 1e-14);
    initialize_gas(gas, 0.001, grid);
    assert(gas(1, 2) == 0.001 && gas(-1, 0) == 0.);
    initialize_wind_velocity(vR, vZ, grid);
    assert(vZ(1, 3) == 0.1 && vZ(1, 0) == 0.);
    // Matching storage dimensions alone do not make a center field a face field.
    Field2D wrong_location(3, 3);
    bool wrong_rejected = false;
    try { initialize_wind_velocity(wrong_location, vZ, grid); }
    catch (const std::invalid_argument&) { wrong_rejected = true; }
    assert(wrong_rejected);
    wrong_rejected = false;
    try { initialize_gas(vR, 1., grid); }
    catch (const std::invalid_argument&) { wrong_rejected = true; }
    assert(wrong_rejected);
    Field2D mismatch(3, 2);
    bool rejected = false;
    try { initialize_gas(mismatch, 1., grid); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    const Grid2D large(AxisGrid::linear(300, 0., 1.), AxisGrid::linear(400, 0., 1.));
    Field2D large_field(large);
    initialize_gas(large_field, 2., large);
    assert(large_field(299, 399) == 2.);
}
