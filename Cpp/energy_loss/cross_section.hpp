#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string_view>

// All physical constants, unit conversions and nuclear data of the project.
namespace physical_constants {

    // Use the Julian year consistently in transport and nuclear reactions.
    inline constexpr double seconds_per_year = 31557600.0;
    inline constexpr double cm_per_kpc = 3.0856775814913673e21;
    inline constexpr double cm_per_km = 1e5;
    inline constexpr double km_s_to_kpc_yr = cm_per_km * seconds_per_year / cm_per_kpc;
    inline constexpr double speed_of_light_cm_s = 3.0e10;
    inline constexpr double nucleon_rest_energy_GeV = 0.9315;

} // namespace physical_constants

// Nuclear data. Cross sections: cm^2; gas number density: cm^-3;
// reaction rates built from them: yr^-1.
namespace cross_section {

    enum class Isotope : std::size_t {
        O16,
        N15,
        N14,
        C13,
        C12,
        B11,
        B10,
        Be10,
        Be9,
        Be7,
        Proton, // Transport only: proton reaction cross sections are not supplied.
        Count
    }; // Count is not a nuclide, just a sentinel for array sizes.

    inline constexpr std::size_t isotope_count = static_cast<std::size_t>(Isotope::Count);

    struct Nuclide {
        std::string_view name;
        int charge;      // Z
        int mass_number; // A
    };

    inline constexpr std::array<Nuclide, isotope_count> nuclides{
        {{"O16", 8, 16},
         {"N15", 7, 15},
         {"N14", 7, 14},
         {"C13", 6, 13},
         {"C12", 6, 12},
         {"B11", 5, 11},
         {"B10", 5, 10},
         {"Be10", 4, 10},
         {"Be9", 4, 9},
         {"Be7", 4, 7},
         {"Proton", 1, 1}}};

    constexpr std::size_t index(Isotope isotope) {
        const auto i = static_cast<std::size_t>(isotope);
        if (i >= isotope_count)
            throw std::out_of_range("Unknown isotope");
        return i;
    }

    constexpr const Nuclide &nuclide(Isotope isotope) {
        return nuclides[index(isotope)];
    }

    // sigma[parent][daughter]. Diagonal: total inelastic destruction.
    // Off-diagonal: partial fragmentation parent -> daughter.
    // Original DRAGON2 values converted from mb to cm^2.
    // The added proton row/column are zero: no proton reactions are modeled.
    inline constexpr std::array<std::array<double, isotope_count>, isotope_count> sigma{
        {{{293.2e-27, 60.77e-27, 30.07e-27, 21.66e-27, 33.16e-27, 26.83e-27, 10.40e-27, 4.07e-27,
           3.35e-27, 8.70e-27}},
         {{0.0, 279.7e-27, 20.15e-27, 46.14e-27, 43.82e-27, 36.48e-27, 13.59e-27, 3.40e-27,
           7.98e-27, 5.34e-27}},
         {{0.0, 0.0, 265.9e-27, 16.26e-27, 52.00e-27, 29.22e-27, 10.16e-27, 2.70e-27, 2.71e-27,
           12.1e-27}},
         {{0.0, 0.0, 0.0, 255.7e-27, 56.73e-27, 41.65e-27, 4.70e-27, 6.55e-27, 7.52e-27, 4.70e-27}},
         {{0.0, 0.0, 0.0, 0.0, 237.1e-27, 52.41e-27, 18.08e-27, 4.05e-27, 7.32e-27, 9.48e-27}},
         {{0.0, 0.0, 0.0, 0.0, 0.0, 211.9e-27, 40.85e-27, 12.95e-27, 15.15e-27, 4.56e-27}},
         {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 197.1e-27, 0.0, 13.95e-27, 21.27e-27}},
         {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 197.0e-27, 36.37e-27, 36.94e-27}},
         {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 181.7e-27, 45.17e-27}},
         {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 150.0e-27}}}};

} // namespace cross_section
