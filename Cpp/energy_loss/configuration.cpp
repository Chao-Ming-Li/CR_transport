#include "cross_section.hpp"
#include "configuration.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {

    std::string trim(const std::string &text) {
        const auto first = text.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            return "";
        return text.substr(first, text.find_last_not_of(" \t\r\n") - first + 1);
    }

} // namespace

ModelParameters read_model_parameters(const std::filesystem::path &path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open configuration: " + path.string());

    ModelParameters p{};
    std::set<std::string> seen;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) { // automatically starts at line_number = 0
        ++line_number;
        line = trim(line.substr(0, line.find('#')));
        if (line.empty())
            continue;
        const auto separator = line.find('='); // the return type is std::string::size_type, which is unsigned.
        if (separator == std::string::npos)
            throw std::runtime_error("Expected key = value at line " + std::to_string(line_number));
        const std::string key = trim(line.substr(0, separator)); // line.substr(a, b) means start at a and take b characters
        const std::string value = trim(line.substr(separator + 1));
        if (value.empty() || !seen.insert(key).second)
            throw std::runtime_error("Empty value or duplicate key: " + key);

        // The stream try to read the value as the type defined for key values. If it succeeds, it will store it in the corresponding member variable. If it fails, it will set the failbit of the stream.
        std::istringstream number(value);
        if (key == "steps") number >> p.steps; 
        else if (key == "nR") number >> p.nR;
        else if (key == "nz") number >> p.nz;
        else if (key == "nRg") number >> p.nRg;
        else if (key == "dt_yr") number >> p.dt_yr;
        else if (key == "R_max_kpc") number >> p.R_max_kpc;
        else if (key == "z_max_kpc") number >> p.z_max_kpc;
        else if (key == "Rg_min_GV") number >> p.Rg_min_GV;
        else if (key == "Rg_max_GV") number >> p.Rg_max_GV;
        else if (key == "D0_cm2_per_s") number >> p.D0_cm2_per_s;
        else if (key == "delta") number >> p.delta;
        else if (key == "vA_km_per_s") number >> p.vA_km_per_s;
        else if (key == "nH_cm3") number >> p.nH_cm3;
        else if (key == "helium_enhancement") number >> p.helium_enhancement;
        else if (key == "source_rate") number >> p.source_rate;
        else if (key == "wind_km_per_s") number >> p.wind_km_per_s;
        else if (key == "source_R_kpc") number >> p.source_R_kpc;
        else if (key == "source_z_kpc") number >> p.source_z_kpc;
        else if (key == "source_width_kpc") number >> p.source_width_kpc;
        else if (key == "source_index") number >> p.source_index;
        else if (key == "mode") {
            p.mode = value;
            continue;
        }
        else if (key == "output_directory") {
            p.output_directory = value;
            continue;
        }
        else if (key == "source_weights") {
            std::istringstream values(value);
            std::string entry;
            if (value.back() == ',')
                throw std::runtime_error("Empty source weight");
            while (std::getline(values, entry, ',')) {
                std::istringstream weight(entry);
                double x;
                if (!(weight >> x) || !(weight >> std::ws).eof() || !std::isfinite(x) || x < 0.0)
                    throw std::runtime_error("Invalid source weight: " + entry);
                p.source_weights.push_back(x);
            }
            continue;
        } else {
            throw std::runtime_error("Unknown configuration key: " + key);
        }
        if (!number || !(number >> std::ws).eof())
            throw std::runtime_error("Invalid value for " + key + ": " + value);
    }
    if (input.bad())
        throw std::runtime_error("Failed reading configuration: " + path.string());
    for (const char *key : {
        "mode", "output_directory", "steps", "dt_yr", "source_weights",
        "nR", "nz", "nRg", "R_max_kpc", "z_max_kpc", "Rg_min_GV", "Rg_max_GV",
        "D0_cm2_per_s", "delta", "vA_km_per_s", "wind_km_per_s",
        "nH_cm3", "helium_enhancement", "source_rate", "source_R_kpc",
        "source_z_kpc", "source_width_kpc", "source_index"
    }) {
        if (!seen.contains(key))
            throw std::runtime_error("Missing configuration key: " + std::string(key));
    }
    validate_model_parameters(p);
    return p;
}

void validate_model_parameters(const ModelParameters &p) {
    // Setup-only checks: prevent invalid grids and singular diffusion coefficients.
    if (p.steps < 0 || p.nR < 2 || p.nz < 2 || p.nRg < 1 || p.nRg == std::numeric_limits<int>::max())
        throw std::runtime_error("Require steps >= 0, nR,nz >= 2 and 1 <= nRg < INT_MAX");
    if (!std::isfinite(p.dt_yr)) throw std::runtime_error("Nonfinite parameter: dt_yr");
    if (!std::isfinite(p.R_max_kpc)) throw std::runtime_error("Nonfinite parameter: R_max_kpc");
    if (!std::isfinite(p.z_max_kpc)) throw std::runtime_error("Nonfinite parameter: z_max_kpc");
    if (!std::isfinite(p.Rg_min_GV)) throw std::runtime_error("Nonfinite parameter: Rg_min_GV");
    if (!std::isfinite(p.Rg_max_GV)) throw std::runtime_error("Nonfinite parameter: Rg_max_GV");
    if (!std::isfinite(p.D0_cm2_per_s)) throw std::runtime_error("Nonfinite parameter: D0_cm2_per_s");
    if (!std::isfinite(p.delta)) throw std::runtime_error("Nonfinite parameter: delta");
    if (!std::isfinite(p.vA_km_per_s)) throw std::runtime_error("Nonfinite parameter: vA_km_per_s");
    if (!std::isfinite(p.nH_cm3)) throw std::runtime_error("Nonfinite parameter: nH_cm3");
    if (!std::isfinite(p.helium_enhancement)) throw std::runtime_error("Nonfinite parameter: helium_enhancement");
    if (!std::isfinite(p.source_rate)) throw std::runtime_error("Nonfinite parameter: source_rate");
    if (!std::isfinite(p.wind_km_per_s)) throw std::runtime_error("Nonfinite parameter: wind_km_per_s");
    if (!std::isfinite(p.source_R_kpc)) throw std::runtime_error("Nonfinite parameter: source_R_kpc");
    if (!std::isfinite(p.source_z_kpc)) throw std::runtime_error("Nonfinite parameter: source_z_kpc");
    if (!std::isfinite(p.source_width_kpc)) throw std::runtime_error("Nonfinite parameter: source_width_kpc");
    if (!std::isfinite(p.source_index)) throw std::runtime_error("Nonfinite parameter: source_index");
    if (p.dt_yr < 0.0 || p.R_max_kpc <= 0.0 || p.z_max_kpc <= 0.0 ||
        p.Rg_min_GV <= 0.0 || p.Rg_max_GV <= p.Rg_min_GV)
        throw std::runtime_error("Invalid timestep or grid bounds");
    if (p.D0_cm2_per_s <= 0.0 || p.delta <= 0.0 || p.delta >= 2.0 || p.vA_km_per_s < 0.0)
        throw std::runtime_error("Require D0 > 0, 0 < delta < 2 and vA >= 0");
    if (p.nH_cm3 < 0.0 || p.helium_enhancement < 0.0 || p.source_rate < 0.0 || p.source_width_kpc <= 0.0)
        throw std::runtime_error("Require nonnegative gas/source rates and positive source width");
    if (p.mode.empty() || p.mode.back() == ',' || p.output_directory.empty())
        throw std::runtime_error("Empty mode/output directory or trailing comma in mode");
}

void write_model_parameters(const std::filesystem::path &path, const ModelParameters &p) {
    std::ofstream output(path);
    output << std::setprecision(17);
    output << "mode = " << p.mode << '\n';
    output << "output_directory = " << p.output_directory << '\n';
    output << "steps = " << p.steps << '\n';
    output << "nR = " << p.nR << '\n';
    output << "nz = " << p.nz << '\n';
    output << "nRg = " << p.nRg << '\n';
    output << "dt_yr = " << p.dt_yr << '\n';
    output << "R_max_kpc = " << p.R_max_kpc << '\n';
    output << "z_max_kpc = " << p.z_max_kpc << '\n';
    output << "Rg_min_GV = " << p.Rg_min_GV << '\n';
    output << "Rg_max_GV = " << p.Rg_max_GV << '\n';
    output << "D0_cm2_per_s = " << p.D0_cm2_per_s << '\n';
    output << "delta = " << p.delta << '\n';
    output << "vA_km_per_s = " << p.vA_km_per_s << '\n';
    output << "nH_cm3 = " << p.nH_cm3 << '\n';
    output << "helium_enhancement = " << p.helium_enhancement << '\n';
    output << "source_rate = " << p.source_rate << '\n';
    output << "wind_km_per_s = " << p.wind_km_per_s << '\n';
    output << "source_R_kpc = " << p.source_R_kpc << '\n';
    output << "source_z_kpc = " << p.source_z_kpc << '\n';
    output << "source_width_kpc = " << p.source_width_kpc << '\n';
    output << "source_index = " << p.source_index << '\n';
    output << "source_weights = ";
    for (std::size_t s = 0; s < p.source_weights.size(); ++s) {
        if (s > 0) output << ',';
        output << p.source_weights[s];
    }
    output << '\n';
    if (!output)
        throw std::runtime_error("Cannot write configuration: " + path.string());
}

// D0 is the converted spatial diffusion coefficient in kpc^2/yr at 1 GV.
void write_run_parameters(const std::filesystem::path &path, const ModelParameters &p, const std::vector<cross_section::Isotope> &isotopes, double D0, const CFLNumbers &cfl) {
    std::ofstream parameters(path);
    parameters << std::setprecision(17)
               << "mode=" << p.mode << "\nsteps=" << p.steps << "\ndt_yr=" << p.dt_yr << "\ntime_yr=" << p.steps * p.dt_yr
               << "\nD0_cm2_per_s_at_1GV=" << p.D0_cm2_per_s
               << "\nD0_kpc2_per_yr_at_1GV=" << D0 << "\ndelta=" << p.delta << "\nvA_km_per_s=" << p.vA_km_per_s
               << "\nnH_cm-3=" << p.nH_cm3 << "\nhelium_enhancement=" << p.helium_enhancement << "\nsource_rate=" << p.source_rate
               << "\nspatial_grid=" << p.nR << 'x' << p.nz << ", R=[0," << p.R_max_kpc << "], z=[0," << p.z_max_kpc << "] kpc"
               << "\nrigidity_grid=" << p.nRg << " logarithmic bins, [" << p.Rg_min_GV << ',' << p.Rg_max_GV << "] GV"
               << "\nspatial_boundaries=reflecting lower, zero outer ghost density"
               << "\nrigidity_boundaries=adiabatic outflow, zero diffusive flux"
               << "\nnuclear_method=backward Euler, straight-ahead, constant cross sections"
               << "\nproton_reactions=not modeled\nradioactive_decay=not modeled\n";
    for (std::size_t s = 0; s < isotopes.size(); ++s) {
        parameters << "source_weight_" << cross_section::nuclide(isotopes[s]).name << '=' << p.source_weights[s] << '\n';
    }
    parameters << "cfl_definition=dt times maximum discrete diagonal loss rate\n"
               << "cfl_advection=" << cfl.advection << '\n'
               << "cfl_advection_limit=" << cfl.advection_limit << '\n'
               << "cfl_spatial_diffusion_R=" << cfl.spatial_diffusion_R << '\n'
               << "cfl_spatial_diffusion_z=" << cfl.spatial_diffusion_z << '\n'
               << "cfl_spatial_diffusion=" << cfl.spatial_diffusion << '\n'
               << "cfl_momentum_diffusion=" << cfl.momentum_diffusion << '\n'
               << "cfl_adiabatic=" << cfl.adiabatic << '\n'
               << "cfl_momentum_combined=" << cfl.momentum << '\n'
               << "nuclear_loss_dt=" << cfl.nuclear_loss << '\n'
               << "cn_nonnegative_rhs_sufficient_limit=2 per ADI axis or combined momentum operator\n"
               << "implicit_indicators=reported without timestep rejection\n"
               << "source_timestep_constraint=none for prescribed additive injection\n";
    if (!parameters)
        throw std::runtime_error("Cannot write run parameters: " + path.string());
}

void write_density_distributions(const std::filesystem::path &directory, const std::vector<cross_section::Isotope> &isotopes, const std::vector<Field3D> &density) {
    std::filesystem::create_directories(directory);
    for (std::size_t s = 0; s < isotopes.size(); ++s) {
        const std::string name(cross_section::nuclide(isotopes[s]).name);
        const std::filesystem::path path = directory / (name + ".bin");
        std::ofstream output(path, std::ios::binary | std::ios::trunc);
        if (!output)
            throw std::runtime_error("Cannot open density output: " + path.string());
        const Field3D &q = density[s];
        const auto row_bytes = static_cast<std::streamsize>(q.nz() * sizeof(double));
        for (int k = 0; k < q.nrigidity(); ++k) {
            for (int i = 0; i < q.nR(); ++i) {
                // Physical z cells are contiguous; skip the ghost cells at both ends.
                output.write(reinterpret_cast<const char *>(&q(k, i, 0)), row_bytes);
            }
        }
        output.close();
        if (!output)
            throw std::runtime_error("Failed writing density output: " + path.string());
    }
}

std::vector<cross_section::Isotope> select_isotopes(const ModelParameters &p) {
    using cross_section::Isotope;
    const std::string &mode = p.mode;
    std::vector<Isotope> isotopes;
    if (mode == "proton") {
        isotopes = {Isotope::Proton};
    } else if (mode == "network") {
        isotopes = {Isotope::O16, Isotope::N15, Isotope::N14, Isotope::C13, Isotope::C12,
                    Isotope::B11, Isotope::B10, Isotope::Be10, Isotope::Be9, Isotope::Be7};
    } else {
        std::istringstream names(mode);
        std::string name;
        while (std::getline(names, name, ',')) {
            int found = -1;
            for (std::size_t s = 0; s < cross_section::isotope_count; ++s) {
                if (cross_section::nuclides[s].name == name) {
                    found = static_cast<int>(s);
                    break;
                }
            }
            if (found < 0)
                throw std::runtime_error("Unknown isotope: " + name);
            const Isotope isotope = static_cast<Isotope>(found);
            if (std::find(isotopes.begin(), isotopes.end(), isotope) != isotopes.end())
                throw std::runtime_error("Duplicate isotope: " + name);
            isotopes.push_back(isotope);
        }
        if (isotopes.empty())
            throw std::runtime_error("Select at least one isotope");
    }

    if (p.source_weights.size() != isotopes.size())
        throw std::runtime_error("source_weights must match the selected isotope count");
    return isotopes;
}

CFLNumbers calculate_cfl_numbers(const ModelParameters &p, const SpatialGrid &grid, const RgGrid &Rg, const Field2D &vR, const Field2D &vZ, const Field2D &gas, const cr_diffusion::DiffusionModel &model, const std::vector<cross_section::Isotope> &isotopes) {
    CFLNumbers cfl;
    double max_outflow_rate = 0.0;
    double max_gas = 0.0;
    std::vector<double> rigidity_diffusion_rate(Rg.size(), 0.0);
    // Match momentum_solver.cpp: zero diffusive boundary flux, two-point
    // interior flux -Rg^2 D_RgRg d(q/Rg^2)/dRg, divided by receiving bin width.
    for (int k = 1; k < Rg.size(); ++k) {
        const double face = Rg.face(k);
        const double conductance = face * face * model.rigidity_coefficient(face) / Rg.center_distance(k - 1);
        const double left = conductance / (Rg.center(k - 1) * Rg.center(k - 1));
        const double right = conductance / (Rg.center(k) * Rg.center(k));
        rigidity_diffusion_rate[k - 1] += left / Rg.width(k - 1);
        rigidity_diffusion_rate[k] += right / Rg.width(k);
    }
    for (int i = 0; i < grid.nR(); ++i) {
        const double radial_volume = grid.R().center(i) * grid.R().width(i);
        double radial_rate = grid.R().face(i + 1) / (radial_volume * grid.radial_centroid_distance(i));
        if (i > 0)
            radial_rate += grid.R().face(i) / (radial_volume * grid.radial_centroid_distance(i - 1));
        for (int j = 0; j < grid.nz(); ++j) {
            const double volume = grid.volume(i, j);
            const double area_z = grid.vertical_face_area(i);
            // Advection lower faces are reflecting, irrespective of input velocity.
            double outflow = grid.radial_face_area(i + 1, j) * std::max(vR(i + 1, j), 0.0)
                           + area_z * std::max(vZ(i, j + 1), 0.0);
            if (i > 0)
                outflow += grid.radial_face_area(i, j) * std::max(-vR(i, j), 0.0);
            if (j > 0)
                outflow += area_z * std::max(-vZ(i, j), 0.0);
            max_outflow_rate = std::max(max_outflow_rate, outflow / volume);
            max_gas = std::max(max_gas, gas(i, j));

            // Match the diffusion stencil, including its zero outer ghost values.
            double vertical_rate = 1.0 / (grid.z().width(j) * grid.z().center_distance(j));
            if (j > 0)
                vertical_rate += 1.0 / (grid.z().width(j) * grid.z().center_distance(j - 1));
            const double divergence = (grid.radial_face_area(i + 1, j) * vR(i + 1, j)
                                     - grid.radial_face_area(i, j) * vR(i, j)
                                     + area_z * (vZ(i, j + 1) - vZ(i, j))) / volume;
            const double adiabatic_rate = -divergence / 3.0;
            for (int k = 0; k < Rg.size(); ++k) {
                const double D = model.spatial_coefficient(Rg.center(k));
                const double radial = p.dt_yr * D * radial_rate;
                const double vertical = p.dt_yr * D * vertical_rate;
                cfl.spatial_diffusion_R = std::max(cfl.spatial_diffusion_R, radial);
                cfl.spatial_diffusion_z = std::max(cfl.spatial_diffusion_z, vertical);
                cfl.spatial_diffusion = std::max(cfl.spatial_diffusion, radial + vertical);
                // Both rigidity boundaries allow outgoing adiabatic flux in main.
                const double cooling_heating = (std::max(adiabatic_rate * Rg.face(k + 1), 0.0)
                                              - std::min(adiabatic_rate * Rg.face(k), 0.0)) / Rg.width(k);
                const double diffusion = p.dt_yr * rigidity_diffusion_rate[k];
                const double adiabatic = p.dt_yr * cooling_heating;
                cfl.momentum_diffusion = std::max(cfl.momentum_diffusion, diffusion);
                cfl.adiabatic = std::max(cfl.adiabatic, adiabatic);
                cfl.momentum = std::max(cfl.momentum, diffusion + adiabatic);
            }
        }
    }
    // The triangular nuclear network has destruction rates on its diagonal.
    // Use the same beta, cross sections and year conversion as nuclear_solver.cpp.
    for (cross_section::Isotope isotope : isotopes) {
        const auto id = cross_section::index(isotope);
        const auto &nuclide = cross_section::nuclides[id];
        for (int k = 0; k < Rg.size(); ++k) {
            const double momentum_per_nucleon = Rg.center(k) * nuclide.charge / nuclide.mass_number;
            const double mass = physical_constants::nucleon_rest_energy_GeV;
            const double beta = momentum_per_nucleon / std::sqrt(momentum_per_nucleon * momentum_per_nucleon + mass * mass);
            const double loss = p.dt_yr * max_gas * beta * physical_constants::speed_of_light_cm_s
                              * physical_constants::seconds_per_year * p.helium_enhancement * cross_section::sigma[id][id];
            cfl.nuclear_loss = std::max(cfl.nuclear_loss, loss);
        }
    }
    cfl.advection = p.dt_yr * max_outflow_rate;
    if (!std::isfinite(cfl.advection) || cfl.advection > cfl.advection_limit) {
        std::ostringstream message;
        message << std::setprecision(17) << "Advection CFL = " << cfl.advection
                << " exceeds " << cfl.advection_limit
                << "; require dt_yr <= " << cfl.advection_limit / max_outflow_rate;
        throw std::runtime_error(message.str());
    }
    return cfl;
}
