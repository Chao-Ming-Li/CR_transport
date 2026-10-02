# Cosmic-ray transport in the Galactic circumgalactic medium

Implementation review: 2026-10-02. This document describes the current working-tree sources, including uncommitted changes. “Implemented” below refers to code present in the sources; it does not imply that the current executable builds or that all numerical properties have been verified.

## Physical motivation and intended model

The goal is to develop a C++ cosmic-ray (CR) transport code for the Milky Way disk and its extended circumgalactic medium (CGM). The intended domain extends from disk scales to hundreds of kiloparsecs. Compared with the conventional tens-of-kpc domains motivating this project, the intended improvement is the ability to resolve the disk while following transport into the CGM. Quantitative comparisons with GALPROP and DRAGON2 remain future work.

The schematic target equation for phase-space distribution of species i is

$$
\frac{\partial f_i}{\partial t}
= \nabla_r\cdot(D_{rr}\nabla_r f_i-\mathbf{v}_r f_i)
+ \nabla_p\cdot(D_{pp}\nabla_p f_i-\mathbf{v}_p f_i)
+ Q_i
+ \sum_{j\ne i}\frac{f_j}{\tau_{j\to i}}
- \frac{f_i}{\tau_i}.
$$

The intended physics includes spatial diffusion and advection, momentum diffusion/reacceleration, continuous energy losses, nuclear fragmentation, and radioactive decay. This is a target equation, not the equation currently solved in full.

Currently we have not considered the energy dependence, and are still implementing the spatial transport: diffusion, advection...
$$
\frac{\partial n}{\partial t}
= -\frac{1}{R}\frac{\partial(Rv_Rn)}{\partial R}
  -\frac{\partial(v_zn)}{\partial z}
  +\frac{1}{R}\frac{\partial}{\partial R}
    \left(RD\frac{\partial n}{\partial R}\right)
  +\frac{\partial}{\partial z}\left(D\frac{\partial n}{\partial z}\right).
$$

Here D is a constant isotropic scalar. Length and time units are kpc and yr; velocity and diffusion units are kpc/yr and kpc^2/yr. The normalization and physical units of n are not yet specified.

## Geometry, grids, and field storage

Implemented in `field.hpp` and `field.cpp`:

- Axisymmetric cylindrical R-z geometry, with no azimuthal dependence. The production setup covers positive R and z with boundary faces at R=0 and z=0; the latter assumes midplane symmetry.
- `AxisGrid::linear`: uniform widths with an exact requested endpoint.
- `AxisGrid::geometric`: a fixed cell count and successive-width ratio, with exact requested endpoints; ratio=1 reduces to a linear grid.
- `AxisGrid::geometric_capped`: widths grow from an initial width, then remain at a maximum width. The upper coordinate is a target: the last cell retains its full width, so the actual boundary can exceed the target by at most the maximum width, up to rounding. This design aims to identify the tiny structures on the disk as well as control the numerical diffusion close to the outer boundary.
- `AxisGrid::from_faces`: explicitly supplied strictly increasing physical faces.
- Two ghost-cell layers per axis; ghost geometry is reflected across the boundary faces.
- `Grid2D` provides annular volumes, radial/vertical face areas, and cached radial volume centroids and centroid distances. The arithmetic radial center used in volume factors differs from the volume centroid used in reconstruction.
- `Field2D` stores density at `Centroid`, vR at `RadialFace`, and vZ at `VerticalFace`. Storage is contiguous, includes ghosts, and has z as the fastest-varying index. Location/dimension checks are available through `matches`.

For cell faces R_- and R_+, the volume is

$$V_{ij}=\pi(R_+^2-R_-^2)\Delta z_j.$$

The conservative radial denominator `R.center(i) * R.width(i)` equals (R_+^2-R_-^2)/2. Thus the implemented flux differences account for cylindrical geometry even on nonuniform grids. A mass diagnostic must use sum(n_ij V_ij) over physical cells only. The grid volume includes the full 2pi azimuth but only the modeled positive-z half; a symmetric full-domain integral requires the corresponding factor of two.

## Boundary conditions

`apply_boundary_conditions` in `CR_transport_solver.cpp` sets reflective boundary at R = 0 and z = 0, and absorbing boundary at $\rm R_{max}$ and $z_{max}$. The assumption is that the spatial distribution is central symmetric and there is no net flux accross inner boundaries. While when $\rm R_{max}$ and $z_{max}$ get large enough, they do not affect the results any more.


## Spatial advection

The active implementation is `cr_advection::Solver` in `CR_advection.hpp` and `CR_advection.cpp`.

- Constant reconstruction gives first-order upwind spatial fluxes.
- PLM reconstructs a linear profile from neighboring densities using actual radial volume-centroid or vertical-center distances.
- Available slope limiters are Minmod, Van Leer, and Monotonized Central (MC). Reconstructed face states are clamped between the donor and adjacent cell averages; this is a face-state bound, not clipping of evolved cell averages.
- Face states are selected by velocity sign. Each shared face flux is computed once per spatial-operator evaluation.
- `Unsplit` updates both spatial directions together, using Euler or SSPRK2. SSPRK2 uses two Euler stages and their convex combination.
- `RadialThenVertical` performs sequential Euler directional updates in the present implementation; the header describes this combination as Euler-only.
- Velocities are held fixed over an advance call. Reusable stage, RHS, and flux fields belong to the solver object.

The production wrapper selects PLM + MC + unsplit SSPRK2; the default `Options` limiter is Van Leer. `PPM` exists as an enum value, but no separate PPM reconstruction is implemented: the current face-value routine treats every non-Constant choice as linear reconstruction. Do not describe PPM as an available implemented method.

There is no automatic CFL calculation, rejection, or substepping in the current `advance` implementation. The caller must choose a stable dt on the actual grid. Face-state limiting alone does not guarantee nonnegative updated cell averages for an arbitrary timestep. Claimed global spatial order, positivity limits, and convergence rates require validation for this version.

## Spatial diffusion

`cr_diffusion::Solver` is implemented in the same active header/source pair. It solves constant-coefficient isotropic diffusion with the Peaceman–Rachford Crank–Nicolson ADI scheme:

$$
(I-\tfrac{\Delta t}{2}L_R)n^*
=(I+\tfrac{\Delta t}{2}L_z)n^k,
$$

$$
(I-\tfrac{\Delta t}{2}L_z)n^{k+1}
=(I+\tfrac{\Delta t}{2}L_R)n^*.
$$

Directional coefficients use shared face conductances divided by cell volume, with radial centroid distances and vertical center distances. The tridiagonal systems use cached Thomas factors. Grid, D, timestep, and boundary options are fixed at construction; construct a new solver when these change. Working fields and line storage are reused.

D=0 or dt=0 disables the diffusion update. The method is designed to be second order in time for the fixed linear operator, but Crank–Nicolson does not preserve positivity at arbitrary large timesteps. Signed results are retained; final physical-cell results are checked for finiteness. Spatially varying or anisotropic diffusion is not implemented.

## Current driver and initial conditions

`main.cpp`, `initialization.cpp`, and `CR_transport_solver.cpp` define a transport demonstration:

| Quantity | Current setting |
| --- | --- |
| R grid | Target [0,200] kpc; first width 0.05 kpc; maximum width 1 kpc; growth ratio 1.03 |
| z grid | Target [0,200] kpc; first width 0.01 kpc; maximum width 1 kpc; growth ratio 1.03 |
| Initial CR density | exp[-(R_centroid-30)^2/9-(z_center-30)^2/9] multiplied by DT |
| Wind | vR=vZ=3e-7 kpc/yr away from the lower reflecting faces |
| Gas | Uniform `ndis_H=0.001`; unused by the current transport update |
| Timestep / step count | DT=10000 yr; NT=10000; total requested evolution 1e8 yr |
| Diffusion | Default D=0 in the wrapper; `main.cpp` does not override it |

`initialize_CR_source` assigns an initial Gaussian once. Despite its name, it does not supply continuous injection Q during evolution. Its DT-dependent amplitude is a demonstration convention, not a documented physical source normalization. `ndis_C` is a single scalar field, not a completed carbon-isotope transport model.

The driver performs advection(dt), then diffusion(dt) each iteration. This is first-order operator splitting when both processes are active, even though the individual SSPRK2 and ADI integrators are second order in time. The unsplit advection option refers only to combining R and z within advection. There is no Strang composition of advection and diffusion in the current driver.

After the first update and every ten iterations thereafter, output appends the entire density storage, including ghost cells, as native binary doubles. With D=0 the filename is `ndis_C_transport_dR50_dz10_r1.03_dt1e4_MC_fix.bin`; otherwise it is `ndis_C_advection_diffusion.bin`. There is no embedded grid/time metadata or initial-state record, and an existing file is appended to on a repeated run.

## Nuclear data and unimplemented physics

`cross_section.hpp` contains isotope metadata and a fixed 10-by-10 cross-section table for O16, N15, N14, C13, C12, B11, B10, Be10, Be9, and Be7. The network data start from O16. The table convention is `sigma[parent][daughter]`: diagonal entries represent total inelastic destruction and off-diagonal entries partial fragmentation. Values are in cm^2; comments attribute them to DRAGON2 values converted from mb.

The header also defines unit/physical constants, including a default helium enhancement factor. It does not currently implement reaction-rate functions, energy-dependent cross sections, decay lifetimes, or a reaction solver. The table is not connected to the production transport loop. Its energy applicability and precise provenance remain to be documented.

Pending physics:

- Momentum/energy grid and variable convention, momentum diffusion, and continuous losses (including a defined treatment of wind-related adiabatic changes).
- Coupled fragmentation production/destruction, radioactive decay, and the associated rate/unit conversions.
- Hadronic/pion-production losses with an explicit distinction between continuous losses and catastrophic removal.
- Physically normalized injection, disk/CGM gas profiles, and species-dependent transport where needed.

## Build and validation status

The Makefile's production source list is `main.cpp`, `CR_advection.cpp`, `CR_transport_solver.cpp`, `field.cpp`, and `initialization.cpp`, using C++20 and Homebrew LLVM/libomp. Legacy files such as `CR_advection_old.cpp` and `CR_advection_RK.cpp` are not part of that target. OpenMP headers/build flags are present, but the active transport sweeps currently contain no OpenMP parallel directives.

A fresh compile/link of those five production sources on 2026-10-02 failed: `cr_advection::Solver::Solver(const Grid2D&, const Options&)` is declared and called but has no definition in the current active sources. Consequently this snapshot is not a buildable production solver. The unused `require` and `finite_nonnegative` helpers also show that header claims of input validation are not currently realized: advection `advance` lacks dimension/timestep validation, and the diffusion constructor lacks its advertised finite/nonnegative D and dt checks. These are implementation gaps, not completed guarantees.

Existing validation sources include:

- `test_grid.cpp`: grid construction, ghost geometry, cylindrical volumes/areas, cached centroids, and invalid inputs.
- `test_diffusion.cpp`: nonuniform-grid conservation, constant preservation, discrete-mode amplification, temporal refinement, absorbing boundaries, and validation expectations.
- `test_solvers.cpp`: reused-versus-fresh solver results and validation/no-op expectations.
- `test_upwind.cpp`: conservative upwind transport, both velocity signs, positivity examples, and CFL/validation expectations.
- Advection benchmarks/comparisons: additional historical tests and measurements, whose interfaces must be reconciled with the active solver.

These are test intentions, not a report that the current snapshot passes. Several tests expect validation/CFL behavior missing from the current implementation. Some Makefile test targets reference absent `grid.cpp` or `CR_advection_unified.cpp`; `test_advection_unified.cpp` also includes an absent `CR_advection_unified.hpp`. Existing executables or object files are not evidence for the current sources.

`NUMERICS.md`, `ADVECTION_UNIFIED.md`, and the comparison documents describe earlier variants. For example, `NUMERICS.md` describes quadratic reconstruction and automatic CFL substeps that are absent from the current active source. Historical benchmark results must not be attributed to this version without rerunning compatible tests.

## Next development priorities

1. Restore the advection constructor and reconcile actual input validation, supported options, and CFL handling with the public contract.
2. Repair test/build dependencies and run fresh tests of the active implementation; establish conservation, boundary flux balance, positivity under stated timestep limits, and spatial/time convergence.
3. Specify the density/source normalization and output metadata; distinguish numerical Gaussian demonstrations from physical disk/CGM simulations.
4. If second-order combined transport is required, implement and verify a suitable symmetric operator composition. Do not infer overall order from the individual integrators.
5. Define momentum/species variables and reaction data conventions, then integrate the missing physical processes with isolated tests before full coupled runs.
