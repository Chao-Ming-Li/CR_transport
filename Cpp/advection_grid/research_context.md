# Cosmic-ray transport: research and code reference

Reviewed against working-tree sources on **2026-10-05**. This describes the source code, not a validated executable: the production build currently has a missing advection constructor.

## Scope and physical model

Goal: model cosmic-ray transport from the Milky Way disk into the circumgalactic medium, extending to hundreds of kpc while retaining fine disk resolution. Comparisons with GALPROP/DRAGON2 are future work.

Currently, one scalar density field evolves by spatial advection and constant isotropic diffusion in axisymmetric cylindrical coordinates:

$$
\frac{\partial n}{\partial t}
= -\frac{1}{R}\frac{\partial(Rv_R n)}{\partial R}
  -\frac{\partial(v_z n)}{\partial z}
  +\frac{1}{R}\frac{\partial}{\partial R}\left(RD\frac{\partial n}{\partial R}\right)
  +D\frac{\partial^2 n}{\partial z^2}.
$$

Units: length **kpc**, time **yr**, velocity **kpc/yr**, diffusion coefficient **kpc²/yr**. Density/source normalization remains undefined. There is no energy or species dimension, continuous injection, reaction network, or loss term in the current transport loop.

## Code map

| Files | Responsibility |
| --- | --- |
| `field.hpp`, `field.cpp` | Grid geometry, ghost coordinates, field storage |
| `CR_advection.hpp`, `CR_advection.cpp` | Advection and diffusion solver classes |
| `CR_transport_solver.cpp` | Boundary filling, timestep checks, evolution loop, binary output |
| `main.cpp`, `initialization.cpp`, `initialization.hpp` | Demonstration grid, initial density, wind, gas, DT/NT |
| `cross_section.hpp` | Isotope metadata and fixed cross-section table; unused by transport |
| `Makefile`, `test_*.cpp` | Build targets and numerical tests; some dependencies are stale |

## Grid and boundaries

- Domain: positive R and z; z=0 represents midplane symmetry.
- Grid choices: uniform (`linear`), geometric with fixed cell count (`geometric`), capped geometric widths (`geometric_capped`), or explicit faces (`from_faces`).
- Capped grids retain the last full cell: the actual upper boundary may exceed the requested target by up to the maximum width, within rounding.
- Density is a cell average stored at `Centroid`; velocities are stored at `RadialFace` and `VerticalFace`. Two ghost layers surround each axis; z is the fastest storage index.
- Radial reconstruction uses volume centroids; flux denominators use the arithmetic radial center times cell width. These are distinct coordinates.
- Lower boundaries: reflected density ghosts and explicitly zero flux. Outer R and upper z boundaries: zero-density ghosts (absorbing). Boundary influence must be checked by enlarging the domain.

Particle-number diagnostics must sum **physical cells only**:

$$
N=\sum_{ij}n_{ij}V_{ij},\qquad
V_{ij}=\pi(R_{i+1}^2-R_i^2)\Delta z_j.
$$

Volumes include the full azimuth but only positive z; a symmetric full-domain integral is twice this value.

## Numerical methods

**Advection.** Production options are **PLM + MC limiter + unsplit SSPRK2**. PLM uses actual centroid/center distances and clamps reconstructed face states between adjacent cell averages. Velocity sign selects the donor cell; each shared face flux is computed once per stage. Evolved densities are not clipped.

Other options: constant reconstruction (first-order upwind), Minmod/Van Leer limiters, Euler integration, and sequential `RadialThenVertical` Euler sweeps. The default limiter in `Options` is Van Leer, overridden to MC by the wrapper. The advection constructor copies vR/vZ and stores a fixed timestep; `advance(density)` reuses them. Construct a new solver to change these inputs. There is no automatic timestep selection or substepping.

**Diffusion.** Peaceman–Rachford Crank–Nicolson ADI alternates implicit and explicit directions:

$$
(I-\tfrac{\Delta t}{2}L_R)n^*=(I+\tfrac{\Delta t}{2}L_z)n^k,
\qquad
(I-\tfrac{\Delta t}{2}L_z)n^{k+1}=(I+\tfrac{\Delta t}{2}L_R)n^*.
$$

Coefficients use face conductances divided by cell volume, with radial centroid distances and vertical center distances. Thomas factors and working storage are reused. Grid, D, and timestep are fixed at construction; lower boundaries reflect and outer boundaries absorb. D=0 disables diffusion. Variable/anisotropic D is not implemented.

**Combined update.** Each iteration applies advection(Δt), then diffusion(Δt): **first-order operator splitting** when both are active. Individual second-order time integrators do not make this composition second order. “Unsplit” refers to R/z advection only.

## Timestep checks

`solve_advection_equation` checks once before evolution because grid, velocities, D, and DT remain fixed.

**Advection:** reject the timestep if

$$
C_{\rm adv}=\Delta t\max_{ij}
\left[\frac{\sum_{\text{outgoing faces}}A_f|v_f|}{V_{ij}}\right]>0.5.
$$

Both directions contribute; reflecting faces are excluded. This is a velocity-based outgoing fraction, not the actual particle fraction, which also depends on reconstructed face density. The code intends 0.5 as a conservative PLM positivity bound. However, face states are clamped to adjacent averages, not explicitly to twice the donor density; a general positivity guarantee on arbitrary stretched grids still needs verification.

**Diffusion:** warn, but continue, if

$$
\frac{\Delta t}{2}\max\left[\max_i(a_{R,i}+b_{R,i}),\max_j(a_{z,j}+b_{z,j})\right]>1.
$$

Here a and b are the nonnegative lower/upper neighbor diffusion coefficients, including absorbing-boundary contributions. Within the bound, each explicit half-step has nonnegative weights and the implicit solves preserve nonnegativity, provided input/boundary densities are nonnegative and D≥0. Exceeding it can cause negative densities; it is not an explicit diffusion stability limit.

For uniform Cartesian cells, this ADI positivity bound reduces to **Δt ≤ min(ΔR², Δz²)/D**. The code uses cylindrical, nonuniform-grid coefficients instead.

## Current demonstration

| Setting | Value |
| --- | --- |
| R grid | Target 0–200 kpc; first width 0.05 kpc; maximum 1 kpc; ratio 1.03 |
| z grid | Target 0–200 kpc; first width 0.01 kpc; maximum 1 kpc; ratio 1.03 |
| Initial density | DT × exp[−(R_centroid−30)²/9 − (z_center−30)²/9] |
| Wind | vR=vZ=3×10⁻⁷ kpc/yr; lower-face velocities zero |
| Gas | Uniform `ndis_H=0.001`; unused in evolution |
| Time | DT=10⁴ yr; NT=10⁴; requested duration 10⁸ yr |
| Diffusion | D=0 by default; `main.cpp` uses this default |

`initialize_CR_source` initializes the Gaussian once; it is not continuous injection. Its DT-dependent amplitude is a demonstration convention. `ndis_C` is not yet a coupled carbon-isotope model.

Output follows updates 1, 11, 21, … and appends **all storage, including ghosts**, as native binary doubles:

- D=0: `ndis_C_transport_dR50_dz10_r1.03_dt1e4_MC_fix.bin`
- D≠0: `ndis_C_advection_diffusion.bin`

There is no initial-state record or grid/time metadata. Repeated runs append to existing files.

## Nuclear data and remaining physics

`cross_section.hpp` lists O16, N15, N14, C13, C12, B11, B10, Be10, Be9, and Be7. Its fixed 10×10 table uses `sigma[parent][daughter]`: diagonal entries are total inelastic destruction; off-diagonal entries are fragmentation production. Values are in cm²; comments attribute them to DRAGON2 values converted from mb. Energy applicability and precise provenance remain undocumented.

Constants include a helium enhancement factor of 1.3. Despite the header's introductory rate comment, no reaction-rate functions, decay lifetimes, or reaction solver are implemented.

Remaining physics: normalized injection and disk/CGM gas profiles; momentum/energy variables; reacceleration and continuous/adiabatic losses; fragmentation and decay; hadronic losses with a defined continuous-versus-catastrophic treatment.

## Build status and next steps

Production uses C++20 and Homebrew LLVM/libomp, compiling the five source files listed in `Makefile`. Legacy advection variants are excluded. Active transport sweeps have no OpenMP parallel directives.

Known gaps in the current sources:

1. **Build blocker:** `cr_advection::Solver` declares a constructor but has no definition. The production link failed on 2026-10-05.
2. **Input validation:** advection `advance` lacks the advertised dimension/timestep checks; diffusion construction lacks finite/nonnegative D and dt checks. The wrapper's current CFL scan also lacks dimension/finite-input checks before accessing velocities. `require` and `finite_nonnegative` helpers are unused.
3. **Test dependencies:** some targets reference absent `grid.cpp`, `CR_advection_unified.cpp`, or `CR_advection_unified.hpp`. Existing binaries do not establish that current sources pass.
4. **Historical notes:** `NUMERICS.md`, `ADVECTION_UNIFIED.md`, and comparison documents describe earlier variants. Quadratic reconstruction and automatic CFL substeps described there are absent from active code.

Priority: restore the constructor and input checks; repair test targets; run fresh conservation, boundary-flux, positivity, and convergence tests. Existing grid/diffusion/reuse/upwind tests provide starting points, not a current pass report. Then define normalization/output metadata, decide whether second-order combined transport is needed, and add the missing physics incrementally.
