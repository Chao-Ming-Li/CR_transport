# Numerical diffusion comparison of the current code

The production schemes were called directly, without modifying CR_advection.cpp.
Reproduce with `make compare-diffusion`. Outputs are `diffusion_metrics.csv`
and `diffusion_profiles.csv` (initial, exact translated, and numerical profiles).

## Setup and interpretation

Pure vertical transport, vZ=0.1, vR=0, final time 500, translation 50.
Gaussian initial condition exp(-(z-30)^2/9), integrated exactly over each cell;
also a unit top hat on [25,35]. Two identical radial rows isolate vertical
transport. Uniform dz=1 and 0.5, dt=0.1, 1, 2, 5 where CFL <=0.5.
The domain is [0,200], extended beyond the production [0,100] to avoid
outflow bias in variance and mass. Initial Gaussian amplitude is normalized
to one rather than multiplied by DT as in initialize_CR_source.

Variance uses the piecewise constant cell distribution, including dz^2/12.
D_effective = (variance_final - variance_initial)/(2*500). Negative values
mean net narrowing, not physical diffusion. This is a profile-dependent
diagnostic, not a universal diffusion coefficient for nonlinear limiters.
Relative L1 is integrated absolute error divided by exact integrated density.
The translation is an integer number of cells on both grids, so the exact
cell-average profile has the same variance and peak as the initial profile.

## Gaussian, dz=1

Exact cell-average peak: 0.964166. Continuous Gaussian variance: 4.5;
the piecewise constant cell-average reference variance is 4.666667.

| dt | Scheme | Variance change | D_effective | Peak | Relative L1 |
|---:|---|---:|---:|---:|---:|
| 0.1 | Upwind | 49.5000 | 0.049500 | 0.288662 | 105.58% |
| 0.1 | TVD / van Leer | 4.5838 | 0.004584 | 0.625633 | 38.40% |
| 0.1 | PLM / MC | 2.5816 | 0.002582 | 0.696256 | 29.06% |
| 1 | Upwind | 45.0000 | 0.045000 | 0.301251 | 102.58% |
| 1 | TVD / van Leer | 1.9903 | 0.001990 | 0.711075 | 26.07% |
| 1 | PLM / MC | 0.0755 | 0.000076 | 0.815046 | 14.81% |
| 2 | Upwind | 40.0000 | 0.040000 | 0.317373 | 99.19% |
| 2 | TVD / van Leer | -0.0201 | -0.000020 | 0.813668 | 16.06% |
| 2 | PLM / MC | -1.1432 | -0.001143 | 0.911045 | 18.54% |
| 5 | Upwind | 25.0000 | 0.025000 | 0.387689 | 84.80% |
| 5 | TVD / van Leer | -1.8335 | -0.001834 | 0.961241 | 27.91% |
| 5 | PLM / MC | -1.9517 | -0.001952 | 0.964166 | 31.19% |

Upwind broadens strongly. Its measured diffusion matches
D_num = |v| dz (1-CFL)/2, with CFL=|v|dt/dz. Smaller timesteps at fixed
spatial resolution increase this diffusion; they do not remove spatial error.

PLM is less diffusive than van Leer here, but not always more accurate:
at dt=2 and 5, van Leer has lower Gaussian L1 error. At dt=5, preserving
the peak conceals substantial narrowing/distortion. Near-zero variance change
also does not imply an exact profile (PLM at dt=1 still has 14.81% L1 error).

## Top hat and refinement

At dz=1, dt=5, top-hat relative L1 errors are 77.95% (upwind),
5.57% (van Leer), and 5.00% (PLM). Variance increases are 25,
0.2791, and 0.25 respectively. No significant overshoot or increase in total
variation occurs in the tested final profiles. PLM minima are negative only
at floating-point roundoff scale (roughly 1e-16); these are not clipped.

At dt=0.1, refining dz from 1 to 0.5 reduces Gaussian L1 error from
105.58% to 84.05% for upwind, 38.40% to 13.33% for van Leer,
and 29.06% to 9.35% for PLM. This fixed-dt comparison is not an order-of-
convergence study: CFL changes with resolution. All 42 runs finish, and
relative mass changes are below 1e-13. Upwind variance increments agree
with the analytical discrete-moment prediction to floating-point accuracy.

## What the implementation means

- advection_upwind uses piecewise constant donor-cell fluxes and forward Euler.
  The production driver currently calls advection_upwind_1 instead. For this
  uniform-grid, positive constant vertical velocity test the two updates are
  algebraically equivalent. They are not generally equivalent for variable
  velocities, radial flow, negative velocities, or nonuniform grids.
- advection_TVD uses van Leer limited linear reconstruction and sequential
  R/Z forward-Euler sweeps. advection_PLM uses MC limited linear reconstruction
  and the same time/splitting structure. Both therefore use piecewise linear
  reconstruction; the comparison is largely van Leer versus MC limiting.
- Neither includes a Hancock half-step predictor or SSPRK2 time integration.
  For smooth regions with second-order spatial reconstruction, forward Euler
  contributes a leading negative diffusion term -v^2 dt/2. Limiter diffusion
  can offset this, explaining the timestep-dependent broadening/narrowing.
  A TVD limiter name alone does not establish TVD for every timestep and grid.
- The fixed +1e-12 in van Leer ratios breaks exact density-scale invariance.
  These results use unit-amplitude data and cannot establish behavior at all
  physical density scales.
- NUMERICS.md describes SSPRK2 and quadratic reconstruction that are absent
  from the current CR_advection.cpp; its old comparison does not apply here.

For reliable second-order time accuracy, the next change would be SSPRK2
with the spatial reconstruction, or a suitable MUSCL-Hancock predictor,
followed by CFL, positivity, and convergence validation. No solver change
was made as part of this comparison. These results isolate vertical numerical
transport and do not validate radial geometry or stretched-grid behavior.
