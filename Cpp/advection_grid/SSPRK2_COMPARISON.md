# Van Leer: forward Euler versus SSPRK2

Reproduce with `make compare-ssprk2`. Source: `compare_ssprk2.cpp`.
Full metrics: `ssprk2_metrics.csv`; profiles: `ssprk2_profiles.csv`.
No production solver source was modified.

## Experiment

Both methods use the current van Leer reconstruction. Forward Euler calls
`advection_TVD`; SSPRK2 calls `advance_TVD_SSPRK2`. Set vR=0, vZ=0.1,
t=500, and use identical radial rows to isolate vertical advection. The
sequential R/Z structure of advection_TVD has no splitting effect when vR=0.
The harness refreshes input ghost cells before each timestep, since the
SSPRK2 wrapper does not do so before its first stage.

Initial Gaussian: exp(-(z-30)^2/9), exact finite-volume cell averages,
unit amplitude. The reference is the same Gaussian translated to z=80.
The [0,200] domain avoids boundary mass loss present on the shorter production
domain. Also test a unit top hat initially on [25,35]. Uniform dz=1 and 0.5;
dt=0.1,1,2,5 where CFL<=0.5. These are vertical tests, not validation of
combined R/z stability, radial geometry, or stretched grids.

Variance is measured using the piecewise constant cell distribution,
including dz^2/12. The exact translation preserves this variance.
D_effective = variance_change / (2*t); negative means net narrowing.
This is a profile-dependent diagnostic, not a constant material diffusivity.
Relative L1 is integrated absolute error divided by exact mass.

## Gaussian results

| dz | dt | Integrator | Variance change | D_effective | Peak | Relative L1 |
|---:|---:|---|---:|---:|---:|---:|
| 1 | 0.1 | vanLeer_FE | 4.583824 | 0.0045838 | 0.625633 | 38.40% |
| 1 | 0.1 | vanLeer_SSPRK2 | 4.912081 | 0.0049121 | 0.616969 | 39.68% |
| 1 | 1 | vanLeer_FE | 1.990269 | 0.0019903 | 0.711075 | 26.07% |
| 1 | 1 | vanLeer_SSPRK2 | 4.936263 | 0.0049363 | 0.616659 | 39.73% |
| 1 | 2 | vanLeer_FE | -0.020136 | -0.0000201 | 0.813668 | 16.06% |
| 1 | 2 | vanLeer_SSPRK2 | 5.014706 | 0.0050147 | 0.615451 | 39.92% |
| 1 | 5 | vanLeer_FE | -1.833509 | -0.0018335 | 0.961241 | 27.91% |
| 1 | 5 | vanLeer_SSPRK2 | 5.603825 | 0.0056038 | 0.609888 | 42.17% |
| 0.5 | 0.1 | vanLeer_FE | 0.532097 | 0.0005321 | 0.838922 | 13.33% |
| 0.5 | 0.1 | vanLeer_SSPRK2 | 0.924145 | 0.0009241 | 0.818269 | 14.99% |
| 0.5 | 1 | vanLeer_FE | -1.672670 | -0.0016727 | 0.968099 | 25.03% |
| 0.5 | 1 | vanLeer_SSPRK2 | 0.953819 | 0.0009538 | 0.817352 | 15.32% |
| 0.5 | 2 | vanLeer_FE | -2.025648 | -0.0020256 | 0.988597 | 31.34% |
| 0.5 | 2 | vanLeer_SSPRK2 | 1.054821 | 0.0010548 | 0.815069 | 16.52% |

## Interpretation

At dz=1, dt=5, SSPRK2 broadens the Gaussian (variance change +5.6038),
while forward Euler narrows it (-1.8335). SSPRK2 L1 error is 42.17%,
versus 27.91% for forward Euler. The exact cell-average peak is 0.964166;
SSPRK2 gives 0.609888 and forward Euler gives 0.961241. A preserved peak
alone does not establish that the shape is accurate.

SSPRK2 removes the leading forward-Euler antidiffusion term -v^2*dt/2.
The van Leer spatial reconstruction still has substantial limiter diffusion
on this coarse Gaussian. Forward Euler partly offsets that diffusion;
removing its time error can therefore increase total L1 error and broadening.
The coarse-grid measurements do not support a claim that SSPRK2 always
improves the numerical solution.

At dz=0.5, dt=1, SSPRK2 instead reduces L1 error from 25.03% to 15.32%,
and replaces narrowing (-1.6727) with broadening (+0.9538). At dt=2,
L1 falls from 31.34% to 16.52%. The ranking depends on resolution and timestep.

At fixed dz=1, SSPRK2 results approach a nonzero spatial error as dt decreases:
L1 is 42.17%, 39.92%, 39.73%, and 39.68% for dt=5,2,1,0.1.
Reducing dt alone cannot remove limiter diffusion. At dt=0.1, refining dz
from 1 to 0.5 reduces SSPRK2 L1 from 39.68% to 14.99% and variance growth
from 4.9121 to 0.9241. This is not an order-of-convergence measurement.

The earlier short-transport benchmark (sigma=0.4, dz=0.0625, shift=2)
had a different resolution and propagation distance; its improvement with
SSPRK2 should not be extrapolated to this production-scale Gaussian.

## Top hat and validation

At dz=1, dt=5, top-hat L1 error is 5.57% for forward Euler versus 37.71%
for SSPRK2; variance growth is 0.2791 versus 5.9366. Forward Euler maintains
a much sharper front in this test.

All 28 runs completed. Final mass errors are below 1e-13, minima are
nonnegative to tolerance 1e-13, and final total variation does not exceed
the initial value (tolerance 1e-12). These are final-state checks, not a proof
of stability for arbitrary data or timesteps.

Use SSPRK2 for consistent second-order time integration; use finer spatial
resolution to address the remaining diffusion. Select based on full-profile
error and broadening together, rather than peak height alone.
