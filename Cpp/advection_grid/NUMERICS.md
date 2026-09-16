# Advection update

The solver evolves nonnegative finite-volume cell averages of density in
axisymmetric cylindrical geometry. Radial fluxes include the face radius and
are divided by `R_center * dR`; this is equivalent to face area / annular volume.
The velocity field is held fixed over each requested interval.

## What changed

- Replaced forward-Euler directional sweeps with an unsplit SSPRK(2,2) update.
  The old Strang sequence used first-order subsolvers and was not second order
  in time. The directional public functions now use the same SSPRK2 engine.
- Replaced van Leer ratios containing a fixed `1e-12` with a limited quadratic
  fit to three cell averages. Reconstruction accounts for nonuniform widths,
  cylindrical volume centroids, and second moments.
- A five-cell curvature check permits the quadratic in smooth regions,
  including resolved extrema: neighboring second differences must have the
  same sign and magnitudes within a factor of two. Elsewhere, face increments
  are bounded by both neighboring differences; unresolved extrema are flattened.
  This is a local smoothness heuristic, not an implementation of a named PPM
  algorithm. Boundary cells retain constant reconstruction.
- Each face state is bounded to `[0, 2 * donor_density]`. Automatic substeps
  keep the summed outgoing area-weighted Courant number at most 0.45, below
  the positivity limit of 0.5. Positivity follows at each Euler stage from
  this bound; SSPRK2 is a convex combination of these stages. Cell averages
  are never clipped, which preserves the conservative flux balance.
- Lower faces have zero flux without changing the supplied velocities. Upper
  faces use interior outflow and zero-density inflow. Shared face fluxes are
  computed once per stage. All public advection functions validate inputs.

The historical `advection_TVD*` names remain for source compatibility, but the
smoothness exception does **not** guarantee strict TVD or a maximum principle
for arbitrary data. It trades that guarantee for better smooth-peak accuracy.
Tests cover top-hat bounds and total variation, positivity, conservation,
small density scaling, both velocity signs, and uniform/stretched grids.
Boundary accuracy remains first order locally. Two flux stages and the wider
reconstruction cost more per step than the old Euler sweeps.

The motivation for preserving smooth extrema is discussed by
[Colella and Sekora](https://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/ColellaSekora.pdf).
Their published limiter differs from the curvature check implemented here.
[Clawpack's solver documentation](https://www.clawpack.org/v5.10.x/pyclaw/solvers.html)
also describes combining higher-order spatial reconstruction with Runge–Kutta
integration.

## Measured comparison

Use the existing Gaussian setup: 100 vertical cells on `[0,100]`,
`n = exp(-(z-30)^2/9)`, `vZ=0.1`, `dt=1`, 500 steps. The reference is the
same sampled Gaussian translated to `z=80`. A two-cell radial grid with
identical rows isolates vertical transport; it has the same vertical solution
as the 100-row production grid. The baseline was compiled from the original
`CR_advection.cpp` before edits and called its directional Z update.

| Metric | Original | Updated |
| --- | ---: | ---: |
| L1 error, sum times dz | 1.41792 | 1.12224 |
| Peak density | 0.713125 | 0.782662 |
| Vertical variance | 6.53811 | 6.25314 |

The exact sampled peak is 0.972604 and variance is 4.5. L1 error decreases
by about 21%; excess variance decreases by about 14%. Both velocity signs
produce the same updated metrics. These gains describe this benchmark,
not every initial profile or grid. Finer resolution still helps substantially.

## Reproduce validation

```sh
make all test-advection test-grid test-field
```

`test-advection` prints the long-transport metrics above and checks regression
thresholds, mass conservation, positive density, boundary flux balance,
directional/unsplit consistency, CFL substeps, and density-scale independence.
Time-step halving yields error ratios 4.09 and 4.05, consistent with second-order
time accuracy. Gaussian grid-refinement ratios are approximately 5.4–7.8 for
the tested vertical cases and 7.6–7.7 for radial transport. These are measured
ratios, not a claim of global third-order accuracy with limiting and boundaries.

The field tests also had stale initializer names and an outdated factor of
1000 in source expectations; these now match the existing initialization code.
The production binary still appends its final state to `ndis_C_t_ssad.bin`.
