# Unified advection

`CR_advection_unified.cpp` implements conservative cylindrical transport through
one face reconstruction and flux-divergence implementation. Its public API is in
`CR_advection_unified.hpp`. Existing solvers and the active simulation are unchanged;
the new implementation is included in the default build.

To use it, include the header and replace the advection call in the timestep loop:

```cpp
cr_advection::Options options;
options.reconstruction = cr_advection::Reconstruction::PLM;
options.limiter = cr_advection::Limiter::MC;
options.integrator = cr_advection::TimeIntegrator::SSPRK2;
options.splitting = cr_advection::Splitting::Unsplit;
cr_advection::advance(ndis_C, vR, vZ, DT, grid, options);
```

| Legacy routine | Reconstruction | Limiter | Integrator | Splitting |
|---|---|---|---|---|
| `advection_upwind` | Constant | ignored | Euler | Unsplit |
| `advection_TVD` | PLM | VanLeer | Euler | RadialThenVertical |
| `advection_PLM` | PLM | MC | Euler | RadialThenVertical |
| `advance_TVD_SSPRK2` | PLM | VanLeer | SSPRK2 | Unsplit |

Minmod is also supported. PPM/WENO are not implemented. The nonconservative,
vertical-only `advection_upwind_1` solves a different equation and is intentionally
not an option in this conservative solver.

R reconstruction uses cylindrical volume centroids; z uses cell centers. MC's
centered slope uses the actual neighbor distances on stretched grids. Van Leer
uses a sign-aware harmonic mean without the old absolute `1e-12` denominator
offset, so agreement with legacy van Leer is approximate, not bit-for-bit.

The lower faces always have zero flux, even if a supplied lower-face velocity is
nonzero. Outer ghost densities are zero. Input and stage ghosts are refreshed.
The caller must choose a stable timestep; positivity checks are not a CFL proof,
particularly for stretched grids. Negative or nonfinite stage densities throw
without committing any changes to the input. Density is never clipped.
SSPRK2 checks both full Euler stages before averaging. Velocities must remain
fixed during each call; time-dependent stage velocities require a different API.

`RadialThenVertical` is supported only with Euler to preserve the legacy split
method explicitly; it is not Strang splitting. Default options select unsplit
PLM/VanLeer/SSPRK2. Work arrays are allocated once per call and reused across
stages; a persistent workspace can be added later if profiling warrants it.

Run `make test-advection-unified` for legacy comparisons and conservation checks
on uniform and stretched grids with both velocity signs.
