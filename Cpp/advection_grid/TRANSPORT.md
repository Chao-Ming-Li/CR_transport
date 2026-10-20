Both transport calculations use an instantaneous burst at t=0:

    n(R,z,0) = exp(-((R-10 kpc)^2 + (z-10 kpc)^2)/(9 kpc^2))
    D = 1e29 cm^2/s = 3.31438e-7 kpc^2/yr
    v_R = 0, v_z = +300 km/s = 3.06814e-7 kpc/yr
    dt = 1e4 yr, duration = 1e8 yr

The source is axisymmetric (a Gaussian ring), with sigma_R=sigma_z=3/sqrt(2)
kpc. Injection occurs only at initialization. There are no energy losses,
spallation, decay, or continuous injection. R=0 and z=0 have zero flux;
R=200 kpc and z=200 kpc absorb particles. Monte Carlo uses Cartesian Brownian
increments and finite-step reflection at z=0; boundary errors require timestep
convergence if studying the boundary layer.

Run the grid solver with `make` then `OMP_NUM_THREADS=4 ./a.out`.
In `../../Monte Carlo`, run `make` then `OMP_NUM_THREADS=4 ./transport`.
The loss/spallation experiments remain available as explicit Makefile targets.

Grid output `ndis_C_v300_D1e29_cp.bin` is overwritten per run, contains native
float64 arrays including two ghost cells on each side, and has shape
(1001, NR+4, NZ+4). Snapshot times are 0, 1e5, ..., 1e8 yr. Physical cells are
[:,2:-2,2:-2]. `transport_grid.csv` gives the cell coordinates and widths for
both identical axes (NR=NZ=number of rows). Grid resolution has changed; readers
must use this geometry rather than older hard-coded dimensions.

Monte Carlo `particles.bin` is overwritten per run and contains native float64
triples (x,y,z) in pc at exactly 1e8 yr. There are initially 10000 particles;
absorbed particles are omitted. Sampling includes the cylindrical Jacobian:
p(R) is proportional to R*exp(-(R-10)^2/9), with uniform azimuth.

For density comparison, divide the grid density by its initial volume integral.
Divide Monte Carlo counts by 10000 and cylindrical bin volumes
pi*(R_hi^2-R_lo^2)*(z_hi-z_lo), with all positions converted to kpc. This retains
any escaped fraction rather than renormalizing the survivors.

Away from the reflecting plane, the expected vertical displacement is 30.6814
kpc and the diffusive variance added in each Cartesian coordinate is 66.2877
kpc^2. The radial second moment grows by 4*D*t = 132.5754 kpc^2. Reflection
makes small corrections to the free-space vertical moments.
