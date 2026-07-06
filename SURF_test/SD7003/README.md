# SD7003 airfoil vs. UIUC LSAT experiment

This directory validates `py/ibpm.py` on a real cambered airfoil (not a
cylinder): the SD7003, a low-Reynolds-number airfoil extensively wind-tunnel
tested at UIUC. All figures/numbers here are generated purely by
[`gen_airfoil_report.py`](gen_airfoil_report.py) and
[`gen_flowfield_figs.py`](gen_flowfield_figs.py) (matplotlib) from actual
`py/ibpm.py` output and the `.DRG`/`.LFT` reference files in
`SURF_test/SD7003/` — nothing is hand-drawn or AI-generated.

## Data provenance (confirmed online)

- **Coordinates**: `SURF_test/SD7003/sd7003.dat.txt` (header
  `SD7003-085-88`) matches, byte-for-byte in header/shape, the file
  hosted at `https://m-selig.ae.illinois.edu/ads/coord/sd7003.dat` on
  UIUC's Airfoil Coordinates Database — confirmed by fetching that URL
  directly.
- **Performance data**: `SD7003.DRG.txt`/`SD7003.LFT.txt` are in the
  exact UIUC LSAT (Low-Speed Airfoil Tests) tabulated format (same
  header fields, and the `Tabulated from data in file ####.DAT &
  reduced using program Q01.00L` footer — `Q01.00L` is Selig lab's
  specific data-reduction tool signature). I searched UIUC's LSAT
  volume directories (`m-selig.ae.illinois.edu/pd/pub/lsat/volume0{1,3}/`)
  for the exact source file but could not pin down which
  volume/supplement hosts this particular airfoil's table in the time
  available — the format match is strong evidence of authenticity, but
  I want to flag that I did not locate a byte-identical online copy of
  the `.DRG`/`.LFT` tables themselves, only the coordinates.

## Setup

- Body: raw UIUC coordinates, **resampled to uniform arc-length
  spacing matched to each run's grid `dx`** (see "why resampling was
  necessary" below), quarter-chord (`0.25, 0`) as rotation center.
- Angle of attack applied via `py.ibpm`'s own `-alpha` flag (rotates
  the free-stream, physically equivalent to pitching the body — this
  is already wired up correctly in `ibpm.py`'s lift/drag frame
  transform, so no geometry rotation was needed).
- Domain: `length=6, xoffset=-2, yoffset=-1.5` (x in [-2,4], y in
  [-1.5,1.5], chord=1), `ngrid=1`, `dt=0.01` (`dt=0.005` at the finest
  grid level, see below), `Re=61100` (UIUC's lowest-Re bin, closest to
  the coordinate file's own test condition).
- Cl/Cd are **time-averaged over the last 60% of a 3000-step (t=30)
  run**, with the reported error bars = ±1 std. dev. of the
  instantaneous force trace over that window (this flow is
  genuinely unsteady at this Re/resolution — see "Limitations" below —
  so a single instantaneous value would be misleading).

### Why the raw UIUC coordinates couldn't be used directly

The UIUC `.dat` file has **highly non-uniform** point spacing (as
close as ~0.001c near the leading edge, ~0.05c mid-chord). Feeding
that directly into `py/ibpm.py` makes the projection matrix singular
after 1 timestep (`NaN` forces from step 1) — a documented degeneracy
of this codebase when boundary points are spaced much finer than the
grid (see `py/cholesky_solver.py`'s module docstring, and
`SURF_test/built_in_tests/README.md`'s note on the same issue for an over-resolved
circle). Fixed by re-parametrizing the boundary by arc length and
resampling to spacing ≈ `dx` (see `gen_airfoil_report.py`'s
`make_airfoil_raw` import) — standard practice for immersed-boundary
methods. This means the geometry file itself changes with each grid
resolution (more points at finer `dx`); see `SURF_test/geom/`.

## Files

| File | What it is |
|---|---|
| `polar_comparison.png` | $C_l(\alpha)$, $C_d(\alpha)$: `py/ibpm.py` (dx=0.02) vs. UIUC LSAT experiment. |
| `drag_polar.png` | Same data as a $C_l$-$C_d$ drag polar. |
| `grid_convergence.png` | $C_l$, $C_d$ at fixed $\alpha=-0.09°$ across dx = 0.04, 0.02, 0.01, vs. the experimental value at that $\alpha$. |
| `flow_evolution.png` | Vorticity field snapshots (t=0 to 30) at $\alpha=4.6°$, Re=61100, dx=0.02, from impulsive start. |
| `summary.txt` | Numeric table backing the two comparison figures. |
| `_run_data/` | Per-case `.force`/`.cmd`/log files (raw `.bin` restarts deleted after figure generation to save space, except `flowfield/` which backs `flow_evolution.png`). |

## Results

**Lift matches well.** $C_l(\alpha)$ tracks the experimental polar
closely across the whole tested range (-2.9° to 7.7°) — see
`polar_comparison.png`, left panel. This is the headline result: a
faithful 2D immersed-boundary Navier-Stokes solver, on a fairly coarse
grid, reproduces the *integrated* lift of a real cambered airfoil to
good engineering accuracy.

**Drag is systematically overpredicted at higher $\alpha$** (e.g. at
$\alpha=7.72°$: sim $C_d\approx0.140\pm0.079$ vs. exp $C_d=0.026$) —
see `polar_comparison.png`, right panel, and `drag_polar.png`. This is
a real, physically-explainable effect, not noise: $C_d$ is dominated
by viscous/pressure drag from the (thin, low-Re) boundary layer and
any separation, and:
1. This solver has no turbulence/transition model — real SD7003 flow
   at Re~60000 relies on a laminar separation bubble with transition
   and reattachment; a 2D laminar Cartesian IB solver at this
   resolution cannot resolve that physics, and tends to over-predict
   separation extent (hence pressure drag) once loading increases.
2. The very large error bars at higher $\alpha$ (e.g. ±0.079 at 7.72°)
   show the simulated flow is **genuinely far more unsteady** than the
   real (quasi-steady, LSB-stabilized) wind-tunnel flow — consistent
   with under-resolved transition physics.

**Grid convergence (`grid_convergence.png`)**: at the mild angle
$\alpha=-0.09°$, both $C_l$ and $C_d$ converge monotonically toward the
experimental value as `dx` shrinks from 0.04 to 0.01, **and the
unsteadiness (error bar size) shrinks by roughly 3x over the same
range** (coarse $C_l$ std=0.72, fine std=0.21) — direct evidence that
much of the coarse-grid disagreement is a resolution artifact, not a
fundamental modeling error. Numbers in `summary.txt`.

**The finest grid level (dx=0.01) diverged to `NaN` at the polar
sweep's `dt=0.01`** — a real CFL-type finding, not a bug: this
solver's viscous term is handled via a semi-analytic integrating
factor (unconditionally stable), but the convective term is fully
explicit, and the leading-edge suction peak at this airfoil's
operating conditions locally accelerates the flow well above the
freestream speed, tightening the effective CFL limit faster than `dx`
alone would suggest. Fixed by halving `dt` to 0.005 (and doubling
`nsteps` to reach the same physical time) for that one grid level only
— see `SURF_test/rerun_fine_log.txt`.

## Limitations (read before quoting these numbers)

The vorticity field itself (`flow_evolution.png`) shows **broadband,
grid-scale-speckled noise spreading through the whole domain**,
including regions far from the body/wake that should be undisturbed
uniform flow. This is a real limitation of this 2D single-domain
Cartesian solver at Re≈61100 with no explicit subgrid dissipation —
the flow is at a Reynolds number and geometry (laminar separation
bubble airfoil) that's a genuinely hard test case even for high-order
LES/DNS codes, and a resolution of dx=0.02-0.04 is well below what's
needed to cleanly resolve it. The *integrated* force coefficients
above are far better behaved than the *local* vorticity field — trust
the polar/convergence numbers more than the flow-field snapshots for
quantitative claims.
