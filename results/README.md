# Validation results: Python port vs. C++ reference

This directory contains evidence that the Python port in `py/` produces
the same physics as the original C++ code in `src/`, generated entirely
by [`generate_validation_report.py`](generate_validation_report.py) — a
script that re-runs *both* implementations, on every valid example
geometry in this repository, and compares their actual output. No figure
or number here is hand-drawn, estimated, or AI-generated; everything
traces back to `build/ibpm` and `python3 -m py.ibpm` output files read
with `numpy`/`matplotlib`.

**If you read one section, read "A previous version of this report
claimed Python was faster than C++ — here's why that was misleading"
below.** That correction is the main thing that changed since the first
version of this report.

## Which cases are being run?

Both implementations are run on **real example geometry files that ship
with this repository** — not synthetic test cases. This repo contains
four `.geom` files; three are used, one is excluded (and the exclusion is
reported by the script, not silently skipped):

| Case | Geometry file | Points | Grid | dt | Notes |
|---|---|---|---|---|---|
| `cylinder` | `examples/cylinder.geom` | 160 | 200x200 | 0.02 | the standard example distributed with this project |
| `ibpm_geom` | `ibpm.geom` | 314 | 400x400 | 0.01 | a finer boundary, originally sized to match a Fortran reference case |
| `cylinder2pa` | `benchmarking/cylinder2Pa.geom` | 160 | 200x200 | 0.02 | same body as `cylinder.geom`, but with an explicit `motion fixed 0 0 0` command — exercises `RigidBody`'s motion-parsing code path, which `cylinder.geom` doesn't |
| ~~`cylinder2paplunge`~~ | `benchmarking/cylinder2PaPlunge.geom` | — | — | — | **excluded**: its `motion PitchPlunge 0 0 0.5 0.2` line supplies only 4 of PitchPlunge's 6 required parameters (a malformed file, not a porting issue), and PitchPlunge motion itself isn't ported to Python yet (see `py/fixed_position.py`'s module docstring). `Geometry.load()` returns `False` for this file in `py/` — it currently cannot be run through the Python port at all. |

All three valid cases use Re=100, RK3 timestepping, 250 steps, run from a
zero initial condition. `ibpm.geom` uses `dt=0.01` instead of `0.02`
because of a real, geometry-independent finding: **at `dt=0.02`, the
finer 400x400 grid blows up to NaN around t=0.3 — identically, to the
last digit, in both C++ and Python** (verified: both diverge at the same
timestep with matching exploding values, e.g. both reach $C_d
\approx -3.4\times10^{78}$ at the same step before going to NaN one step
later). This is a genuine CFL/stability limit of the explicit-viscous
fractional-step scheme at that resolution, not a disagreement between the
two codes, so `dt=0.01` (stable for the full run) is used for that case
instead. This is exactly the kind of behavior a correct-but-unstable port
should show: **matching failure**, not silently-different failure.

Earlier debugging (see `py/cholesky_solver.py`'s module docstring)
separately found that `cylinder.geom` on a *64x64* grid is a different,
unrelated degeneracy (an over-resolved boundary makes the projection
matrix singular, producing `NaN` in both C++ and Python) — the 200x200
resolution used here does not hit that case.

## Files

| File | What it is |
|---|---|
| `generate_validation_report.py` | The script that produces everything below. Re-run it with `python3 results/generate_validation_report.py` (requires a local C++ build at `build/ibpm`; see the top-level `README`). Takes ~2 minutes total for all three cases. |
| `validation_metrics.csv` | The full per-case error/runtime table, in machine-readable form. |
| `figures/<case>/*.png` | Five per-case figures (force comparison, force error, vorticity comparison, vorticity parity, flow evolution) — described below. |
| `figures/all_cases_force_comparison.png` | **Stitched**: all three cases' $C_d(t)$ curves, stacked in one figure. |
| `figures/all_cases_vorticity_comparison.png` | **Stitched**: all three cases' final-snapshot vorticity fields (C++ / Python / diff), stacked in one figure. |
| `figures/runtime_phase_breakdown.png` | **Stitched**: runtime for all three cases, broken into 3 phases (see below). |
| `_run_data/` | Raw output (restart `.bin` files, `.force` files, per-run timestamped logs) from the runs. Not committed to git (regenerated each time the script runs; see `.gitignore`). |

## A previous version of this report claimed Python was faster than C++ — here's why that was misleading

The first version of this report measured *total cold-start wall-clock
time* and found Python finishing in ~40% of the C++ time. That number was
real, but comparing it head-on was misleading, and it does **not** mean
"the Python port's numerics are faster than the C++ code's." Here's what
was actually happening:

`src/EllipticSolver2d.cc` builds its FFTW plan with
`fftw_plan_r2r_2d(..., FFTW_EXHAUSTIVE)`. `FFTW_EXHAUSTIVE` doesn't just
pick a reasonable FFT algorithm — it *times every strategy it knows* for
that exact transform size and picks the fastest, which is extremely slow
the first time a given size is requested (and fast forever after, since
FFTW caches that decision in-process). This plan gets built **inside the
constructor** of every `PoissonSolver`/`HelmholtzSolver`, and
`NavierStokesModel` (1x) and each RK3 substep's `ProjectionSolver` (3x)
each construct one — 4 total per run. We confirmed this directly by
timestamping every line of program output:

```
cylinder case, C++:
  0.00s  process start
  9.05s  NavierStokesModel finishes constructing (this line is where "Using
         Cholesky solver..." first prints) <- 9.05s spent, almost entirely
         in ONE FFTW_EXHAUSTIVE plan search
 11.06s  all 3 ProjectionSolvers finished (3 more FFTW_EXHAUSTIVE plans,
         but now nearly instant -- FFTW reuses its cached "wisdom" for
         same-size transforms)
 14.00s  250 timesteps complete
```

Python's `scipy.fft` has no equivalent planning phase at all, so it never
pays this cost — which is why the *first* version of this report's "total
time" comparison favored Python so heavily. That's a real difference
between the two builds, but it's a statement about FFTW's algorithm-search
behavior, not about which language/implementation computes a timestep
faster.

**The fix**: this version of the script times three phases separately
(using the same progress messages both programs print, so no code was
changed to measure this):

1. **model construction** — dominated by the one-time FFTW_EXHAUSTIVE
   search in C++; near-instant in Python.
2. **solver factorization** — building/factoring the projection operator
   for each RK3 substep (Cholesky factorization + `scipy`/FFTW calls that
   reuse cached wisdom in C++). Comparable between the two.
3. **timestepping** — the actual N-step time-integration loop. **This is
   the number that reflects the Python port's computational performance**,
   and it is close to parity with C++ across all three cases:

| Case | C++ ms/step | Python ms/step | Python/C++ |
|---|---|---|---|
| cylinder | 11.5 | 12.1 | 1.05x |
| ibpm_geom | 48.5 | 47.6 | 0.98x |
| cylinder2pa | 11.2 | 12.0 | 1.08x |

i.e. Python is **within ~5-8% of C++ per timestep** (sometimes very
slightly slower, once very slightly faster — consistent with run-to-run
noise, not a systematic gap in either direction). See
`figures/runtime_phase_breakdown.png` for the visual breakdown, and the
`ms/step` rows in `validation_metrics.csv` for exact numbers from the
most recent run.

**Takeaway for a presentation:** don't quote a single "Python vs. C++
runtime" number without saying which phase it's for — the honest, useful
claim is "the Python port's per-timestep cost is within ~10% of the
compiled C++ code," not "Python is 2-3x faster," which was an artifact of
one C++ file's FFTW planning mode.

## Summary table (see `validation_metrics.csv` for the full version)

| Case | $C_d$ max rel. error | $C_l$ max abs. error | Vorticity max abs. error | Vorticity max error / field peak |
|---|---|---|---|---|
| cylinder | 0.000e+00 | 5.19e-14 | 6.14e-12 | 2.66e-13 |
| ibpm_geom | 0.000e+00 | 1.12e-13 | 9.04e-12 | 3.50e-13 |
| cylinder2pa | 0.000e+00 | 5.19e-14 | 6.14e-12 | 2.66e-13 |

**Note on $C_l$:** for these symmetric geometries (no angle of attack),
lift is exactly zero by symmetry, so both implementations produce
$C_l \approx 10^{-14}$ (floating-point roundoff of an exact zero, not a
physical result). *Relative* error against a ~$10^{-14}$ reference is
meaningless (dominated by noise), so *absolute* error is reported instead
— and even that is itself only floating-point-roundoff-sized.

**Why $C_d$'s error is exactly `0.000e+00`, not just "small":** the two
`.force` files are byte-identical in that column across all 251 rows in
every case (`np.abs(cd_py - cd_cpp).max() == 0.0`, verified, not rounded
for display). The two codes execute numerically identical floating-point
operations, in the identical order, for the entire force-computation
path — evidence the arithmetic wasn't just reorganized during porting.

## Figures (per case, under `figures/<case>/`)

### `force_coefficients_vs_time.png`
$C_d(t)$ and $C_l(t)$ from both implementations, overlaid (Python dashed
on top of C++ solid). The headline "do these two codes agree" plot: the
drag curve shows the expected impulsive-start transient, decaying to a
roughly steady value, and the two lines are visually indistinguishable
throughout. The lift panel's y-axis is scaled to $10^{-14}$ — that panel
is a picture of floating-point noise, not a physical signal.

### `force_error_vs_time.png`
The dense, per-timestep validation signal: $|C_{d,\mathrm{py}} -
C_{d,\mathrm{cpp}}|$ and $|C_{l,\mathrm{py}} - C_{l,\mathrm{cpp}}|$ vs.
time, log-scale. Answers "does disagreement grow over time" (e.g. from
accumulating floating-point drift) — it does not: $C_d$ error is exactly
zero at every timestep (the flat line at $10^{-18}$ is a plotting floor
so an exact zero is visible on a log axis, not a measured value), and
$C_l$ error stays flat around $10^{-14}$–$10^{-15}$ with no upward trend.

### `vorticity_field_comparison.png`
A 3x3 grid: rows are snapshots at $t=0$ and the final time; columns are
C++, Python, and their difference (its own, much smaller color scale).
The C++ and Python columns are visually identical; the difference column
shows $O(10^{-12})$ deviations — about 13 orders of magnitude below the
field's own peak magnitude, i.e. floating-point roundoff, not a
systematic discrepancy.

### `vorticity_parity_plot.png`
Python vorticity vs. C++ vorticity, one point per grid node, at the final
snapshot, with a $y=x$ reference line. Perfect agreement puts every point
exactly on the line, which is what's shown ($R^2 = 1.0000000000$ to 10
decimal places) — and it holds uniformly across the *entire* range of
values, not just on average.

### `flow_evolution_python.png`
Not a comparison — a standalone physical-interpretation figure (Python
output only) showing the vortex pair developing behind the body as the
impulsively-started flow evolves. Useful to show the port produces
physically sensible results, independent of the C++ comparison.

## Stitched (all-cases) figures

### `all_cases_force_comparison.png`
All three cases' $C_d(t)$ curves stacked vertically — one glance to see
that the C++/Python agreement holds across every geometry tested, not
just the one headline example.

### `all_cases_vorticity_comparison.png`
All three cases' final-snapshot vorticity fields (C++ / Python / diff),
stacked. Note the per-case final *time* differs (labeled in each
subplot title) since `ibpm_geom` uses a smaller `dt` for stability (see
above) — the step count (250) is the same, the physical time reached
isn't.

### `runtime_phase_breakdown.png`
Stacked bars (grey = model construction, light blue = solver
factorization, dark teal = timestepping) for C++ (solid) and Python
(hatched), one pair per case. This is the figure that makes the FFTW
finding visible at a glance: C++'s grey segment dwarfs Python's in every
case, while the dark teal ("timestepping") segments are comparable —
sometimes C++ is taller, sometimes Python is.

## Regenerating these results

```bash
cd build && make          # build the C++ reference, if not already built
cd ..
pip install -r py/requirements.txt matplotlib
python3 results/generate_validation_report.py
```

The script re-runs all three cases, for both implementations, from
scratch each time (no caching of prior results) — about 2 minutes total
on the machine this was last run on. It will always produce
numbers/figures consistent with the current state of `src/` and `py/`.
