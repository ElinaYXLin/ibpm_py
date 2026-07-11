# SURF_test/vortall/inner -- residuals and inner-workings audit

`../README.md`, `../three_way_summary.txt`, and `../shedding_summary.txt`
validate vortall's *output* (force coefficients, shedding period, Strouhal
number, qualitative vorticity field). They never check whether the
linear-algebra machinery producing that output is actually solving the
equations it claims to solve at each timestep. This directory does that: it
re-runs the vortall cylinder case (Re=100, nx=450, ny=200, same geometry/dt
as `run_vortall.py`) through the internal `py/` API with instrumentation
added, and directly measures the residuals of the equations the solver is
supposed to satisfy at every one of 10,000 timesteps (t=0 to t=200 -- through
setup, the exponential instability growth, and into saturated periodic
shedding).

## Files

| File | What it is |
|---|---|
| `compute_residuals.py` | Runs the instrumented simulation, writes `data/residuals.csv`. |
| `plot_residuals.py` | Reads that CSV, writes the four figures below and `residual_summary.txt`. |
| `data/residuals.csv` | Raw per-step data (10,001 rows): every residual/quantity below, every step. |
| `data/run_meta.txt` | Grid/physics parameters and timing for the run that produced the CSV. |
| `residual_summary.txt` | Numeric summary (max/median of each residual over the run). |
| `figures/01_solver_residuals.png` | Divergence / no-slip / Poisson residuals vs. t. |
| `figures/02_vorticity_smoothness.png` | `\|\|Laplacian(omega)\|\|` vs. t (sanity check, not a residual). |
| `figures/03_conservation.png` | Circulation / enstrophy / kinetic energy vs. t. |
| `figures/04_force_reproducibility.png` | This run's Cd/Cl vs. the archived `_run_data/vortall.force`. |

Regenerate with:
```bash
python3 SURF_test/vortall/inner/compute_residuals.py --nsteps 10000
python3 SURF_test/vortall/inner/plot_residuals.py
```
(~5.7 minutes for the 10,000-step run on this machine, ~34ms/step; see `data/run_meta.txt`.)

## What was measured, and why each one matters

Four residuals of equations the solver is supposed to satisfy exactly (up to
floating-point roundoff), each targeting a different module:

1. **Divergence of the flux field**, `||div(q)||` (`figures/01`, top panel).
   Computed with a finite-difference formula written from scratch in
   `compute_residuals.py` -- it does **not** call `vector_operations.Curl` or
   any other `py/` code. A curl-derived flux field is divergence-free by
   construction on this staggered grid; if this were not at machine
   precision, `Curl()` (`flux.py`/`vector_operations.py`) or the projection
   step (`ib_solver.py`/`projection_solver.py`) would be broken.
   **Result: max 6.9e-18 over the whole run** (`float64` machine epsilon is
   ~2.2e-16) -- at the theoretical floor. **Clean.**

2. **No-slip constraint residual**, `||C(omega) - b||` (`figures/01`, middle
   panel) -- literally the equation `ProjectionSolver.solve()`
   (`projection_solver.py`) is supposed to satisfy on the cylinder surface
   every substep. Implicates `CholeskySolver`
   (`cholesky_solver.py`, the solver actually selected here, since the
   cylinder is stationary) or `Regularizer` (`regularizer.py`, which builds
   the interpolation/regularization operators `C`/`B` used inside it) if
   large. **Result: max 1.6e-14, median ~2e-15, flat for the entire 200 time
   units** including through the exponential instability growth and into
   saturated shedding. **Clean** -- the no-slip BC is enforced to
   floating-point precision at every single step, not just on average.

3. **Streamfunction Poisson-equation residual**, `||Laplacian(psi) -
   (-omega)||` (`figures/01`, bottom panel) -- whether
   `NavierStokesModel.vorticityToStreamfunction` (`navier_stokes_model.py`)
   actually returns a streamfunction that solves the Poisson equation it
   claims to. Implicates `PoissonSolver`/`EllipticSolver`
   (`elliptic_solver.py`) if large. **Result: max 5.5e-12, median ~4e-12** --
   about two orders of magnitude above the previous two (expected: this is
   an FFT-based direct elliptic solve, which accumulates more rounding error
   per call than the simpler no-slip/divergence checks, but is still 8
   orders of magnitude below anything that would indicate a real problem).
   **Clean.**

4. **Vorticity-field smoothness**, `||Laplacian(omega)||` (`figures/02`) --
   not a residual of any equation (there's nothing to compare it to), just a
   direct evaluation that flags grid-scale checkerboard noise or blow-up,
   the kind of thing an actual `Regularizer`/projection bug would produce.
   **Result:** smooth, bounded growth as the wake develops (~2e3 to ~5e3),
   no discontinuities or spikes anywhere in the 200-time-unit window.
   **Clean.**

Plus three conservation/physical-sanity quantities (`figures/03`) -- not
expected to be constant, just physically sensible: **circulation** stays at
~0 (symmetric flow) until the instability visibly breaks symmetry around
t~110, then oscillates with growing then constant amplitude exactly as
vortex shedding predicts; **enstrophy** and **kinetic energy** both rise
during the initial transient, plateau, then start oscillating once shedding
saturates around t~130 -- all standard, expected behavior for this flow, no
red flags.

## The reproducibility check turned up something worth explaining (not a bug)

`figures/04` compares this script's independently-written driver against
the archived `SURF_test/vortall/1-baseline/_run_data/vortall.force` (produced by
`run_vortall.py` via `SURF_test/cost/run_ibpm_case.py`). Both run the exact
same deterministic algorithm from the same zero initial condition, so in
principle they should be identical -- but they're two different Python
processes/driver code paths, and per `../README.md`, this Re=100 case's
wake instability is seeded only by floating-point roundoff, so any
last-bit-level difference between the two code paths is expected to be
chaotically amplified over time, not stay at zero.

That is exactly what `figures/04`'s bottom panel shows: `|Cl_this -
Cl_archived|` starts at ~1e-16 (machine epsilon) at t=0, and grows roughly
exponentially, crossing into the 1e-6-ish range by around **t~110-125** --
matching, independently, the *same* instability growth timeline
`../README.md` reports from a completely different method (watching `Cl`
itself cross 0.5 by t~135). `|Cd_this - Cd_archived|` is flatter, sitting
near 5e-6 for almost the entire run -- this ceiling is fully explained by
`OutputForce` (`output_force.py`) writing `vortall.force` with `"%.5e"`
(6-significant-figure) text formatting, i.e. it's a text-precision floor,
not evidence of a growing discrepancy. Neither difference exceeds 1e-4 (a
threshold chosen to be far above print-truncation noise but far below an
O(1) real phase split) anywhere in the t=0-200 window -- see
`residual_summary.txt`.

In short: this is independent, unplanned corroboration of the chaotic-growth
story already documented in `../README.md`, arrived at via a completely
different measurement (two independent driver processes' pointwise force
difference, rather than watching one run's `Cl` amplitude) -- not a bug.

## Conclusion

**No re-audit needed for vortall's own code path.** Every equation the
solver claims to satisfy in this run -- flux divergence-free-ness, the
cylinder's no-slip boundary condition, and the streamfunction's Poisson
equation -- holds to floating-point precision at **every one of 10,000
consecutive timesteps**, not just in aggregate or at the final step. That
directly exercises and clears: `vector_operations.Curl`/`Laplacian`,
`ib_solver.IBSolver`/`NonlinearIBSolver`, `projection_solver.ProjectionSolver`,
`cholesky_solver.CholeskySolver`, `elliptic_solver.PoissonSolver`, and
`regularizer.Regularizer`'s one-time (stationary-body) evaluation. Combined
with `../three_way_summary.txt`'s independent macroscopic agreement (C++ vs.
Python periodic-state statistics to 4-5 significant figures) and this
report's microscopic, per-timestep residual check, the mentor's skepticism
about vortall specifically should be resolved: there is no numerical
evidence of a porting bug anywhere in the code path this test exercises.

**But there IS a real, concrete gap worth re-auditing -- just not in
vortall.** Checking every test driver in this repository
(`SURF_test/vortall/1-baseline/run_vortall.py`, `SURF_test/airfoil_driver.py`,
`SURF_test/built_in_tests/generate_validation_report.py`,
`SURF_test/cost/backends.py`) shows **every single one passes `-ngrid 1`**
-- no test anywhere in this repository runs the multi-domain grid path. And
of the geometry files available, only `benchmarking/cylinder2PaPlunge.geom`
specifies genuine time-varying body motion (`motion PitchPlunge 0 0 0.5
0.2`) -- and per
`SURF_test/built_in_tests/generate_validation_report.py`'s own docstring,
that file is **malformed** (`Geometry.load()` returns `False` in both
implementations) and is excluded from every comparison. `cylinder2Pa.geom`'s
`"motion fixed 0 0 0"` only exercises `RigidBody`'s motion-*parsing* code
(confirmed by reading `FixedPosition.isStationary()`, `fixed_position.py:37`
-- it hardcodes `return True`), not actual moving-body dynamics.

The practical consequence: **`ConjugateGradientSolver`
(`conjugate_gradient_solver.py`) -- the solver `IBSolver.createSolver`
(`ib_solver.py:153-165`) switches to the moment any body actually
moves -- has never been exercised by an end-to-end numerical run anywhere in
this repository.** Neither has `Regularizer.update()`'s *repeated*,
per-timestep re-evaluation (`navier_stokes_model.py:110-116`, triggered only
when `geTimeDependent()` is true) -- every test today calls it exactly once,
at initialization, because every geometry currently runnable is stationary.
Nor has `elliptic_solver_2d.py`'s coarse-to-fine boundary-condition transfer
between grid levels, since `ngrid=1` everywhere.

### Most suspicious module: `regularizer.py`

Not because this audit found anything wrong with it here -- its one-time,
stationary-body evaluation is exactly what the near-machine-precision
no-slip residual above validates -- but because:
1. It has a **documented history of a real bug**: commit `0c68193` ("Fix
   Regularizer RAM blowup: match C++'s per-point locality, not just its loop
   structure") found and fixed an actual correctness/resource issue in this
   exact file.
2. Its **repeated-update code path is completely untested end-to-end**
   (see above) -- the one-time path this audit validated is not the same
   code path a moving-body run would exercise every timestep.
3. It is the module most directly responsible for the fidelity of the
   immersed-boundary coupling (`toBoundary`/`toFlux`), so a bug here would
   silently produce plausible-looking but physically wrong forces on a
   moving body -- exactly the kind of failure mode that would not show up
   in vortall (stationary cylinder) but would in, e.g., a future pitching-
   airfoil or flapping-wing case.

**Recommendation:** before trusting results from any *moving-body* or
*multi-domain-grid* case with this port, build/repair a runnable
moving-body test case (e.g. fix `benchmarking/cylinder2PaPlunge.geom`'s
malformed `motion PitchPlunge` line, which needs 6 parameters and only
supplies 4 -- see `py/pitch_plunge.py` for the expected signature) and
run an audit of this same shape (per-timestep no-slip/divergence/Poisson
residuals, plus a Cholesky-vs-ConjugateGradient cross-check by forcing both
solvers on the same stationary case, since `CholeskySolver` is already
validated here and the two should agree to `ConjugateGradientSolver`'s
convergence tolerance) against `regularizer.py`'s repeated-update path,
`conjugate_gradient_solver.py`, and (if a multi-domain case is added)
`elliptic_solver_2d.py`'s coarse/fine boundary transfer -- none of which
this audit, or any existing test in this repository, currently covers.
