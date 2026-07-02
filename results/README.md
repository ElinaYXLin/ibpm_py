# Validation results: Python port vs. C++ reference

This directory contains evidence that the Python port in `py/` produces
the same physics as the original C++ code in `src/`, generated entirely
by [`generate_validation_report.py`](generate_validation_report.py) — a
script that re-runs *both* implementations and compares their actual
output. No figure or number here is hand-drawn, estimated, or
AI-generated; everything traces back to `build/ibpm` and
`python3 -m py.ibpm` output files read with `numpy`/`matplotlib`.

## What case is being run?

**Both implementations are run on the actual example geometry file,
`examples/cylinder.geom`** (a 160-point circle, diameter 1) — not a
synthetic or simplified test case. The specific configuration matches
`examples/ibpm.cmd`, the command the original project's example run was
generated with:

```
-geom cylinder.geom      (loaded from examples/cylinder.geom)
-nx 200 -ny 200 -ngrid 1 -length 4.0 -xoffset -2.0 -yoffset -2.0
-Re 100 -dt 0.02 -nsteps 250 -scheme rk3
```

i.e. flow past an impulsively-started cylinder at Reynolds number 100, on
a 200x200 grid, for 250 timesteps (5 time units), started from a zero
initial condition — exactly the standard example distributed with this
project, run at its default resolution. This is called out explicitly
because an *under-resolved* version of this same case (a coarser grid,
e.g. 64x64) is numerically degenerate — the boundary is then ~3x
over-resolved relative to the grid spacing, which makes the projection
matrix singular and produces `NaN` output in **both** the C++ and Python
code (verified independently — not a porting bug). See
`py/cholesky_solver.py`'s module docstring for the full explanation. The
200x200 case used for this validation does not hit that degeneracy.

## Files

| File | What it is |
|---|---|
| `generate_validation_report.py` | The script that produces everything below. Re-run it with `python3 results/generate_validation_report.py` (requires a local C++ build at `build/ibpm`; see the top-level `README`). |
| `validation_metrics.csv` | The summary error/runtime table, in machine-readable form (for pulling into a paper/LaTeX table). |
| `figures/*.png` | The plots, described one by one below. |
| `_run_data/` | Raw output (restart `.bin` files, `.force` files, logs) from the two runs the script performs. Not committed to git (regenerated each time the script runs) — see `.gitignore`. |

## Summary table

| Quantity | Error type | Max error | Mean error | Max error, relative to field peak |
|---|---|---|---|---|
| $C_d$ (drag coefficient) | relative, `\|py-cpp\|/\|cpp\|`, over t>0 | 0.000e+00 | 0.000e+00 | — |
| $C_l$ (lift coefficient) | absolute (see note) | 5.19e-14 | 1.38e-14 | — |
| Vorticity field, t=0 | absolute | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| Vorticity field, t=2 | absolute | 4.16e-12 | 2.51e-14 | 1.78e-13 |
| Vorticity field, t=4 | absolute | 6.14e-12 | 2.76e-14 | 2.66e-13 |
| **Runtime, C++** (`build/ibpm`) | wall-clock, cold start | 12.16 s | | |
| **Runtime, Python** (`py/ibpm.py`) | wall-clock, cold start | 4.56 s | | |
| **Runtime ratio** (Python / C++) | | 0.375 | | |

*(Exact numbers regenerate slightly run-to-run for the runtime rows only,
since those are live wall-clock measurements; see `validation_metrics.csv`
for the numbers from the most recent run. The force/field error rows are
deterministic to the last bit and do not change between runs.)*

**Note on $C_l$:** for this symmetric geometry (a circular cylinder with
no angle of attack), lift is exactly zero by symmetry, so both
implementations produce $C_l \approx 10^{-14}$ (floating-point roundoff of
an exact zero, not a physical result). *Relative* error against a
reference value of ~$10^{-14}$ is meaningless (it would be dominated by
noise), so the table reports *absolute* error instead — and even that is
itself only floating-point-roundoff-sized.

**Why $C_d$'s error is exactly 0.000e+00, not just "small":** the two
`.force` files are byte-identical in that column across all 251 rows
(`np.abs(cd_py - cd_cpp).max() == 0.0`, verified, not rounded for
display). This means the two codes are executing numerically identical
floating-point operations, in the identical order, for the entire
force-computation path — evidence that the arithmetic wasn't just
reorganized during porting.

## Figures

### `force_coefficients_vs_time.png`
$C_d(t)$ and $C_l(t)$ from both implementations, overlaid (Python dashed
on top of C++ solid). This is the headline "do these two codes agree"
plot: the drag curve shows the expected impulsive-start transient (a
spike as the flow starts, decaying to a roughly steady value around
$C_d \approx 2.1$ by $t=5$), and the two lines are visually
indistinguishable across the whole run. The lift panel shows both
implementations correctly producing $C_l \approx 0$ (as expected for a
non-lifting symmetric body) — note the y-axis is scaled to $10^{-14}$,
i.e. this panel is really a picture of floating-point noise, not a
physical signal.

### `force_error_vs_time.png`
The dense, per-timestep validation signal: $|C_{d,\mathrm{py}} -
C_{d,\mathrm{cpp}}|$ and $|C_{l,\mathrm{py}} - C_{l,\mathrm{cpp}}|$ vs.
time, log-scale. This is the plot that answers "does disagreement grow
over time" (e.g. from accumulating floating-point drift between two
differently-ordered implementations) — it does not: the $C_d$ error is
exactly zero at every one of the 251 timesteps (see note above; the flat
line at $10^{-18}$ is a plotting floor so an exact zero is visible on a
log axis, not a measured value), and the $C_l$ error stays flat around
$10^{-14}$–$10^{-15}$ for the full 5 time units, with no upward trend. If
the Python port had a subtle bug causing errors to accumulate step by
step, this is the plot that would show it.

### `vorticity_field_comparison.png`
A 3x3 grid: rows are snapshots at $t=0, 2, 4$ (steps 0, 100, 200 — the
restart files both runs write); columns are C++, Python, and their
difference. This is the main "the whole flow field matches, not just an
integrated force" figure. The C++ and Python columns are visually
identical (as they should be, given the force agreement above); the
difference column uses its own (much smaller) color scale to show that
what differences exist are $O(10^{-12})$ — about 13 orders of magnitude
below the field's own peak magnitude ($\approx 23$) — i.e.
indistinguishable from floating-point roundoff, not a systematic
discrepancy. The cylinder outline (radius 0.5, from `cylinder.geom`) is
overlaid for reference.

### `vorticity_parity_plot.png`
A scatter of Python vorticity vs. C++ vorticity, one point per grid node
(39,601 points), at the final snapshot ($t=4$), with a $y=x$ reference
line. This is a standard validation plot for papers: perfect agreement
would put every point exactly on the line, which is what's shown here
($R^2 = 1.0000000000$ to 10 decimal places). It's a complementary view to
the spatial difference plot above — it makes clear the agreement holds
uniformly across the *entire* range of vorticity values (near-zero
far-field points and large near-body values alike), not just on average.

### `runtime_comparison.png`
Total wall-clock time for a cold-start run (250 steps, including the
one-time Cholesky factorization of the projection operator, which is the
dominant fixed cost for this problem size). Measured by literally timing
`subprocess.run([...])` around each program. On this machine, the Python
port is faster (~0.375x the C++ time) — this is a genuine, reproducible
measurement, not a typo, likely because `py/cholesky_solver.py` and
`py/vector_operations.py` route their heavy array math through
`numpy`/`scipy`'s BLAS-backed operations, whereas the equivalent C++ code
uses hand-written loops. **This is a single-machine, single-build
comparison, not a general "Python beats C++" claim** — it depends on the
BLAS library linked on this machine, the C++ compiler/optimization flags
used for `build/ibpm`, and the specific operation mix of this test case;
report it as "the port is not slower" rather than a general performance
claim, unless you benchmark further.

### `flow_evolution_python.png`
Not a Python-vs-C++ comparison — a standalone physical-interpretation
figure (Python output only) showing the vortex pair developing behind the
cylinder as the impulsively-started flow evolves from $t=0$ (still, zero
vorticity) through $t=2$ and $t=4$ (a forming recirculating wake, with the
classic pair of counter-rotating vortex sheets separating off the top and
bottom of the cylinder). Useful for a presentation/paper to show the
port produces physically sensible results, independent of the numerical
comparison to C++.

## Regenerating these results

```bash
cd build && make          # build the C++ reference, if not already built
cd ..
pip install -r py/requirements.txt matplotlib
python3 results/generate_validation_report.py
```

The script re-runs both implementations from scratch each time (no
caching of prior results), so it will always produce numbers/figures
consistent with the current state of `src/` and `py/`.
