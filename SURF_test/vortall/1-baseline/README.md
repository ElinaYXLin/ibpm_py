# Cylinder flow at Re=100 vs. `VORTALL.mat`

This directory validates the Python IBPM port (`py/ibpm.py`) against
`VORTALL.mat` — a reference vorticity dataset (89351 × 151: a
(449×199) vorticity field snapshotted 151 times) for 2D flow past a
circular cylinder at Re=100. All figures here are generated purely by
[`gen_vortall_report.py`](gen_vortall_report.py)/matplotlib from actual
simulation output and the `.mat` file — nothing is hand-drawn or
AI-generated.

## How the comparison was set up

`VORTALL.mat`'s snapshot count (89351 = 449 × 199) matches exactly what
`py/ibpm.py`'s own binary restart format produces for a grid with
`nx=450, ny=200` (interior vorticity nodes are stored as
`(nx-1)×(ny-1) = 449×199`, see `VALIDATION_README.md`). Combined with
the standard domain used for this classic dataset (length 9 in x,
length 4 in y, dx=0.02 → x ∈ [-1, 8], y ∈ [-2, 2], cylinder of diameter
1 centered at the origin), this pins down the run, driven by
[`run_vortall.py`](run_vortall.py):

```bash
python3 SURF_test/vortall/run_vortall.py
```

which runs both `build/ibpm` (into `_run_data_cpp/`) and `py/ibpm.py`
(into `_run_data/`) out to step 14000 (t=280), each resuming *only* from
its own last restart file if run again (e.g. after an interruption).

### Correction: an earlier version of this run seeded Python from C++'s own output

An earlier version of this dataset advanced the Python trajectory's
final 100 steps with `-ic SURF_test/vortall/_run_data_cpp/vortall13900.bin`
— i.e., using the **C++** run's own restart file as Python's initial
condition, instead of Python's own. That meant 13900 of the reported
"Python" run's 14000 steps (99.3%) were actually C++-computed, and
Python only independently integrated the last 100 steps (2 of 280 time
units) from a state it did not itself produce. It also silently
overwrote a fully-independent Python trajectory that a prior, correct
multi-chunk run (`_run_data_log.txt`/`log2.txt`/`log3.txt`) had already
completed on its own.

This made the headline "Python vs. C++ agree to 5e-12" claim close to
tautological — two runs sharing 99.3% identical history and diverging
only over a 100-step tail from an *identical* starting state will
trivially agree at the end — and it left `_run_data/vortall.force` with
only 101 rows (t=278–280) even though `shedding_summary.txt` and
`force_coefficients_saturated.png` (generated earlier, from the correct
run, and never regenerated afterward) report statistics computed over
t=200–280. The two were quietly inconsistent with each other.

**Fixed** by [`run_vortall.py`](run_vortall.py): its `_last_checkpoint()`
only ever globs restart files inside the implementation's *own* output
directory, so it is structurally impossible for one implementation's
continuation to pick up the other's restart file, even by accident. All
figures/tables in this directory were regenerated from a from-scratch
rerun (both implementations independently integrated from t=0 through
t=280, no cross-seeding) using this script.

**Local environment note:** on this machine, `python3 -m py.ibpm`
fails with `ModuleNotFoundError: __path__ attribute not found on 'py'`
— an unrelated PyPI package literally named `py` (an old pytest
dependency) is installed in this Python's `site-packages` and shadows
this repo's `py/` namespace package. This is an environment quirk, not
a bug in the repo. Worked around by forcibly registering this repo's
`py/` directory as `sys.modules["py"]` before importing, instead of
using `-m`; see the runner used for these runs, reconstructed here:

```python
import sys, types, pathlib
repo_root = pathlib.Path("/path/to/ibpm_py-main")
sys.path.insert(0, str(repo_root))
pkg = types.ModuleType("py")
pkg.__path__ = [str(repo_root / "py")]
sys.modules["py"] = pkg
from py.ibpm import main
sys.exit(main(["py.ibpm", ...your -flags...]))
```

The reshape order for `VORTALL.mat` columns was determined empirically:
`VORTALL[:,k].reshape(449, 199, order='C')` reproduces a clean,
recognizable cylinder wake (row-major storage, matching this port's own
`Scalar._data[i, j]` convention with `i` the x-index); `order='F'`
produces a garbled non-physical field and was discarded after a visual
sanity check (see the four-way disambiguation this script's development
ran, not kept in this directory).

## Why the run needed to go out to t=280

With a symmetric geometry and a zero initial condition, `py/ibpm.py`
(like the C++ original) starts with an exactly symmetric flow. The
Re=100 wake instability is real but the only thing that seeds it is
floating-point roundoff, so it grows *exponentially slowly* from
~1e-14: `Cl` stayed below 1e-6 through t=80, crossed 0.5 by about
t=135, and reached its full periodic limit-cycle amplitude (`Cl` peak =
0.846, unchanging cycle-to-cycle) by about t=160 in this run (the exact
growth timescale depends on the fine details of the floating-point
roundoff pattern, so it's expected to differ slightly, e.g. between the
C++ and Python runs — see "Correction", below). t=200 is used as the
"safely saturated" cutoff for all statistics/plots in this directory,
comfortably past that. This is expected, physically-correct behavior
for an unperturbed impulsive start, not a bug — see `shedding_summary.txt`
and `force_coefficients_saturated.png`.

**Caveat on the comparison:** `VORTALL.mat`'s snapshots are *not*
time-aligned with this run (its own snapshot 0 already shows a fully
saturated wake, so it was not started the same way / it's an excerpt
of a longer run). The comparison below is therefore of the periodic
*state* (same wavelength, same vortex spacing/strength, same
qualitative pattern), not of matching timestamps.

## Files

| File | What it is |
|---|---|
| `run_vortall.py` | Drives both `build/ibpm` and `py/ibpm.py` independently out to t=280, each resuming only from its own restart files (see "Correction" above). Run this first. |
| `gen_vortall_report.py` | Script that produces the Python-only figures. |
| `gen_three_way.py` | Script that produces the VORTALL/C++/Python three-way figures below (needs a C++ `build/ibpm` reference run in `_run_data_cpp/`, see next section). |
| `vorticity_comparison.png` | `VORTALL.mat` snapshot 150 vs. this run's saturated-shedding vorticity field (t=280), same colorbar/domain. Visually near-identical vortex street. |
| `vorticity_comparison_3way.png` | Same, but with a third row: this repo's **C++** reference build (`build/ibpm`), run on the identical grid/Re/domain, also advanced to t=280. |
| `python_vs_cpp_diff.png` | C++ vs. Python vs. their pointwise difference, at t=280. |
| `three_way_summary.txt` | Peak-vorticity and Python-vs-C++ agreement numbers backing the "why are Python and VORTALL different" answer below. |
| `flow_evolution_python.png` | Six vorticity snapshots (t=200..280) from this run's periodic regime, showing the wake pattern advecting/shedding self-consistently over multiple periods. |
| `force_coefficients_saturated.png` | `Cd(t)`, `Cl(t)` over t=200-280: `Cl` oscillates with constant amplitude cycle-to-cycle — direct evidence the flow has reached its periodic limit cycle, not still transiently evolving. |
| `shedding_summary.txt` | Measured shedding period/Strouhal number from `Cl` peak spacing, with a literature comparison. |
| `_run_data/`, `_run_data_cpp/` | `.force`/`.cmd`/`.cholesky`/`run_log.txt` files from the Python and C++ runs respectively (`.bin` and `.cholesky` are gitignored — see `.gitignore` — so those aren't committed; `.force`/`.cmd`/`run_log.txt` are small and are). The restart `.bin` snapshots used to make the figures above were deleted afterward to save local disk space — rerun `python3 SURF_test/vortall/run_vortall.py` to regenerate a full trajectory from scratch; the `.cholesky` cache means a rerun on this grid/geometry/dt skips the slow factorization step. |

## Why do the Python and VORTALL.mat panels look different?

Short answer: **they're not different because of a Python-port bug** —
this repo's C++ reference build (`build/ibpm`, compiled from `src/`,
the actual reference implementation) was run independently on the exact
same grid/Re/domain, from its own zero initial condition (no
cross-seeding with the Python run; see "Correction" above). Its final
snapshot does **not** match Python's pointwise (max |diff| = 7.4,
comparable to the field's own ~24 peak magnitude — see
`python_vs_cpp_diff.png` and `three_way_summary.txt`), and that's
*expected*: the Re=100 wake instability here is seeded only by
floating-point roundoff (see previous section), so two independently-run
trajectories are chaotically sensitive to fp-level differences during
the transient growth phase and end up at different, uncorrelated *phases*
of the same periodic cycle by t=280 — like two identical pendulum clocks
started a fraction of a second apart, still ticking at the same rate
but no longer showing the same second. The right comparison for two
independently-seeded chaotic runs is therefore the periodic *state* they
converge to, not a pointwise snapshot at a fixed time — and there, C++
and Python agree to 4-5 significant figures (shedding period, Strouhal
number, peak `Cl`, mean `Cd`; see `three_way_summary.txt`), exactly as
the exact (non-chaotic) agreement already shown for the 200×200 case in
`SURF_test/built_in_tests/README.md`. So whatever's different from `VORTALL.mat` is
different from *both* C++ and Python (which agree with each other in
periodic state) equally, and traces to a real mismatch with the
(unknown, undocumented) parameters used to generate `VORTALL.mat`
itself, not to this codebase. Two concrete, measured differences:

1. **Peak vorticity magnitude**: `VORTALL.mat`'s snapshot 150 has
   max |ω| ≈ 18.1; this repo's C++ and Python runs reach max |ω| ≈ 23.7
   and ≈ 24.5 respectively at t=280 (the two codes' own instantaneous
   peaks differ by a few percent from each other for the same phase-drift
   reason as above — a snapshot-to-snapshot vorticity extremum is a
   pointwise, phase-sensitive quantity, unlike the periodic-state
   statistics in `three_way_summary.txt`) — both around 30-35% higher
   than `VORTALL.mat`. Plausible causes (not independently confirmed, since
   `VORTALL.mat` ships with no metadata beyond the array itself):
   the reference dataset may have been produced with a coarser `dt`
   or different multi-domain (`ngrid`) far-field treatment (this run
   uses `ngrid=1`, a single fixed-boundary domain, which is the
   simplest/least-accurate far-field BC this code supports), or the
   dataset may have been lightly smoothed/interpolated when packaged
   for teaching use.
2. **Snapshot phase isn't aligned**: `VORTALL.mat`'s snapshot 0 is
   *already* in the saturated periodic regime (no transient), so its
   151 snapshots are some arbitrary excerpt of a longer run — there's
   no way to line up "snapshot 150" with "t=280" of a fresh run by
   timestamp. The comparison is necessarily of the periodic *state*
   (wavelength, vortex spacing, alternating-sign pattern), not
   matching instants.

Despite both of these, the qualitative agreement is strong: same
wavelength, same vortex-core spacing, same staggered alternating-sign
pattern, same downstream decay envelope — see
`vorticity_comparison_3way.png`.

## Result

- **Qualitative match**: `vorticity_comparison_3way.png` shows the same
  staggered, alternating-sign vortex street, same approximate vortex
  core spacing and downstream decay envelope, across all three fields
  (`VORTALL.mat`, C++, Python).
- **Python == C++ (this codebase)**: NOT a pointwise match at t=280 (max
  |diff| = 7.4 — expected, see "Why do the panels look different?"
  above) but agreement to 4-5 significant figures in every periodic-state
  statistic — shedding period, Strouhal number, peak `Cl`, mean `Cd`
  (`three_way_summary.txt`) — which is the correct evidence the port is
  faithful for a chaotically-sensitive case, consistent with the exact
  (non-chaotic) agreement in the 200×200 validation in `SURF_test/built_in_tests/README.md`.
- **Quantitative check vs. literature**: this run's Strouhal number,
  St ≈ 0.215 (period ≈ 4.655, see `shedding_summary.txt`), is higher
  than the unbounded-domain literature value (St ≈ 0.164-0.17 at
  Re=100). This is explained by blockage: the domain height here is 4
  (cylinder diameter 1 → 25% blockage ratio), and confined-flow
  blockage is known to raise shedding frequency. `VORTALL.mat` —
  generated on the same 449×199/domain-size grid — should have the
  same blockage effect, so this is a property of the shared domain
  size, not a Python-port discrepancy.

## Building the C++ reference (`build/ibpm`)

This checkout was missing `build/Makefile` entirely (not just its
gitignored build artifacts). It was reconstructed from `test/Makefile`'s
expectations (`BUILDDIR=../build`, target `libibpm.a`) and `src/`'s own
file layout (`ibpm.cc`/`checkgeom.cc` are the two driver `main()`s,
everything else compiles into the library). `make` from the repo root
now builds `build/ibpm` and `build/checkgeom` successfully (FFTW3 found
at `/usr/local/include`/`/usr/local/lib` on this machine).
