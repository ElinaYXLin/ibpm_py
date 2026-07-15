# Reproducibility check: the hardest case (SD8000, coarse grid, Re=60800)

Same check as [`../README.md`](../README.md) (NACA0012, Re=500), but on the
hardest configuration in `SURF_test/airfoils/`, to stress-test the
fixed-algorithm DST change (`FFTW_EXHAUSTIVE` &rarr; `FFTW_ESTIMATE |
FFTW_UNALIGNED`) under the worst conditions this repo has documented.

## Why this is the hardest setup

`SURF_test/airfoils/README.md`'s "mentor question" investigation and
`SURF_test/gen_port_fidelity_diagnostic.py` both single out **SD8000's
coarse grid (dx=0.04, Re=60800, &alpha;=-0.81&deg;)** as the one case in the
whole airfoil suite where this solver is not just under-resolved but
**numerically unstable**: the original C++ run at this resolution ends in
a catastrophic single-step blow-up (Cl spiking to ~2000). Worse,
[`SURF_test/airfoils/7-chaos_sensitivity/README.md`](../../SURF_test/airfoils/7-chaos_sensitivity/README.md)
found the *timing* of that blow-up is chaotically sensitive to last-bit
floating-point roundoff — in one pair of reruns, Python and C++ swapped
which one blew up first, traced there to `FFTW_EXHAUSTIVE`'s own
run-to-run non-reproducibility (it re-times candidate algorithms and can
select a different one under machine load, a last-bit-level perturbation
big enough to flip which implementation crosses the instability threshold
first).

That makes it the most demanding possible test of *this specific* change.
A clean/laminar case (like NACA0012 at Re=500) would reproduce trivially
regardless of whether the DST algorithm were pinned or not, because there's
no chaotic amplification to expose a difference. This case has both
properties needed to actually test the fix: (1) it's chaotic, so any
residual non-determinism gets amplified into a large, visible divergence
rather than staying buried at the roundoff floor, and (2) it's *already
documented* to have been non-reproducible for exactly the reason this
change targets.

Every other airfoil in `SURF_test/airfoils/` (SD7003, ClarkY, GM15) sits at
a similar Re (~40,600-61,100) and shows the same qualitative broadband
speckle at production resolution (`dx=0.02`), but SD8000's coarse grid
(`dx=0.04`) is the one pushed all the way to outright blow-up — a strictly
harder failure mode than "noisy but bounded."

## Test case

`SURF_test/run_all_airfoils.py`'s `CONV_LEVELS[0]` ("coarse") for SD8000,
unmodified: dx=0.04 (nx=150, ny=75), domain length=6/xoffset=-2/
yoffset=-1.5, ngrid=1, Re=60800, &alpha;=-0.81&deg; (`conv_alpha`), dt=0.01,
nsteps=3000 (t=0 to 30). `restart=250` added (the original conv_coarse run
used `restart=0`, no snapshots) to get 13 vorticity snapshots per run for
the figures below, same convention as `../README.md`'s NACA0012 check.

Run via [`run_static_suite_sd8000.py`](run_static_suite_sd8000.py):
`py_static/ibpm.py` 5 times in sequence &rarr; `py_run1/` .. `py_run5/`,
then `build_static/ibpm` 5 times in sequence &rarr; `cpp_run1/` ..
`cpp_run5/`. Wall time: ~16s/run (py), ~12s/run (cpp) — this coarse grid
(nx=150 vs. NACA0012's nx=300) is ~4x fewer points, so faster despite the
same nsteps.

## Results

### Reproducibility: still perfect, even through a chaotic blow-up

Every `.bin` restart snapshot and the `.force` trace were compared
byte-for-byte across all 5 runs of each implementation:

- **py_run2 .. py_run5 are byte-identical to py_run1** — all 13 snapshots,
  full `.force` trace, every run.
- **cpp_run2 .. cpp_run5 are byte-identical to cpp_run1** — same.
- **py_static blows up to NaN at exactly step 2827 (t=28.27) in every
  single one of its 5 runs** — not just "close," the *identical* step,
  confirmed via `.force` (all 5 traces first go `nan` at line 2827).
  `cpp_static` stays bounded (no NaN, max|force| ~2.46) through t=30 in
  all 5 of its runs.

This is the headline result: previously (`7-chaos_sensitivity/`), which
implementation blew up first — and when — was itself irreproducible,
because `FFTW_EXHAUSTIVE`'s timed plan search was a source of run-to-run
last-bit perturbation feeding directly into this chaotic instability. With
both implementations now pinned to `FFTW_ESTIMATE` (no timed search, same
codelet chosen every time), that source of irreproducibility is gone:
**the blow-up itself is now a fully deterministic, exactly reproducible
event**, not a chaotic coin-flip. See `reproducibility_diff.png`.

### Python vs. C++: chaotic amplification, exactly as previously documented

Unlike NACA0012 (roundoff-level agreement throughout), here the two
implementations' trajectories genuinely diverge — this *is* the expected
behavior for a chaotic system, not a fidelity gap (`py_vs_cpp_field_diff.txt`):

| t | max&#124;&Delta;&omega;&#124; | peak&#124;&omega;&#124; (cpp) | note |
|---|---|---|---|
| 0.0 | 0.0 | 0.0 | |
| 2.5 | 2.5e-12 | 96.6 | roundoff floor |
| 5.0 | 3.3e-10 | 91.0 | |
| 7.5 | 3.6e-7 | 144.5 | |
| 10.0 | 3.9e-4 | 107.4 | |
| 12.5 | 3.1e-2 | 113.2 | |
| 15.0 | 1.1e+2 | 116.5 | difference now same order as the field itself |
| 20.0 | 1.98e+2 | 135.1 | |
| 25.0 | 3.57e+2 | 160.4 | |
| 27.5 | 5.28e+2 | 341.9 | |
| 30.0 | n/a | 447.2 | py is NaN past its blow-up at t=28.27 |

The difference grows from the floating-point floor (~1e-12 at t=2.5) to
the same order as the vorticity field itself (~100+) by t=15 — textbook
exponential amplification of last-bit roundoff by a chaotic system,
**exactly the mechanism `SURF_test/gen_port_fidelity_diagnostic.py`'s
Panel C already documented** (bit-identical early, exponential divergence
from the roundoff floor). See `py_vs_cpp_diff.png`.

## Anomalies

**None beyond the case's own known, pre-documented instability.** py_static
diverging from cpp_static is expected and previously documented — not new.
What *is* new and worth flagging as the actual finding here: the blow-up
step (2827) and the full divergence trajectory are now **exactly
reproducible**, run after run, in a case that was previously shown to be
sensitive enough to reshuffle which implementation blew up first. No
unexpected behavior — every result above is consistent with what fixing
the DST algorithm should do and nothing more.

## Figures

- `flow_evolution_py_vs_cpp.png` — 2-row vorticity-field evolution,
  `py_run1` (top) vs. `cpp_run1` (bottom), t=0 to 30. Both show broadband
  speckle from t&approx;2.5 onward (same "hardest" phenomenon documented
  throughout `SURF_test/airfoils/`); py's t=30 panel is a "blown up (NaN)"
  placeholder.
- `reproducibility_diff.png` — vorticity difference, run5 minus run1, one
  row per implementation, fixed &plusmn;1e-10 color axis. All panels are
  flat/zero, including py's NaN panel (labeled "both NaN, bit-identical" —
  confirmed at the byte level, not just via `isnan`).
- `py_vs_cpp_diff.png` — vorticity difference, `py_run1` minus `cpp_run1`,
  at each pre-blow-up snapshot, color axis rescaled per panel (see title
  for actual max&#124;&Delta;&omega;&#124;) since the magnitude spans
  ~1e-15 to ~5e2 across the run.

## Files

Each `py_run{1..5}/flowfield/` and `cpp_run{1..5}/flowfield/` contains
`flow.cmd`, `flow.force`, and 13 `flow?????.bin` snapshots. `.cholesky`
cache files are not committed (regenerable, same convention as `../`).
`run_timing_summary.txt` has per-run wall-clock times.
