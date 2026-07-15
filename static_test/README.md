# Reproducibility check: fixed-algorithm DST (py_static / cpp_static)

Confirms that pinning the FFTW discrete-sine-transform planner flag to a
fixed, non-adaptive choice (`FFTW_EXHAUSTIVE` &rarr; `FFTW_ESTIMATE |
FFTW_UNALIGNED`, in [`../py_static/_fftw_dst_shim.c`](../py_static/_fftw_dst_shim.c)
and [`../cpp_static/EllipticSolver2d.cc`](../cpp_static/EllipticSolver2d.cc))
makes each implementation's own output bit-identical run-to-run, and checks
how closely `py_static` and `cpp_static` agree with each other.

## What was confirmed unchanged first

`diff -rq` of `py_static/` vs `py/` and `cpp_static/` vs `src/` shows no
differences outside the two documented DST-flag lines plus one required
supporting change (`py_static/_fftw_native.py`'s cached-library filename,
so `py_static` doesn't silently reuse `py/`'s already-built
`FFTW_EXHAUSTIVE` shared object from `build/`).

## Test case

The standard NACA0012 flow-evolution case from
[`../SURF_test/low_re/run_naca0012.py`](../SURF_test/low_re/run_naca0012.py)
(default parameters, unmodified): dx=0.02 (nx=300, ny=150), domain
length=6/xoffset=-2/yoffset=-1.5, ngrid=1, Re=500, alpha=5&deg;, dt=0.01,
nsteps=3000 (t=0 to 30), restart every 500 steps &rarr; 7 vorticity
snapshots per run (t=0,5,10,...,30).

Run via [`run_static_suite.py`](run_static_suite.py):
[`run_ibpm_case_static.py`](run_ibpm_case_static.py) (imports `py_static`
instead of `py`) 5 times in sequence &rarr; `py_run1/` .. `py_run5/`, then
`build_static/ibpm` (a from-scratch build of `cpp_static/`, via
[`../build_static/Makefile`](../build_static/Makefile)) 5 times in sequence
&rarr; `cpp_run1/` .. `cpp_run5/`. Wall time: ~51s/run (py), ~46s/run
(cpp), ~8 minutes total for all 10 (see `run_timing_summary.txt`) — much
faster than the original `FFTW_EXHAUSTIVE` runs would have been, since
`FFTW_ESTIMATE` skips the timed codelet search entirely.

## Results

**Reproducibility (the point of this check): perfect, no anomalies.**
Every `.bin` restart snapshot and the `.force` time-history file were
compared byte-for-byte (`cmp`) across all 5 runs of each implementation:

- **py_run2 .. py_run5 are byte-identical to py_run1** — all 7 snapshots, every run.
- **cpp_run2 .. cpp_run5 are byte-identical to cpp_run1** — all 7 snapshots, every run.
- Max absolute vorticity-field difference between any two same-implementation
  runs: **exactly 0.0** (not just "small" — literally zero, confirmed with
  `numpy.max(numpy.abs(...))`).

This is the expected result now that both implementations plan their DST
with `FFTW_ESTIMATE` (a fixed heuristic, no timed search) instead of
`FFTW_EXHAUSTIVE` (which searches and can pick different codelets run to
run, machine-load permitting). See `reproducibility_diff.png`.

**Python vs. C++: agree at floating-point-roundoff level, not bit-identical.**
`py_vs_cpp_field_diff.txt`:

| t | max&#124;&Delta;&omega;&#124; | rms&#124;&Delta;&omega;&#124; | peak&#124;&omega;&#124; (cpp) |
|---|---|---|---|
| 0.0 | 0.0 | 0.0 | 0.0 |
| 5.0 | 7.99e-13 | 1.66e-14 | 68.82 |
| 10.0 | 7.64e-13 | 1.65e-14 | 68.61 |
| 15.0 | 6.13e-13 | 1.56e-14 | 68.61 |
| 20.0 | 5.44e-13 | 1.51e-14 | 68.61 |
| 25.0 | 1.38e-12 | 2.17e-14 | 68.61 |
| 30.0 | 1.24e-12 | 2.04e-14 | 68.61 |

`.force` (integrated Cl/Cd over all 3000 steps) is byte-identical between
`py_run1` and `cpp_run1`. The field-level difference (max ~1e-12 against a
peak vorticity of ~68.6, i.e. ~14-15 significant digits of agreement) is
consistent with ordinary floating-point roundoff accumulating differently
across two independently-compiled binaries (Python/ctypes/numpy vs. g++),
not a fidelity gap — see `py_vs_cpp_diff.png`.

## Anomalies

**None found.** Both implementations reproduced exactly across all 5 runs;
Python and C++ track each other to roundoff precision throughout the run,
with no growth pattern suggesting divergence. Nothing here required
investigation beyond what's reported above.

## Figures

- `flow_evolution_py_vs_cpp.png` — standard 2-row vorticity-field evolution,
  `py_run1` (top) vs. `cpp_run1` (bottom), t=0 to 30.
- `reproducibility_diff.png` — vorticity difference, run5 minus run1, one
  row per implementation. Fixed color axis of &plusmn;1e-10 to give the
  panels a visible scale; they render flat because the true difference is
  exactly zero everywhere (see reproducibility results above).
- `py_vs_cpp_diff.png` — vorticity difference, `py_run1` minus `cpp_run1`,
  at each snapshot, on a roundoff-scale (&plusmn;2e-12) color axis, showing
  the faint structure of the two implementations' floating-point-level
  disagreement.

## Files

Each `py_run{1..5}/flowfield/` and `cpp_run{1..5}/flowfield/` contains
`flow.cmd` (exact command line), `flow.force` (per-step Cl/Cd), and the 7
`flow?????.bin` vorticity/flux restart snapshots. `.cholesky` cache files
(solver setup, regenerable, ~1.4MB/run) are not committed, matching this
repo's existing convention for other `flowfield/` directories.
`run_timing_summary.txt` has per-run wall-clock times.
