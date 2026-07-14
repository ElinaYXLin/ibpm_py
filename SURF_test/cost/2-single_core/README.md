# Single-core cost comparison: forcing Python onto one CPU core

This directory re-runs the exact same benchmark as
[`../1-multi-core/`](../1-multi-core/) (same cylinder/Re=100/RK3 case,
same 4 grid resolutions, same phase-timing/`getrusage`/`psutil`
methodology, reused directly by importing `../1-multi-core/run_benchmark.py`
rather than reimplementing it), but forces every subprocess onto a single
CPU core's worth of parallelism, so the Python-vs-C++ comparison holds a
fixed "one core" assumption for both implementations, per a mentor request.

**No file under `py/` was modified to do this.** See "Why only Python
needed anything at all" below for why that's possible.

## What was actually done

[`run_benchmark_single_core.py`](run_benchmark_single_core.py) sets five
environment variables before launching each subprocess:

```python
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"
```

`subprocess.Popen` (called inside the imported, unmodified
`run_one()`) inherits the parent process's environment by default, so
setting these here is sufficient -- no code in `py/`, `src/`,
`backends.py`, or `run_benchmark.py` itself needed to change. This is the
standard, well-established numpy/scipy-ecosystem mechanism for
reproducible single-core benchmarking.

## Why only Python needed anything at all

Confirmed directly by inspection before running anything:

- **C++ (`build/ibpm`) never spawns a worker thread.** `grep -rn "#pragma omp" src/` -> nothing.
  No `pthread`, no `fftw_init_threads`/`fftw_plan_with_nthreads` calls anywhere in `src/`.
  `build/Makefile` links `-lfftw3` only (not `-lfftw3_threads`), and
  `nm build/ibpm | grep fftw.*thread` finds nothing. C++ was already an
  inherent single-core baseline; the env vars above don't change its
  behavior (and running it with them set anyway, for a same-conditions
  comparison, was harmless and confirmed to change nothing -- see below).
- **`py/ibpm.py`'s own FFTW usage is equally single-threaded** by the same
  argument: the custom shim (`py/_fftw_native.py` + the compiled
  `_fftw_dst_shim.so`) issues plain `fftw_plan_r2r_2d`/`fftw_execute`
  calls via ctypes, no thread-init calls, matching C++.
- **The one real risk was underneath, in numpy's BLAS backend, not in this
  repo's code at all.** `np.show_config()` in this environment shows numpy
  linked against `scipy-openblas` with `MAX_THREADS=64`. Every `np.dot`
  call in `py/cholesky_solver.py` (`computeFactorization`/`Minv`),
  `py/conjugate_gradient_solver.py`, `InnerProduct`, etc. dispatches into
  OpenBLAS, which decides on its own, per call, whether a given array size
  is worth spreading across multiple threads -- exactly the kind of hidden
  parallelism a single-core comparison needs to rule out explicitly rather
  than assume away.

## Result: it was real, and large -- at peak, not on average

The clearest single number is **peak instantaneous CPU%** from the
`psutil` time series (each backend's own worst moment during the run, not
an average):

| backend | nx=ny | multi-core peak CPU% | single-core peak CPU% |
|---|---|---|---|
| cpp    | 100 | 112.2 | 111.4 |
| cpp    | 200 | 111.7 | 114.0 |
| cpp    | 300 | 112.5 | 115.2 |
| cpp    | 400 | 114.7 | 115.5 |
| python | 100 | **722.3** | 110.9 |
| python | 200 | **610.5** | 113.6 |
| python | 300 | **726.0** | 112.8 |
| python | 400 | **787.1** | 114.8 |

C++ never exceeds ~115% at any grid size, in either condition (100% = one
core saturated; the extra ~10-15% is ordinary OS/measurement noise from
other threads/processes on the machine, not real parallelism). **Python,
without the thread-limiting env vars, spiked to 610-787% -- using 6-8
cores simultaneously at its peak moment** (almost certainly during the
Cholesky `computeMatrixM`/`computeFactorization` phase, the part of the
code that actually calls `np.dot`). With the env vars set, Python's peak
drops to 111-115%, matching C++ almost exactly. This confirms the
suspected mechanism precisely: OpenBLAS was silently parallelizing across
most of this machine's cores whenever Python ran, undetectable from
wall-clock time alone.

## But wall-clock time barely changed -- the extra cores weren't buying real speed

| nx=ny | wall time ratio (multi-core) | wall time ratio (single-core) | CPU-seconds ratio (multi-core) | CPU-seconds ratio (single-core) |
|---|---|---|---|---|
| 100 | 1.17x | 1.14x | 1.49x | 1.09x |
| 200 | 1.06x | 1.06x | 1.17x | 1.05x |
| 300 | 1.03x | 1.04x | 1.08x | 1.04x |
| 400 | 1.01x | 0.99x | 1.05x | 0.99x |

(ratio = Python / C++; see `tables/cost_ratio_python_vs_cpp.md` in each
directory for the full table including peak-RSS and timestepping-only
ratios, which are likewise essentially unchanged between conditions.)

Despite using up to 8 cores' worth of instantaneous parallelism in the
multi-core run, **Python's wall-clock time relative to C++ moved by at
most 0.03x** across every grid size -- nowhere close to the 6-8x
speedup you'd expect if that parallelism were doing useful work. This is
the classic small-problem-size threading-overhead story: the boundary-point
counts here (79-314) make the Cholesky matrices modest (`size = 2 x
numPoints`, up to ~628), and OpenBLAS's decision to spread a
several-hundred-element dot product across 6-8 threads costs more in
thread-spawn/synchronization overhead than it saves in compute --
multi-core mode was burning extra CPU-seconds (up to 1.49x C++'s at
nx=100, vs. 1.09x once forced single-core) without any corresponding
wall-clock benefit. **The single-core measurement is therefore not just
"more comparable" -- it's a strictly more honest per-core cost number**,
since the multi-core CPU-seconds figure was partly measuring wasted
parallelization overhead rather than Python's actual algorithmic cost per
core.

Peak RAM is essentially unchanged between conditions (both ~148.5MB at
nx=100-300 python, ~140.9-148.5MB range -- a few MB of noise, not a real
difference), confirming the earlier `../1-multi-core/README.md` RAM
finding is independent of threading.

## Files

| File | What it is |
|---|---|
| `run_benchmark_single_core.py` | Sets the thread-limiting env vars, then reuses `../1-multi-core/run_benchmark.py`'s `BACKENDS`/`run_one`/`ResourceMonitor`/grid sweep by direct import -- not a reimplementation. |
| `gen_cost_report.py` | Identical copy of `../1-multi-core/gen_cost_report.py` (it's already fully self-relative via `pathlib.Path(__file__).resolve().parent`, so it needed zero changes to work here -- verified). |
| `raw/cost_results.json`, `raw/<backend>_nx<N>/` | Raw per-run data for this (single-core) condition, same format as `../1-multi-core/raw/`. |
| `tables/`, `*.png` | Same tables/figures as `../1-multi-core/`, regenerated from this condition's data. |

## Regenerating

```bash
cd build && make          # if not already built
python3 SURF_test/cost/2-single_core/run_benchmark_single_core.py
python3 SURF_test/cost/2-single_core/gen_cost_report.py
```
