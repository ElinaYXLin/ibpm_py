# run_benchmark_single_core.py
#
# Same benchmark as ../1-multi-core/run_benchmark.py (same case, same grid
# sweep, same phase-timing/getrusage/psutil methodology -- reused directly
# by import, not reimplemented), but forces every subprocess it launches to
# run on a single CPU core's worth of parallelism, so the comparison holds
# a fixed "one core" assumption for both backends. No file under py/ is
# touched to achieve this -- see the module docstring below for why that's
# possible.
#
# Why only Python needed this in the first place (see ../README.md for the
# full writeup): build/ibpm links only -lfftw3 (not -lfftw3_threads), has no
# OpenMP pragmas and no fftw_plan_with_nthreads/fftw_init_threads calls
# anywhere in src/ -- confirmed by grepping the source and by `nm build/ibpm
# | grep fftw.*thread` (nothing) -- so C++ never spawns a worker thread and
# was already a single-core baseline. py/ibpm.py's own FFTW shim
# (py/_fftw_native.py + the compiled _fftw_dst_shim.so) is equally
# single-threaded by the same argument. The one place multi-core
# parallelism could sneak into the Python side is underneath numpy: this
# environment's numpy is linked against OpenBLAS ("scipy-openblas",
# MAX_THREADS=64 per `np.show_config()`), and every `np.dot` call in
# py/cholesky_solver.py / py/conjugate_gradient_solver.py / InnerProduct
# etc. dispatches into it -- OpenBLAS decides on its own, per-call, whether
# a given array size is worth spreading across multiple threads.
#
# Fix: OpenBLAS (and other common BLAS backends, kept here defensively in
# case this is ever run on a machine with a different one installed) reads
# its thread cap from environment variables at process-launch time. Setting
# them in *this* process, before any subprocess is spawned, is sufficient --
# subprocess.Popen (called inside the imported run_one()) inherits the
# parent's os.environ by default. This is the standard, well-established
# numpy/scipy-ecosystem mechanism for reproducible single-core benchmarking;
# no py/ code (or even backends.py/run_benchmark.py) needs to change.
#
# Usage: python3 SURF_test/cost/2-single_core/run_benchmark_single_core.py
# (run from the repository root; requires build/ibpm, same as the
# multi-core benchmark)

from __future__ import annotations

import os

# Must be set before importing/calling anything that might spawn a BLAS
# thread pool (numpy itself, or the subprocess this script launches).
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import json
import pathlib
import sys
import time

COST_DIR = pathlib.Path(__file__).resolve().parent
MULTI_CORE_DIR = COST_DIR.parent / "1-multi-core"
sys.path.insert(0, str(MULTI_CORE_DIR))
# Reuses BACKENDS, REPO_ROOT, make_circle_geom, GRID_SIZES, NSTEPS, dt_for,
# and run_one (phase-timing + getrusage + psutil ResourceMonitor) verbatim
# from the multi-core benchmark -- the methodology is identical; only the
# environment each subprocess launches into differs (set above).
import run_benchmark as mc  # noqa: E402

GEOM_DIR = MULTI_CORE_DIR / "geom"  # geometry doesn't depend on threading; reuse it
RAW_DIR = COST_DIR / "raw"
RAW_DIR.mkdir(exist_ok=True)


def main():
    results = []
    results_path = RAW_DIR / "cost_results.json"

    print(f"Thread-limiting env vars for this run: "
          f"{ {v: os.environ[v] for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')} }",
          flush=True)

    for backend in mc.BACKENDS:
        if not backend.implemented:
            print(f"=== skipping backend '{backend.name}': {backend.note} ===", flush=True)
            continue
        for nx in mc.GRID_SIZES:
            ny = nx
            dt = mc.dt_for(nx)
            dx = 4.0 / nx
            geom_path = GEOM_DIR / f"cylinder_dx{dx:.5f}.geom"
            npts = mc.make_circle_geom(dx, geom_path)

            name = f"{backend.name}_nx{nx}"
            outdir = RAW_DIR / name
            cmd = backend.build_cmd(geom_path, outdir, "run", nx, ny, dt, mc.NSTEPS)

            print(f"=== {backend.name}: nx=ny={nx} dx={dx:.4f} dt={dt} npts={npts} (single-core) ===",
                  flush=True)
            t0 = time.time()
            r = mc.run_one(cmd, outdir)
            elapsed = time.time() - t0
            r.update(backend=backend.name, label=backend.label, nx=nx, ny=ny, dx=dx, dt=dt,
                     nsteps=mc.NSTEPS, npts=npts)
            results.append(r)
            print(f"    wall={r['wall_time']:.2f}s  (model={r['phase_model']:.2f}s "
                  f"setup={r['phase_setup']:.2f}s steps={r['phase_steps']:.2f}s)  "
                  f"cpu={r['cpu_total']:.2f}s  peakRSS={r['peak_rss_bytes']/1e6:.1f}MB  "
                  f"[{elapsed:.1f}s elapsed]", flush=True)

            results_path.write_text(json.dumps(results, indent=2))

    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
