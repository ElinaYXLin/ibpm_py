"""
run_t146_sweep.py

Re-runs stage 1's steady (non-oscillating) NACA0012 angle-of-attack sweep
(`../1-paper_based/run_kurt_suite.py`, GRIDS["dx0.020"]["steady"]) with the
ONE change motivated by `old/README.md`: duration extended from t=30
(nsteps=3000, this repo's original convention) to t=146 (nsteps=14600),
matching Kurtulus's actual steady-case duration once her paper's t=100s is
converted to this solver's non-dimensional time via c/U_inf (see this
folder's README for the conversion and why it lands on 146, not 100).

Everything else is IDENTICAL to stage 1's fig1_mean_coefficients sweep:
domain 6c (x in [-2,4]), dx=0.02 (nx=300, ny=150), dt=0.01, Re=1000,
same 43-angle set (0-40 step 1, plus 50, 60), both py_static and cpp_static.
This isolates the effect of duration alone -- stage 4's original
run_faithful2.py changed duration AND domain AND grid at once for a single
angle; this sweep changes only duration, across all angles, so any shift
in the mean-coefficient curves can be attributed to duration specifically.

Resumable (a run whose flow.force already reaches nsteps is skipped).

Usage: python3 run_t146_sweep.py [njobs]
Output: runs/t146_sweep/steady_<impl>_a<NN>/
"""
import itertools
import pathlib
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
STAGE4 = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "4-single_case_faithful"
RUNS = STAGE4 / "runs" / "t146_sweep"
CPP_BIN = REPO / "build_static" / "ibpm"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"
GEOM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"

DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 1000
NX, NY = 300, 150
DT = 0.01
NSTEPS = 14600  # t=146
ANGLES = list(range(0, 41)) + [50, 60]  # same 43-angle set as stage 1
IMPLS = ["py", "cpp"]
SNAP_ANGLES = {0, 9, 12}  # keep vorticity restarts for the same angles stage 1 kept


def is_done(outdir, nsteps):
    f = outdir / "flow.force"
    if not f.exists():
        return False
    try:
        last = None
        with open(f) as fh:
            for line in fh:
                if line.strip():
                    last = line
        return last is not None and int(last.split()[0]) >= nsteps
    except Exception:
        return False


def run_one(impl, angle):
    outdir = RUNS / f"steady_{impl}_a{angle:02d}"
    if is_done(outdir, NSTEPS):
        return (impl, angle, "skip", 0.0)
    outdir.mkdir(parents=True, exist_ok=True)
    restart = NSTEPS // 12 if angle in SNAP_ANGLES else 0
    common = [
        "-geom", str(GEOM), "-name", "flow", "-outdir", str(outdir),
        "-nx", str(NX), "-ny", str(NY), "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(angle), "-Re", str(RE),
        "-dt", str(DT), "-nsteps", str(NSTEPS),
        "-tecplot", "0", "-restart", str(restart), "-force", "1",
    ]
    cmd = ([str(CPP_BIN)] if impl == "cpp" else [sys.executable, "-u", str(PY_RUNNER)]) + common
    t0 = time.time()
    with open(outdir / "run_log.txt", "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    status = "ok" if (proc.returncode == 0 and is_done(outdir, NSTEPS)) else "FAIL"
    return (impl, angle, status, dt)


def main():
    njobs = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    jobs = list(itertools.product(IMPLS, ANGLES))
    print(f"{len(jobs)} jobs, {njobs}-way parallel", flush=True)

    done = fail = skip = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=njobs) as ex:
        futs = [ex.submit(run_one, *j) for j in jobs]
        for fut in as_completed(futs):
            impl, angle, status, dt = fut.result()
            if status == "skip":
                skip += 1
            elif status == "ok":
                done += 1
            else:
                fail += 1
                print(f"  FAIL: {impl} a{angle}", flush=True)
            n = done + fail + skip
            if n % 10 == 0 or status == "FAIL":
                print(f"[{n}/{len(jobs)}] done={done} skip={skip} fail={fail} "
                      f"elapsed={time.time()-t0:.0f}s", flush=True)
    print(f"DONE. ran={done} skipped={skip} failed={fail} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
