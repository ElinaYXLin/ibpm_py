"""
run_kurt_suite.py

Reproduces the test cases of Kurtulus (2019), "Unsteady aerodynamics of a
pitching NACA 0012 airfoil at low Reynolds number" (Int J Micro Air Veh 11:1-21,
DOI 10.1177/1756829319890609) with this repo's immersed-boundary solver, using
BOTH py_static (py_static/ibpm.py) and cpp_static (build_static/ibpm).

The paper is a NUMERICAL (ANSYS Fluent) study, not experimental. It runs a
NACA0012 at Re=1000: a non-oscillatory ("steady") baseline plus a sinusoidally
PITCHING airfoil (+/-1 deg about the quarter chord) at f=1 Hz and f=4 Hz, for
mean angles of attack 0-60 deg.

Non-dimensionalization (solver uses chord=1, U_inf=1):
  - Re = 1000.
  - Pitch amplitude A = 1 deg = 0.0174533 rad.
  - Pitch frequency (non-dim, cycles per convective time) = f * c/U_inf with the
    paper's c=0.1 m, U_inf=0.146 m/s: 1 Hz -> 0.684932, 4 Hz -> 2.739726. These
    give reduced frequencies k = 2*pi*f*c/U_inf = 4.30 and 17.2, matching the
    paper.
  - Pitch phase = pi so the effective instantaneous AoA is alpha0 + A*sin(2 pi f t)
    (matching the paper); the mean AoA alpha0 is applied via the -alpha flag
    (free-stream tilt), the pitch via `motion pitchplunge`.

Sweep (user-requested "full" set):
  - Angles: 0..40 deg step 1, plus 50, 60  (the paper's own increment).
  - Motions: steady, f1hz, f4hz.
  - Implementations: py_static, cpp_static.
  - Grid: dx=0.02 (full sweep) and dx=0.01 (KEY angles only -- grid check).

Runs are executed with a bounded process pool and are resumable (a run whose
.force file already ends at the target step is skipped).

Usage:
  python3 SURF_test/low_re/NACA0012/kurt_comp/run_kurt_suite.py [dx002|dx001|all] [njobs]
Output: kurt_comp/runs/<grid>/<motion>_<impl>_a<NN>/
"""
import itertools
import pathlib
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "1-paper_based"
GEOMDIR = KURT / "geom"
RUNS = KURT / "runs"
CPP_BIN = REPO / "build_static" / "ibpm"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"

DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 1000
PITCH_AMP = 0.0174532925199433  # 1 deg in radians
PITCH_PHASE = 3.14159265358979  # pi -> effective AoA = alpha0 + A sin(2 pi f t)
FREQ = {"f1hz": 0.684931506849315, "f4hz": 2.73972602739726}  # f*c/U_inf, c=0.1 U=0.146

# angle set: 0..40 step 1, then 50, 60
ANGLES_FULL = list(range(0, 41)) + [50, 60]
ANGLES_KEY = [0, 5, 9, 15, 40]  # dx=0.01 grid-check subset

# snapshot (restart) angles -- keep vorticity fields for wake-contour figures
SNAP_ANGLES = {0, 9, 12}

GRIDS = {
    "dx0.020": dict(geom_dx="0.0200", nx=300, ny=150,
                    steady=dict(dt=0.01, nsteps=3000),
                    pitch=dict(dt=0.005, nsteps=6000),
                    angles=ANGLES_FULL),
    "dx0.010": dict(geom_dx="0.0100", nx=600, ny=300,
                    steady=dict(dt=0.005, nsteps=6000),
                    pitch=dict(dt=0.0025, nsteps=12000),
                    angles=ANGLES_KEY),
}
MOTIONS = ["steady", "f1hz", "f4hz"]
IMPLS = ["py", "cpp"]


def ensure_geom(grid_key, motion):
    """Return a .geom path for this grid+motion, creating pitching variants."""
    g = GRIDS[grid_key]
    raw = REPO / "SURF_test" / "geom" / f"naca0012_dx{g['geom_dx']}.txt"
    if motion == "steady":
        return REPO / "SURF_test" / "geom" / f"naca0012_dx{g['geom_dx']}.geom"
    GEOMDIR.mkdir(parents=True, exist_ok=True)
    gp = GEOMDIR / f"naca0012_dx{g['geom_dx']}_{motion}.geom"
    freq1 = FREQ[motion]
    gp.write_text(
        "body NACA0012\n"
        f"  raw {raw}\n"
        "  center 0.25 0.0\n"
        f"  motion pitchplunge {PITCH_AMP:.16g} {freq1:.16g} {PITCH_PHASE:.16g} 0 0 0\n"
        "end\n"
    )
    return gp


def target_steps(grid_key, motion):
    g = GRIDS[grid_key]
    return (g["steady"] if motion == "steady" else g["pitch"])["nsteps"]


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


def run_one(grid_key, motion, impl, angle):
    g = GRIDS[grid_key]
    params = g["steady"] if motion == "steady" else g["pitch"]
    geom = ensure_geom(grid_key, motion)
    outdir = RUNS / grid_key / f"{motion}_{impl}_a{angle:02d}"
    if is_done(outdir, params["nsteps"]):
        return (grid_key, motion, impl, angle, "skip", 0.0)
    outdir.mkdir(parents=True, exist_ok=True)
    restart = params["nsteps"] // 12 if angle in SNAP_ANGLES else 0
    common = [
        "-geom", str(geom), "-name", "flow", "-outdir", str(outdir),
        "-nx", str(g["nx"]), "-ny", str(g["ny"]), "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(angle), "-Re", str(RE),
        "-dt", str(params["dt"]), "-nsteps", str(params["nsteps"]),
        "-tecplot", "0", "-restart", str(restart), "-force", "1",
    ]
    if impl == "cpp":
        cmd = [str(CPP_BIN)] + common
    else:
        cmd = [sys.executable, "-u", str(PY_RUNNER)] + common
    t0 = time.time()
    with open(outdir / "run_log.txt", "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    status = "ok" if (proc.returncode == 0 and is_done(outdir, params["nsteps"])) else "FAIL"
    return (grid_key, motion, impl, angle, status, dt)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    njobs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    grid_keys = {"dx002": ["dx0.020"], "dx001": ["dx0.010"], "all": ["dx0.020", "dx0.010"]}[which]

    jobs = []
    for gk in grid_keys:
        for motion, impl, angle in itertools.product(MOTIONS, IMPLS, GRIDS[gk]["angles"]):
            jobs.append((gk, motion, impl, angle))
    print(f"{len(jobs)} jobs, {njobs}-way parallel", flush=True)

    done = fail = skip = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=njobs) as ex:
        futs = [ex.submit(run_one, *j) for j in jobs]
        for fut in as_completed(futs):
            gk, motion, impl, angle, status, dt = fut.result()
            if status == "skip":
                skip += 1
            elif status == "ok":
                done += 1
            else:
                fail += 1
                print(f"  FAIL: {gk} {motion} {impl} a{angle}", flush=True)
            n = done + fail + skip
            if n % 10 == 0 or status == "FAIL":
                print(f"[{n}/{len(jobs)}] done={done} skip={skip} fail={fail} "
                      f"elapsed={time.time()-t0:.0f}s", flush=True)
    print(f"DONE. ran={done} skipped={skip} failed={fail} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
