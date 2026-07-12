"""
run_small_dt_validation.py

Runs py/ibpm.py (now with the real-FFTW3, FFTW_EXHAUSTIVE-planned sine
transform -- see py/elliptic_solver_2d.py / py/_fftw_native.py) and this
repo's C++ reference build (build/ibpm, unchanged) on SD7003 and SD8000,
at the SAME grid as the production polar sweep (dx=0.02, nx=300, ny=150)
but at a MUCH SMALLER timestep than the rest of this test suite uses
(dt=0.001 vs. the usual dt=0.01 -- 10x smaller), with restart snapshots
written periodically so intermediate field state (not just the final
force trace) can be directly compared between the two implementations.

This exists to validate the FFTW-authenticity change in
py/elliptic_solver_2d.py: does the Python port, now calling the real FFTW3
library with FFTW_EXHAUSTIVE (the same planner, the same search over
candidate sine-transform algorithms, as src/EllipticSolver2d.cc), still
agree with the C++ reference to the same bit-level fidelity established
elsewhere in this test suite (see SD8000/2-c++included/
port_fidelity_diagnostic.png)? A small dt reduces time-discretization
error, which sharpens that comparison (less truncation-error "noise" to
distinguish from genuine FFTW-backend differences).

Usage:  python3 SURF_test/run_small_dt_validation.py
Output: SURF_test/airfoils/LSAT-{SD7003,SD8000}/_run_data_smalldt/       (Python)
        SURF_test/airfoils/LSAT-{SD7003,SD8000}/_run_data_smalldt_cpp/   (C++)
"""
import json
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
CPP_BIN = REPO / "build" / "ibpm"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

DX = 0.02
NX, NY = 300, 150
DT = 0.001          # 10x smaller than this suite's usual dt=0.01
NSTEPS = 2000        # t = 2.0
RESTART_EVERY = 200   # 11 snapshots (steps 0,200,...,2000) for intermediate-value comparison

CASES = {
    "SD7003": dict(Re=61100, alpha=-0.09),   # same alpha as the grid-convergence study
    "SD8000": dict(Re=60800, alpha=-0.81),
}

RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
RUNNER.write_text(
    "import sys, types, pathlib\n"
    "repo_root = pathlib.Path(%r)\n" % str(REPO) +
    "sys.path.insert(0, str(repo_root))\n"
    "pkg = types.ModuleType('py')\n"
    "pkg.__path__ = [str(repo_root / 'py')]\n"
    "sys.modules['py'] = pkg\n"
    "from py.ibpm import main\n"
    "sys.exit(main(['py.ibpm'] + sys.argv[1:]))\n"
)


def _grid_args(geom, name, outdir, alpha, Re):
    return [
        "-geom", str(geom), "-name", name, "-outdir", str(outdir),
        "-nx", str(NX), "-ny", str(NY), "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(alpha), "-Re", str(Re),
        "-dt", str(DT), "-nsteps", str(NSTEPS), "-tecplot", "0",
        "-restart", str(RESTART_EVERY), "-force", "1",
    ]


def run_python(geom, name, outdir, alpha, Re):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-u", str(RUNNER)] + _grid_args(geom, name, outdir, alpha, Re)
    log_path = outdir / f"{name}_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"python run failed ({name}): see {log_path}")
    return elapsed


def run_cpp(geom, name, outdir, alpha, Re):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(CPP_BIN)] + _grid_args(geom, name, outdir, alpha, Re)
    log_path = outdir / f"{name}_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"C++ run failed ({name}): see {log_path}")
    return elapsed


def main():
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found. Build it first: cd build && make", file=sys.stderr)
        sys.exit(1)

    timings = {}
    for name, cfg in CASES.items():
        geom = GEOMDIR / f"{name.lower()}_dx{DX:.4f}.geom"
        if not geom.exists():
            raise FileNotFoundError(f"{geom} missing -- run SURF_test/run_all_airfoils.py first")

        py_outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt"
        cpp_outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt_cpp"

        print(f"=== {name}: Python (native FFTW3, FFTW_EXHAUSTIVE), dt={DT}, nsteps={NSTEPS} ===",
              flush=True)
        t_py = run_python(geom, "run", py_outdir, cfg["alpha"], cfg["Re"])
        print(f"  done in {t_py:.1f}s", flush=True)
        timings.setdefault(name, {})["python_s"] = t_py

        print(f"=== {name}: C++ build/ibpm, dt={DT}, nsteps={NSTEPS} ===", flush=True)
        t_cpp = run_cpp(geom, "run", cpp_outdir, cfg["alpha"], cfg["Re"])
        print(f"  done in {t_cpp:.1f}s", flush=True)
        timings[name]["cpp_s"] = t_cpp

    meta = dict(dx=DX, nx=NX, ny=NY, dt=DT, nsteps=NSTEPS, restart_every=RESTART_EVERY,
                cases=CASES, timings=timings)
    (REPO / "SURF_test" / "small_dt_run_meta.json").write_text(json.dumps(meta, indent=2))
    print("ALL DONE")


if __name__ == "__main__":
    main()
