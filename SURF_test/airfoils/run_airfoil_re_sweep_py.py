"""
run_airfoil_re_sweep_py.py

Python (py/ibpm.py) counterpart of run_airfoil_re_sweep.py -- same
geometry/grid/alpha/Re values/dt/nsteps, so gen_airfoil_re_sweep_figs.py
can show Python-vs-C++ agreement at every swept Re, closing the
port-fidelity gap the original (C++-only) sweep left open.

Usage: python3 SURF_test/airfoils/run_airfoil_re_sweep_py.py <SD7003|SD8000> [Re]
Output: SURF_test/airfoils/<name>/4-Re_sweep/_run_data/Re<value>/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

CASES = {
    "SD7003": dict(geom=GEOMDIR / "sd7003_dx0.0200.geom", alpha=4.60),
    "SD8000": dict(geom=GEOMDIR / "sd8000_dx0.0200.geom", alpha=5.36),
}

RE_VALUES = [200, 500, 1000, 5000, 10000, 20000, 40000]
DT = 0.01
NSTEPS = 3000
RESTART = 500


def run_one(name, cfg, Re):
    outdir = REPO / "SURF_test" / "airfoils" / name / "4-Re_sweep" / "_run_data" / f"Re{int(Re)}"
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-u", str(RUNNER),
           "-geom", str(cfg["geom"]), "-name", "run", "-outdir", str(outdir),
           "-nx", "300", "-ny", "150", "-ngrid", "1",
           "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
           "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(cfg["alpha"]), "-Re", str(Re),
           "-dt", str(DT), "-nsteps", str(NSTEPS), "-tecplot", "0",
           "-restart", str(RESTART), "-force", "1"]
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    return proc.returncode == 0, elapsed


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in CASES:
        print("Usage: run_airfoil_re_sweep_py.py <SD7003|SD8000> [Re]", file=sys.stderr)
        sys.exit(1)
    name = sys.argv[1]
    cfg = CASES[name]
    values = [float(sys.argv[2])] if len(sys.argv) > 2 else RE_VALUES
    for Re in values:
        ok, elapsed = run_one(name, cfg, Re)
        print(f"{name} Re={Re}: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
