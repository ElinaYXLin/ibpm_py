"""
run_airfoil_re_sweep.py

Pushes SD7003/SD8000's Reynolds number DOWN from their usual Re~60-61k
(SURF_test/airfoils/LSAT-{SD7003,SD8000}/2-c++included/) to see where -- if
anywhere -- the broadband vorticity speckle documented there clears up.
Companion to SURF_test/vortall/run_cylinder_re_sweep.py, which pushes the
opposite direction (UP from the cylinder's clean Re=100 baseline) --
together the two sweeps converge on the same question from either side.

Same geometry/grid/alpha as the existing flowfield case (dx=0.02,
alpha=4.60 for SD7003 / 5.36 for SD8000) but C++ build/ibpm only, same
rationale as the cylinder sweep: this is a resolution/Re *behavior*
question, not a re-check of port fidelity (already established at the
baseline Re elsewhere in this suite).

Usage: python3 SURF_test/airfoils/run_airfoil_re_sweep.py <SD7003|SD8000>
Output: SURF_test/airfoils/<name>/4-Re_sweep/_run_data_cpp/Re<value>/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
CPP_BIN = REPO / "build" / "ibpm"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

CASES = {
    "SD7003": dict(geom=GEOMDIR / "sd7003_dx0.0200.geom", alpha=4.60, Re_baseline=61100),
    "SD8000": dict(geom=GEOMDIR / "sd8000_dx0.0200.geom", alpha=5.36, Re_baseline=60800),
}

RE_VALUES = [200, 500, 1000, 5000, 10000, 20000, 40000]
DT = 0.01
NSTEPS = 3000  # t = 30
RESTART = 500  # 7 snapshots: t=0,5,...,30


def run_one(name, cfg, Re):
    outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "4-Re_sweep" / "_run_data_cpp" / f"Re{int(Re)}"
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(CPP_BIN),
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
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found. Build it first: cd build && make", file=sys.stderr)
        sys.exit(1)
    if len(sys.argv) < 2 or sys.argv[1] not in CASES:
        print("Usage: run_airfoil_re_sweep.py <SD7003|SD8000> [Re]", file=sys.stderr)
        sys.exit(1)
    name = sys.argv[1]
    cfg = CASES[name]
    values = [float(sys.argv[2])] if len(sys.argv) > 2 else RE_VALUES
    for Re in values:
        ok, elapsed = run_one(name, cfg, Re)
        print(f"{name} Re={Re}: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
