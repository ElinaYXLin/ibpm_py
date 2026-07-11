"""
run_cylinder_re_sweep_py.py

Python (py/ibpm.py) counterpart of run_cylinder_re_sweep.py -- same
geometry/grid/Re values/dt/nsteps/restart cadence, so
gen_cylinder_re_sweep_figs.py can show a Python-vs-C++ comparison at
every Re, not just the C++-only result the first pass produced. This
closes the "make sure Python has perfect fidelity" gap: the earlier
Re-sweep/grid-refine experiments were deliberately C++-only (framed as a
resolution/Re *behavior* question, not a port-fidelity recheck) -- this
adds the Python side back so that claim is actually checked at every one
of the new Re/dx points, not just the pre-existing Re=100/61k baselines.

Usage: python3 SURF_test/vortall/run_cylinder_re_sweep_py.py [Re]
Output: SURF_test/vortall/2-Re_sweep/_run_data/Re<value>/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM = REPO / "examples" / "cylinder.geom"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTBASE = REPO / "SURF_test" / "vortall" / "2-Re_sweep" / "_run_data"

GRID = dict(nx=450, ny=200, ngrid=1, length=9, xoffset=-1, yoffset=-2)
RE_VALUES = [100, 500, 1000, 3000, 10000]
DT = 0.02
NSTEPS = 1500
DT_OVERRIDE = {1000: 0.005, 3000: 0.005, 10000: 0.005}
NSTEPS_OVERRIDE = {1000: 6000, 3000: 6000, 10000: 6000}


def run_one(Re, dt, nsteps, restart):
    outdir = OUTBASE / f"Re{int(Re)}"
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-u", str(RUNNER),
           "-geom", str(GEOM), "-name", "cyl", "-outdir", str(outdir),
           "-nx", str(GRID["nx"]), "-ny", str(GRID["ny"]), "-ngrid", str(GRID["ngrid"]),
           "-length", str(GRID["length"]), "-xoffset", str(GRID["xoffset"]),
           "-yoffset", str(GRID["yoffset"]), "-Re", str(Re), "-dt", str(dt),
           "-nsteps", str(nsteps), "-tecplot", "0", "-restart", str(restart), "-force", "1"]
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    return proc.returncode == 0, elapsed


def main():
    Re = float(sys.argv[1]) if len(sys.argv) > 1 else None
    values = [Re] if Re else RE_VALUES
    for Re in values:
        dt = DT_OVERRIDE.get(int(Re), DT)
        nsteps = NSTEPS_OVERRIDE.get(int(Re), NSTEPS)
        restart = round(nsteps / 6)
        ok, elapsed = run_one(Re, dt, nsteps, restart)
        print(f"Re={Re} (dt={dt}, nsteps={nsteps}): {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
