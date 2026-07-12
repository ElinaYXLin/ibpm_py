"""
run_cylinder_re_sweep.py

Pushes the cylinder case's Reynolds number UP from 1-baseline/'s Re=100
(clean, matches VORTALL.mat) to see where -- if anywhere -- the broadband
vorticity speckle seen throughout SURF_test/airfoils/ (Re~40-61k) first
appears. Same geometry/grid as 1-baseline/ (examples/cylinder.geom,
nx=450, ny=200, dx=0.02, domain length=9/xoffset=-1/yoffset=-2) but a much
shorter transient (t=30 instead of t=280 -- this is a resolution/Re
*behavior* study, not a statistically-converged shedding average) and
C++ build/ibpm only (this is the same, unmodified upstream solver
1-baseline/ already validated Python against at Re=100 -- see
1-baseline/vorticity_comparison_3way.png -- so re-validating port fidelity
at every new Re here would be redundant; the open question is solver
*behavior* vs. resolution/Re, not port correctness).

Usage: python3 SURF_test/vortall/run_cylinder_re_sweep.py
Output: SURF_test/vortall/2-Re_sweep/_run_data_cpp/Re<value>/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM = REPO / "examples" / "cylinder.geom"
CPP_BIN = REPO / "build" / "ibpm"
OUTBASE = REPO / "SURF_test" / "vortall" / "2-Re_sweep" / "_run_data_cpp"

GRID = dict(nx=450, ny=200, ngrid=1, length=9, xoffset=-1, yoffset=-2)
RE_VALUES = [100, 500, 1000, 3000, 10000]
DT = 0.02
NSTEPS = 1500  # t = 30
RESTART = 250  # 7 snapshots: t=0,5,...,30

# Re>=1000 diverged within the first ~30 steps at dt=0.02 -- the
# impulsive-start transient (sharp initial vorticity gradient at the
# boundary) is too violent for this dt once viscous damping drops enough;
# not a late-time chaotic blowup like the ngrid=3 case, a genuine CFL
# violation from step 1. Fixed with a 4x smaller dt, same physical time.
DT_OVERRIDE = {1000: 0.005, 3000: 0.005, 10000: 0.005}
NSTEPS_OVERRIDE = {1000: 6000, 3000: 6000, 10000: 6000}


def run_one(Re, dt, nsteps, restart):
    outdir = OUTBASE / f"Re{int(Re)}"
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(CPP_BIN),
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
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found. Build it first: cd build && make", file=sys.stderr)
        sys.exit(1)

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
