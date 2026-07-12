"""
run_cylinder_grid_refine.py

At a FIXED Reynolds number (picked from 2-Re_sweep/'s results as roughly
where the vorticity field stops being clean -- see 3-grid_refine/README.md
for the actual value used and why), refines dx across three levels
(0.04, 0.02, 0.01) to test whether resolution alone -- at fixed Re --
controls whether the field is clean or speckled. If refining dx cleans up
the SAME Re that looked speckled at dx=0.02, that's direct evidence for
delta/dx (boundary-layer thickness relative to grid spacing) as the
controlling parameter, not Re by itself.

Circle geometry generated fresh at each dx via the `circle_n` command
(examples/cylinder.geom's own syntax) with point count ~ circumference/dx,
so boundary point spacing tracks the grid spacing at every level (same
convention SD7003/SD8000 use via make_airfoil_raw.py).

Usage: python3 SURF_test/vortall/run_cylinder_grid_refine.py <Re>
Output: SURF_test/vortall/3-grid_refine/_run_data_cpp/dx<value>/
"""
import math
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
CPP_BIN = REPO / "build" / "ibpm"
OUTBASE = REPO / "SURF_test" / "vortall" / "3-grid_refine"
GEOMDIR = OUTBASE / "geom"

DOMAIN = dict(length=9, xoffset=-1, yoffset=-2)
# At Re=5000, the cylinder's blunt-body impulsive start diverges within
# the first ~15-18 steps at dt=0.02/dx=0.02 or dt=0.01/dx=0.01 (same
# early-blowup CFL signature as the Re-sweep's Re>=1000 cases, but the
# cylinder's stronger impulsive-start vorticity layer needs an even
# smaller dt than the equivalent-Re airfoil case did). Fixed with dt
# reduced further at every level, same physical time t=30.
DX_LEVELS = [
    dict(tag="coarse", dx=0.04, nx=225, ny=100, dt=0.02, nsteps=1500),
    dict(tag="medium", dx=0.02, nx=450, ny=200, dt=0.005, nsteps=6000),
    dict(tag="fine", dx=0.01, nx=900, ny=400, dt=0.0025, nsteps=12000),
]
RESTART_FRAC = 1 / 6  # 7 snapshots (t=0,5,...,30) regardless of nsteps


def make_geom(dx, tag):
    GEOMDIR.mkdir(parents=True, exist_ok=True)
    n = round(math.pi / dx)  # circumference = pi * diameter = pi * 1
    geom_path = GEOMDIR / f"cylinder_dx{dx:.4f}.geom"
    geom_path.write_text(f"# Cylinder, diameter 1, {n} points (ds~{math.pi/n:.4f}, dx={dx})\n\n"
                          f"body Cylinder\n  circle_n 0 0 0.5 {n}\nend\n")
    return geom_path, n


def run_one(Re, lvl):
    geom_path, n = make_geom(lvl["dx"], lvl["tag"])
    outdir = OUTBASE / "_run_data_cpp" / f"dx{lvl['dx']}"
    outdir.mkdir(parents=True, exist_ok=True)
    restart = round(lvl["nsteps"] * RESTART_FRAC)
    cmd = [str(CPP_BIN),
           "-geom", str(geom_path), "-name", "cyl", "-outdir", str(outdir),
           "-nx", str(lvl["nx"]), "-ny", str(lvl["ny"]), "-ngrid", "1",
           "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
           "-yoffset", str(DOMAIN["yoffset"]), "-Re", str(Re), "-dt", str(lvl["dt"]),
           "-nsteps", str(lvl["nsteps"]), "-tecplot", "0", "-restart", str(restart), "-force", "1"]
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    return proc.returncode == 0, elapsed, n


def main():
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr)
        sys.exit(1)
    Re = float(sys.argv[1]) if len(sys.argv) > 1 else 3000
    tag = sys.argv[2] if len(sys.argv) > 2 else None
    levels = [l for l in DX_LEVELS if tag is None or l["tag"] == tag]
    for lvl in levels:
        ok, elapsed, n = run_one(Re, lvl)
        print(f"Re={Re} dx={lvl['dx']} (n={n} pts): {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
