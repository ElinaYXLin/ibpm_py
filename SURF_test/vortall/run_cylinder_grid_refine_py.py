"""
run_cylinder_grid_refine_py.py

Python counterpart of run_cylinder_grid_refine.py -- same Re, same dx
levels/dt/nsteps, reusing the exact same geometry files already written
to 3-grid_refine/geom/ by the C++ pass.

Usage: python3 SURF_test/vortall/run_cylinder_grid_refine_py.py <Re> [tag]
Output: SURF_test/vortall/3-grid_refine/_run_data/dx<value>/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTBASE = REPO / "SURF_test" / "vortall" / "3-grid_refine"
GEOMDIR = OUTBASE / "geom"

DOMAIN = dict(length=9, xoffset=-1, yoffset=-2)
DX_LEVELS = [
    dict(tag="coarse", dx=0.04, nx=225, ny=100, dt=0.02, nsteps=1500),
    dict(tag="medium", dx=0.02, nx=450, ny=200, dt=0.005, nsteps=6000),
    dict(tag="fine", dx=0.01, nx=900, ny=400, dt=0.0025, nsteps=12000),
]
RESTART_FRAC = 1 / 6


def run_one(Re, lvl):
    geom_path = GEOMDIR / f"cylinder_dx{lvl['dx']:.4f}.geom"
    if not geom_path.exists():
        raise FileNotFoundError(f"{geom_path} missing -- run run_cylinder_grid_refine.py (C++) first")
    outdir = OUTBASE / "_run_data" / f"dx{lvl['dx']}"
    outdir.mkdir(parents=True, exist_ok=True)
    restart = round(lvl["nsteps"] * RESTART_FRAC)
    cmd = [sys.executable, "-u", str(RUNNER),
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
    return proc.returncode == 0, elapsed


def main():
    Re = float(sys.argv[1]) if len(sys.argv) > 1 else 3000
    tag = sys.argv[2] if len(sys.argv) > 2 else None
    levels = [l for l in DX_LEVELS if tag is None or l["tag"] == tag]
    for lvl in levels:
        ok, elapsed = run_one(Re, lvl)
        print(f"Re={Re} dx={lvl['dx']}: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
