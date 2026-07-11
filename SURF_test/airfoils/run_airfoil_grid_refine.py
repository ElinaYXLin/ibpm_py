"""
run_airfoil_grid_refine.py

Companion to vortall/run_cylinder_grid_refine.py: at a FIXED Reynolds
number (Re=5000 -- a point already run by run_airfoil_re_sweep.py at
dx=0.02, reused here rather than rerun), refines dx across the three
levels this suite already has geometry for (0.04, 0.02, 0.01) to test
whether resolution alone controls whether the field is clean or
speckled, independent of Re.

Usage: python3 SURF_test/airfoils/run_airfoil_grid_refine.py <SD7003|SD8000>
Output: SURF_test/airfoils/<name>/5-grid_refine/_run_data_cpp/dx<value>/
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
    "SD7003": dict(alpha=4.60),
    "SD8000": dict(alpha=5.36),
}

RE_FIXED = 5000
DX_LEVELS = [
    dict(tag="coarse", dx=0.04, nx=150, ny=75, dt=0.01, nsteps=3000),
    dict(tag="medium", dx=0.02, nx=300, ny=150, dt=0.01, nsteps=3000),
    dict(tag="fine", dx=0.01, nx=600, ny=300, dt=0.005, nsteps=6000),
]


def run_one(name, cfg, lvl):
    geom_path = GEOMDIR / f"{name.lower()}_dx{lvl['dx']:.4f}.geom"
    outdir = REPO / "SURF_test" / "airfoils" / name / "5-grid_refine" / "_run_data_cpp" / f"dx{lvl['dx']}"
    outdir.mkdir(parents=True, exist_ok=True)
    restart = round(lvl["nsteps"] / 6)  # 7 snapshots t=0,5,...,30
    cmd = [str(CPP_BIN),
           "-geom", str(geom_path), "-name", "run", "-outdir", str(outdir),
           "-nx", str(lvl["nx"]), "-ny", str(lvl["ny"]), "-ngrid", "1",
           "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
           "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(cfg["alpha"]), "-Re", str(RE_FIXED),
           "-dt", str(lvl["dt"]), "-nsteps", str(lvl["nsteps"]), "-tecplot", "0",
           "-restart", str(restart), "-force", "1"]
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    return proc.returncode == 0, elapsed


def main():
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr)
        sys.exit(1)
    if len(sys.argv) < 2 or sys.argv[1] not in CASES:
        print("Usage: run_airfoil_grid_refine.py <SD7003|SD8000> [tag]", file=sys.stderr)
        sys.exit(1)
    name = sys.argv[1]
    cfg = CASES[name]
    tag = sys.argv[2] if len(sys.argv) > 2 else None
    levels = [l for l in DX_LEVELS if tag is None or l["tag"] == tag]
    for lvl in levels:
        ok, elapsed = run_one(name, cfg, lvl)
        print(f"{name} dx={lvl['dx']}: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
