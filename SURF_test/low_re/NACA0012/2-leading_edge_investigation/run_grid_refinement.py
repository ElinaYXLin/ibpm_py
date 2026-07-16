"""
run_grid_refinement.py

Test 1 of the leading-edge (LE) vorticity-speck investigation (see
../README.md "Leading-edge vorticity investigation" section).

The existing grid-convergence sweep (../run_gridconv.py) already gives
Cd(dx) at dx=0.04/0.02/0.01/0.005/0.0025, but was run with -restart 0, so
no field snapshots exist to check whether the LE vorticity "speck" visible
in the raw vorticity data (see flow_evolution.png, alpha=5, dx=0.02)
shrinks and localizes as dx drops below the NACA0012 leading-edge radius
of curvature r_LE = 1.1019*(0.12)^2 = 0.01587 =~ 0.016 (chord=1).

This script reruns Re=500, alpha=0 (matching the existing gridconv
convention) at dx=0.01 and dx=0.005 with restart snapshots enabled, at
the SAME physical times (t=0,5,...,30) for both dx so they're directly
comparable frame-by-frame, using py/ibpm.py only (fidelity to C++ is
already established at these dx in ../fidelity_summary.txt; rerunning
both implementations here would double the (already substantial: dx=0.005
took ~53 min in the original sweep) compute for no new information).

Usage: python3 SURF_test/low_re/NACA0012/leading_edge_investigation/run_grid_refinement.py [dx1 dx2 ...]
Output: SURF_test/low_re/NACA0012/leading_edge_investigation/_run_data/gridconv_dx<dx>_snap/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO / "SURF_test"))
from make_airfoil_raw import make_raw_for_dx  # noqa: E402

RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTBASE = REPO / "SURF_test" / "low_re" / "NACA0012" / "leading_edge_investigation" / "_run_data"
DAT_PATH = REPO / "SURF_test" / "low_re" / "NACA0012" / "naca0012.dat.txt"
GEOM_DIR = REPO / "SURF_test" / "geom"
DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.5)
RE = 500.0
ALPHA = 0.0
T_FINAL = 30.0
T_SNAP = 5.0  # snapshot every t=5, matching the rest of the suite's flow_evolution convention

GRID_DX = [0.01, 0.005]


def ensure_geom(dx):
    geom_path = GEOM_DIR / f"naca0012_dx{dx:.4f}.geom"
    txt_path = GEOM_DIR / f"naca0012_dx{dx:.4f}.txt"
    if not geom_path.exists():
        n, perimeter = make_raw_for_dx(str(DAT_PATH), dx, str(txt_path))
        geom_path.write_text(
            f"body NACA0012\n  raw {txt_path.relative_to(REPO)}\n  center 0.25 0.0\nend\n"
        )
        print(f"  generated {geom_path.relative_to(REPO)}: {n} points, perimeter={perimeter:.4f}")
    return geom_path


def dt_nsteps_restart_for(dx):
    dt = dx / 2  # matches ../run_gridconv.py's dt/dx=0.5 rule for dx<0.04
    nsteps = int(round(T_FINAL / dt))
    restart = int(round(T_SNAP / dt))
    return dt, nsteps, restart


def run_one(outdir, dx):
    geom_path = ensure_geom(dx)
    dt, nsteps, restart = dt_nsteps_restart_for(dx)
    nx = int(round(DOMAIN["length"] / dx))
    ny = int(round(3.0 / dx))
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-u", str(RUNNER),
           "-geom", str(geom_path), "-name", "flow", "-outdir", str(outdir),
           "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
           "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
           "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(ALPHA), "-Re", str(RE),
           "-dt", str(dt), "-nsteps", str(nsteps), "-tecplot", "0",
           "-restart", str(restart), "-force", "1"]
    log_path = outdir / "run_log.txt"
    print(f"  dx={dx}: nx={nx} ny={ny} dt={dt} nsteps={nsteps} restart={restart} -> {outdir.relative_to(REPO)}", flush=True)
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"run failed ({outdir}): see {log_path}")
    print(f"    done in {elapsed:.0f}s", flush=True)


def main():
    dxs = [float(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else GRID_DX
    for dx in dxs:
        outdir = OUTBASE / f"gridconv_dx{dx}_snap"
        if (outdir / "flow03000.bin").exists() or list(outdir.glob("flow*.bin")):
            print(f"dx={dx}: snapshots already present in {outdir.relative_to(REPO)}, skipping")
            continue
        run_one(outdir, dx)
    print("done")


if __name__ == "__main__":
    main()
