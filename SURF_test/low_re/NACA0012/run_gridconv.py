"""
run_gridconv.py

Drives the NACA0012 Re=500, alpha=0 grid-convergence sweep (both
py/ibpm.py and C++ build/ibpm) that backs grid_convergence.png /
fidelity_summary.txt. The original dx=0.04/0.02/0.01 points were produced
by a one-off script that wasn't preserved in this repo; this rewrite
reproduces their exact settings (verified: rerunning dx=0.01 here
reproduces its committed run.force to the mean Cd/Cl already on record)
and generalizes them so further grid-refinement levels can be added by
just extending GRID_DX below.

Convention, matching the existing dx=0.04/0.02/0.01 points:
  - domain: length=6, xoffset=-2, yoffset=-1.5 (same as the polar runs)
  - nx = length/dx, ny = (yoffset range = 3)/dx
  - dt = dx/2 (matches the dt/dx=0.5 ratio already used at dx=0.02 and
    dx=0.01; dx=0.04 used a fixed dt=0.01 since it didn't need anything
    finer, but dt=dx/2=0.02 would also be CFL-safe there, so this rule
    is only applied going forward for dx < 0.04)
  - nsteps = 30/dt, so every dx reaches the same t=30
  - boundary points: raw NACA0012 point file resampled to ds=dx via
    make_airfoil_raw.py (SURF_test/geom/naca0012_dx<dx>.geom/.txt),
    same convention as the existing geometry files

Usage: python3 SURF_test/low_re/NACA0012/run_gridconv.py [dx1 dx2 ...]
       (defaults to running any of GRID_DX not yet present)
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO / "SURF_test"))
from make_airfoil_raw import make_raw_for_dx  # noqa: E402

CPP_BIN = REPO / "build" / "ibpm"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTBASE = REPO / "SURF_test" / "low_re" / "NACA0012"
DAT_PATH = OUTBASE / "naca0012.dat.txt"
GEOM_DIR = REPO / "SURF_test" / "geom"
DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.5)
RE = 500.0
ALPHA = 0.0
T_FINAL = 30.0

# Grid spacings to have data for, each half the previous (existing:
# 0.04, 0.02, 0.01; extend this list to add finer levels).
GRID_DX = [0.04, 0.02, 0.01, 0.005, 0.0025]


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


def dt_nsteps_for(dx):
    # dt/dx = 0.5 for dx < 0.04 (matches the existing dx=0.02, dx=0.01
    # points); dx=0.04 keeps its existing dt=0.01 (see module docstring).
    dt = 0.01 if dx >= 0.04 else dx / 2
    nsteps = int(round(T_FINAL / dt))
    return dt, nsteps


def run_one(cmd_prefix, outdir, dx):
    geom_path = ensure_geom(dx)
    dt, nsteps = dt_nsteps_for(dx)
    nx = int(round(DOMAIN["length"] / dx))
    ny = int(round(3.0 / dx))  # yoffset=-1.5 -> y in [-1.5, 1.5], height 3
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = cmd_prefix + [
        "-geom", str(geom_path), "-name", "run", "-outdir", str(outdir),
        "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(ALPHA), "-Re", str(RE),
        "-dt", str(dt), "-nsteps", str(nsteps), "-tecplot", "0", "-restart", "0", "-force", "1",
    ]
    log_path = outdir / "run_log.txt"
    print(f"  dx={dx}: nx={nx} ny={ny} dt={dt} nsteps={nsteps} -> {outdir.relative_to(REPO)}")
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"run failed ({outdir}): see {log_path}")
    print(f"    done in {elapsed:.0f}s")
    return outdir / "run.force"


def mean_force(force_path, frac=0.6):
    import numpy as np
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    seg = d[int(len(d) * (1 - frac)):]
    return float(seg[:, 3].mean()), float(seg[:, 2].mean())  # cl, cd


def main():
    dxs = [float(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else GRID_DX
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr)
        sys.exit(1)

    for dx in dxs:
        for impl, cmd_prefix, subdir in (
            ("cpp", [str(CPP_BIN)], "_run_data_cpp"),
            ("py", [sys.executable, "-u", str(RUNNER)], "_run_data"),
        ):
            outdir = OUTBASE / subdir / f"gridconv_dx{dx}"
            fpath = outdir / "run.force"
            if fpath.exists():
                cl, cd = mean_force(fpath)
                print(f"[{impl}] dx={dx}: already present, Cl={cl:+.4f} Cd={cd:+.4f} (skipping)")
                continue
            print(f"[{impl}] dx={dx}: running...")
            fpath = run_one(cmd_prefix, outdir, dx)
            cl, cd = mean_force(fpath)
            print(f"[{impl}] dx={dx}: Cl={cl:+.4f}  Cd={cd:+.4f}")

    print("done")


if __name__ == "__main__":
    main()
