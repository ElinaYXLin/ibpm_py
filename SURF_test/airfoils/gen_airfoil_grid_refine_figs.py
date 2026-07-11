"""
gen_airfoil_grid_refine_figs.py

Figures for run_airfoil_grid_refine.py: at FIXED Re=5000 (a point in the
Re-sweep's transitional zone -- see 4-Re_sweep/domain_rms_vs_Re.png),
refine dx across coarse/medium/fine to test whether resolution alone
(not Re) controls whether the field is clean or speckled.

Usage: python3 SURF_test/airfoils/gen_airfoil_grid_refine_figs.py <SD7003|SD8000>
Output: SURF_test/airfoils/<name>/5-grid_refine/grid_refine_comparison.png
"""
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

CASES = {
    "SD7003": dict(alpha=4.60, dat=REPO / "SURF_test" / "airfoils" / "SD7003" / "sd7003.dat.txt"),
    "SD8000": dict(alpha=5.36, dat=REPO / "SURF_test" / "airfoils" / "SD8000" / "sd8000.dat.txt"),
}
LEVELS = [
    dict(tag="coarse", dx=0.04, nx=150, ny=75, final_step=3000),
    dict(tag="medium", dx=0.02, nx=300, ny=150, final_step=3000),
    dict(tag="fine", dx=0.01, nx=600, ny=300, final_step=6000),
]
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
VMAX = 8.0


def load_dat_pts(path):
    lines = pathlib.Path(path).read_text().splitlines()
    pts = []
    for l in lines[1:]:
        l = l.strip()
        if l:
            x, y = l.split()
            pts.append((float(x), float(y)))
    return np.array(pts)


def rms(f):
    return float(np.sqrt(np.nanmean(f.astype(np.float64) ** 2)))


def draw(ax, field, title, pts, alpha_deg):
    nx = field.shape[0] + 1
    dx = DOMAIN["length"] / nx
    xs = DOMAIN["xoffset"] + np.arange(1, nx) * dx
    ys = DOMAIN["yoffset"] + np.arange(1, field.shape[1] + 1) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=41, cmap="RdBu_r", extend="both")
    th = -np.deg2rad(alpha_deg)
    c, sn = np.cos(th), np.sin(th)
    xc, yc = pts[:, 0] - 0.25, pts[:, 1]
    xr, yr = xc * c - yc * sn + 0.25, xc * sn + yc * c
    ax.fill(xr, yr, color="0.15", zorder=5)
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in CASES:
        print("Usage: gen_airfoil_grid_refine_figs.py <SD7003|SD8000>", file=sys.stderr)
        sys.exit(1)
    name = sys.argv[1]
    cfg = CASES[name]
    outdir = REPO / "SURF_test" / "airfoils" / name / "5-grid_refine"
    rundir = outdir / "_run_data_cpp"
    pts = load_dat_pts(cfg["dat"])

    available = [lvl for lvl in LEVELS if (rundir / f"dx{lvl['dx']}" / f"run{lvl['final_step']:05d}.bin").exists()]
    if not available:
        print(f"{name}: no grid-refine snapshots found")
        return

    fig, axes = plt.subplots(1, len(available), figsize=(4.2 * len(available), 4.5))
    axes = np.atleast_1d(axes)
    dxs, rmses, maxes = [], [], []
    for ax, lvl in zip(axes, available):
        s = State(filename=str(rundir / f"dx{lvl['dx']}" / f"run{lvl['final_step']:05d}.bin"))
        w = s.omega._data[0].copy()
        draw(ax, w, f"{name} Re=5000, dx={lvl['dx']} ({lvl['tag']})\nt=30", pts, cfg["alpha"])
        dxs.append(lvl["dx"]); rmses.append(rms(w)); maxes.append(float(np.max(np.abs(w))))
    fig.suptitle(f"{name}: grid refinement at FIXED Re=5000 -- does finer dx alone clean up the field?", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(outdir / "grid_refine_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir / 'grid_refine_comparison.png'}")
    for dx, r, m in zip(dxs, rmses, maxes):
        print(f"  dx={dx}: RMS={r:.4f}  max={m:.3f}")


if __name__ == "__main__":
    main()
