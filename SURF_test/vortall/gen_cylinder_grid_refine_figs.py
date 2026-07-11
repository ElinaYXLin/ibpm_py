"""
gen_cylinder_grid_refine_figs.py

Figures for run_cylinder_grid_refine.py: at FIXED Re=5000, refine dx
across coarse/medium/fine. Companion to
airfoils/gen_airfoil_grid_refine_figs.py, same style/question.

Usage: python3 SURF_test/vortall/gen_cylinder_grid_refine_figs.py
Output: SURF_test/vortall/3-grid_refine/grid_refine_comparison.png
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

OUTDIR = REPO / "SURF_test" / "vortall" / "3-grid_refine"
RUNDIR = OUTDIR / "_run_data_cpp"
DOMAIN = dict(length=9, xoffset=-1, yoffset=-2)
LEVELS = [
    dict(tag="coarse", dx=0.04, final_step=1500),
    dict(tag="medium", dx=0.02, final_step=6000),
    dict(tag="fine", dx=0.01, final_step=12000),
]
VMAX = 8.0


def rms(f):
    return float(np.sqrt(np.nanmean(f.astype(np.float64) ** 2)))


def draw(ax, field, title):
    nx = field.shape[0] + 1
    dx = DOMAIN["length"] / nx
    xs = DOMAIN["xoffset"] + np.arange(1, nx) * dx
    ys = DOMAIN["yoffset"] + np.arange(1, field.shape[1] + 1) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=41, cmap="RdBu_r", extend="both")
    circ = plt.Circle((0, 0), 0.5, color="0.15", zorder=5)
    ax.add_patch(circ)
    ax.set_xlim(-1, 8); ax.set_ylim(-2, 2); ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)


def main():
    available = [lvl for lvl in LEVELS
                 if (RUNDIR / f"dx{lvl['dx']}" / f"cyl{lvl['final_step']:05d}.bin").exists()]
    if not available:
        print("no cylinder grid-refine snapshots found")
        return

    fig, axes = plt.subplots(1, len(available), figsize=(4.8 * len(available), 4.5))
    axes = np.atleast_1d(axes)
    for ax, lvl in zip(axes, available):
        s = State(filename=str(RUNDIR / f"dx{lvl['dx']}" / f"cyl{lvl['final_step']:05d}.bin"))
        w = s.omega._data[0].copy()
        draw(ax, w, f"cylinder Re=5000, dx={lvl['dx']} ({lvl['tag']})\nt=30")
        print(f"  dx={lvl['dx']}: RMS={rms(w):.4f}  max={float(np.max(np.abs(w))):.3f}")
    fig.suptitle("Cylinder: grid refinement at FIXED Re=5000 -- does finer dx alone clean up the field?", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(OUTDIR / "grid_refine_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 'grid_refine_comparison.png'}")


if __name__ == "__main__":
    main()
