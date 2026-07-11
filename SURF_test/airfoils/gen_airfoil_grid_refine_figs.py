"""
gen_airfoil_grid_refine_figs.py

Figures for run_airfoil_grid_refine.py / run_airfoil_grid_refine_py.py:
at FIXED Re=5000, refine dx across coarse/medium/fine, comparing
py/ibpm.py vs. C++ build/ibpm at every level (2 rows x N dx columns).

Usage: python3 SURF_test/airfoils/gen_airfoil_grid_refine_figs.py <SD7003|SD8000>
Output: SURF_test/airfoils/<name>/5-grid_refine/grid_refine_comparison.png,
        SURF_test/airfoils/<name>/5-grid_refine/fidelity_summary.txt
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
    "SD7003": dict(alpha=4.60),
    "SD8000": dict(alpha=5.36),
}
LEVELS = [
    dict(tag="coarse", dx=0.04, nx=150, ny=75, final_step=3000),
    dict(tag="medium", dx=0.02, nx=300, ny=150, final_step=3000),
    dict(tag="fine", dx=0.01, nx=600, ny=300, final_step=6000),
]
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
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
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in CASES:
        print("Usage: gen_airfoil_grid_refine_figs.py <SD7003|SD8000>", file=sys.stderr)
        sys.exit(1)
    name = sys.argv[1]
    outdir = REPO / "SURF_test" / "airfoils" / name / "5-grid_refine"
    py_dir = outdir / "_run_data"
    cpp_dir = outdir / "_run_data_cpp"

    available = [lvl for lvl in LEVELS
                 if (py_dir / f"dx{lvl['dx']}" / f"run{lvl['final_step']:05d}.bin").exists()
                 and (cpp_dir / f"dx{lvl['dx']}" / f"run{lvl['final_step']:05d}.bin").exists()]
    if not available:
        print(f"{name}: no matching Python+C++ grid-refine snapshots found")
        return

    fig, axes = plt.subplots(2, len(available), figsize=(4.2 * len(available), 8.4))
    axes = np.atleast_2d(axes)
    dxs, py_rms_l, cpp_rms_l, py_max_l, cpp_max_l, rel_diff_l = [], [], [], [], [], []
    for col, lvl in enumerate(available):
        w_py = State(filename=str(py_dir / f"dx{lvl['dx']}" / f"run{lvl['final_step']:05d}.bin")).omega._data[0].copy()
        w_cpp = State(filename=str(cpp_dir / f"dx{lvl['dx']}" / f"run{lvl['final_step']:05d}.bin")).omega._data[0].copy()
        draw(axes[0, col], w_py, f"{name} dx={lvl['dx']} ({lvl['tag']})\npy/ibpm.py, t=30")
        draw(axes[1, col], w_cpp, f"{name} dx={lvl['dx']} ({lvl['tag']})\nC++ build/ibpm, t=30")
        dxs.append(lvl["dx"])
        py_rms_l.append(rms(w_py)); cpp_rms_l.append(rms(w_cpp))
        py_max_l.append(float(np.max(np.abs(w_py)))); cpp_max_l.append(float(np.max(np.abs(w_cpp))))
        rel_diff_l.append(rms(w_py - w_cpp) / max(rms(w_cpp), 1e-12))
    fig.suptitle(f"{name}: grid refinement at FIXED Re=5000 -- py/ibpm.py vs. C++ build/ibpm fidelity", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(outdir / "grid_refine_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir / 'grid_refine_comparison.png'}")

    with open(outdir / "fidelity_summary.txt", "w") as f:
        f.write(f"{name} grid refinement (Re=5000): py/ibpm.py vs. C++ build/ibpm fidelity at t=30\n\n")
        f.write(f"{'dx':>6} {'RMS_py':>10} {'RMS_cpp':>10} {'max_py':>10} {'max_cpp':>10} {'RMS(py-cpp)/RMS_cpp':>22}\n")
        for dx, rp, rc, mp, mc, rd in zip(dxs, py_rms_l, cpp_rms_l, py_max_l, cpp_max_l, rel_diff_l):
            f.write(f"{dx:6.2f} {rp:10.4f} {rc:10.4f} {mp:10.3f} {mc:10.3f} {rd:22.4%}\n")
    print(f"wrote {outdir / 'fidelity_summary.txt'}")
    for dx, rp, rc, rd in zip(dxs, py_rms_l, cpp_rms_l, rel_diff_l):
        print(f"  dx={dx}: RMS_py={rp:.4f}  RMS_cpp={rc:.4f}  rel_diff={rd:.2%}")


if __name__ == "__main__":
    main()
