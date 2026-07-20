"""
gen_faithful_figs.py

Vorticity field comparisons for the faithful single-case run (alpha=12deg,
steady, Re=1000, domain 34c x 30c, t=146) against Kurtulus (2019)'s own
Figures 2 (instantaneous) and 3 (mean), cropped in
../1-paper_based/paper_figs/.

Usage: python3 gen_faithful_figs.py
Output: figures/*.png
"""
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "4-single_case_faithful"
RUNS = HERE / "runs"
FIGS = HERE / "figures"
FIGS.mkdir(exist_ok=True)
PAPER_FIGS = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "1-paper_based" / "paper_figs"

sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py_static")
pkg.__path__ = [str(REPO / "py_static")]
sys.modules["py_static"] = pkg
from py_static.state import State  # noqa: E402

NX, NY = 1700, 1500
LENGTH, XOFFSET, YOFFSET = 34, -15, -15
DX = LENGTH / NX
xs = XOFFSET + np.arange(1, NX) * DX
ys = YOFFSET + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")
V = 8.0

MEAN_STEPS = [7300, 8760, 10220, 11680, 13140, 14600]  # t=73..146, within the paper's averaging window
FINAL_STEP = 14600


def load_omega(run_dir, step):
    return State(filename=str(run_dir / f"flow{step:05d}.bin")).omega._data[0].copy()


def load_pts():
    p = REPO / "SURF_test" / "low_re" / "NACA0012" / "1-basics" / "naca0012.dat.txt"
    return np.genfromtxt(p, skip_header=1)


def draw(ax, field, title, pts):
    ax.contourf(X, Y, np.clip(field, -V, V), levels=41, cmap="jet", extend="both")
    ax.fill(pts[:, 0], pts[:, 1], color="0.1", zorder=5)
    ax.set_xlim(-1, 4); ax.set_ylim(-1.2, 1.2); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def main():
    pts = load_pts()

    for impl in ("py", "cpp"):
        run_dir = RUNS / impl
        inst = load_omega(run_dir, FINAL_STEP)
        mean = np.mean([load_omega(run_dir, s) for s in MEAN_STEPS], axis=0)

        fig, axes = plt.subplots(2, 2, figsize=(11, 5.2))
        draw(axes[0, 0], inst, f"{impl}_static, instantaneous (t=146)", pts)
        paper_inst = PAPER_FIGS / "steady_a12.png"
        if paper_inst.exists():
            axes[0, 1].imshow(plt.imread(paper_inst)); axes[0, 1].set_title(
                "Kurtulus (2019) Fig. 2, instantaneous (t=100s)", fontsize=9)
        axes[0, 1].set_xticks([]); axes[0, 1].set_yticks([])
        for s in axes[0, 1].spines.values():
            s.set_visible(False)

        draw(axes[1, 0], mean, f"{impl}_static, mean (t=73-146)", pts)
        paper_mean = PAPER_FIGS / "steady_mean_a12.png"
        if paper_mean.exists():
            axes[1, 1].imshow(plt.imread(paper_mean)); axes[1, 1].set_title(
                "Kurtulus (2019) Fig. 3, mean (t=50-100s)", fontsize=9)
        axes[1, 1].set_xticks([]); axes[1, 1].set_yticks([])
        for s in axes[1, 1].spines.values():
            s.set_visible(False)

        fig.suptitle(f"Faithful single-case comparison: {impl}_static vs. Kurtulus (2019),\n"
                     "NACA0012, alpha=12deg, steady, Re=1000, domain 34c x 30c, t=146 "
                     "(jet: blue=-, green=0, red=+)", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        out = FIGS / f"wake_faithful_{impl}.png"
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
