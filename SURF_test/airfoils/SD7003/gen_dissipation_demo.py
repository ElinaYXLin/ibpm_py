"""
gen_dissipation_demo.py (E5: "explicit dissipation" question)

NOT a solver modification. py/ibpm.py and build/ibpm are DNS-style
solvers with no subgrid/turbulence model by design (faithful ports of the
unmodified cwrowley/ibpm method) -- actually adding a dissipation term
would mean forking the timestepper (touching the RK3 substep / elliptic
solve chain), which is a much larger, riskier change than this repo's
other diagnostic experiments and would no longer be "the same solver"
whose fidelity this whole test suite otherwise validates.

Instead, this is a cheap, honest proxy: take the ALREADY-COMPUTED,
speckled Re=40000 vorticity snapshot (4-Re_sweep/_run_data_cpp/Re40000)
and apply a small Gaussian spatial filter post-hoc, purely as a
visualization -- "if this field had subgrid dissipation damping
grid-scale content, roughly this is what it would look like." This does
NOT re-run the physics and proves nothing about solver correctness; it
only illustrates that the speckle is high-spatial-frequency content that
a low-pass filter removes, consistent with (not new evidence beyond) the
delta/dx resolution story the Re-sweep already established.

Usage: python3 SURF_test/airfoils/SD7003/gen_dissipation_demo.py
Output: SURF_test/airfoils/SD7003/6-explicit_dissipation/filter_demo.png
"""
import pathlib
import sys
import types

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

SRC = REPO / "SURF_test" / "airfoils" / "SD7003" / "4-Re_sweep" / "_run_data_cpp" / "Re40000" / "run03000.bin"
OUTDIR = REPO / "SURF_test" / "airfoils" / "SD7003" / "6-explicit_dissipation"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
NX, NY = 300, 150
DX = DOMAIN["length"] / NX
xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")
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


def draw(ax, field, title, pts, alpha_deg):
    ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=41, cmap="RdBu_r", extend="both")
    th = -np.deg2rad(alpha_deg)
    c, sn = np.cos(th), np.sin(th)
    xc, yc = pts[:, 0] - 0.25, pts[:, 1]
    xr, yr = xc * c - yc * sn + 0.25, xc * sn + yc * c
    ax.fill(xr, yr, color="0.15", zorder=5)
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pts = load_dat_pts(REPO / "SURF_test" / "airfoils" / "SD7003" / "sd7003.dat.txt")
    w = State(filename=str(SRC)).omega._data[0].copy()

    sigmas = [0, 1, 2, 4]  # grid cells
    fig, axes = plt.subplots(1, len(sigmas), figsize=(4 * len(sigmas), 4.2))
    for ax, sigma in zip(axes, sigmas):
        f = w if sigma == 0 else gaussian_filter(w, sigma=sigma)
        label = "raw (Re=40000, t=30)" if sigma == 0 else f"Gaussian-filtered, sigma={sigma} cells"
        draw(ax, f, label, pts, 4.60)
    fig.suptitle("SD7003 Re=40000: post-hoc spatial filtering of the already-computed field\n"
                 "(NOT a solver rerun -- illustrates that the speckle is high-spatial-frequency "
                 "content a low-pass filter removes)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(OUTDIR / "filter_demo.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 'filter_demo.png'}")

    for sigma in sigmas:
        f = w if sigma == 0 else gaussian_filter(w, sigma=sigma)
        print(f"  sigma={sigma}: RMS={np.sqrt(np.mean(f**2)):.4f}  max={np.max(np.abs(f)):.3f}")


if __name__ == "__main__":
    main()
