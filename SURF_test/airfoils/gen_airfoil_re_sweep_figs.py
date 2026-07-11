"""
gen_airfoil_re_sweep_figs.py

Figures for run_airfoil_re_sweep.py's Re sweep (E1+E2 combined): grid of
(Re rows) x (vorticity | velocity magnitude columns) at t=30, plus a
quantitative domain-RMS/max vorticity vs Re companion plot. Companion to
vortall/gen_cylinder_re_sweep_figs.py, same style, opposite Re direction.

Usage: python3 SURF_test/airfoils/gen_airfoil_re_sweep_figs.py <SD7003|SD8000>
Output: SURF_test/airfoils/<name>/4-Re_sweep/re_sweep_comparison.png,
        SURF_test/airfoils/<name>/4-Re_sweep/domain_rms_vs_Re.png
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
from py.vector_operations import FluxToXVelocity, FluxToYVelocity  # noqa: E402
from py.scalar import Scalar  # noqa: E402

CASES = {
    "SD7003": dict(alpha=4.60, dat=REPO / "SURF_test" / "airfoils" / "SD7003" / "sd7003.dat.txt"),
    "SD8000": dict(alpha=5.36, dat=REPO / "SURF_test" / "airfoils" / "SD8000" / "sd8000.dat.txt"),
}
RE_VALUES = [200, 500, 1000, 5000, 10000, 20000, 40000]
FINAL_STEP = 3000  # t=30
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
NX, NY = 300, 150
DX = DOMAIN["length"] / NX
xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")


def load_dat_pts(path):
    lines = pathlib.Path(path).read_text().splitlines()
    pts = []
    for l in lines[1:]:
        l = l.strip()
        if not l:
            continue
        x, y = l.split()
        pts.append((float(x), float(y)))
    return np.array(pts)


def vorticity(s):
    return s.omega._data[0].copy()


def velocity_mag(s):
    u = Scalar(s.q.getGrid())
    v = Scalar(s.q.getGrid())
    FluxToXVelocity(s.q, u)
    FluxToYVelocity(s.q, v)
    return np.sqrt(u._data[0] ** 2 + v._data[0] ** 2)


def rms(f):
    return float(np.sqrt(np.nanmean(f.astype(np.float64) ** 2)))


def draw(ax, field, title, pts, alpha_deg, vmax, cmap):
    ax.contourf(X, Y, np.clip(field, -vmax if cmap == "RdBu_r" else 0, vmax),
                levels=41, cmap=cmap, extend="both")
    th = -np.deg2rad(alpha_deg)
    c, sn = np.cos(th), np.sin(th)
    xc, yc = pts[:, 0] - 0.25, pts[:, 1]
    xr, yr = xc * c - yc * sn + 0.25, xc * sn + yc * c
    ax.fill(xr, yr, color="0.15", zorder=5)
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in CASES:
        print("Usage: gen_airfoil_re_sweep_figs.py <SD7003|SD8000>", file=sys.stderr)
        sys.exit(1)
    name = sys.argv[1]
    cfg = CASES[name]
    outdir = REPO / "SURF_test" / "airfoils" / name / "4-Re_sweep"
    rundir = outdir / "_run_data_cpp"
    pts = load_dat_pts(cfg["dat"])

    available = [Re for Re in RE_VALUES if (rundir / f"Re{Re}" / f"run{FINAL_STEP:05d}.bin").exists()]
    if not available:
        print(f"{name}: no completed Re-sweep snapshots found")
        return

    fig, axes = plt.subplots(len(available), 2, figsize=(9, 2.6 * len(available)))
    axes = np.atleast_2d(axes)
    re_list, rms_list, max_list = [], [], []
    for row, Re in enumerate(available):
        s = State(filename=str(rundir / f"Re{Re}" / f"run{FINAL_STEP:05d}.bin"))
        w = vorticity(s)
        vmag = velocity_mag(s)
        draw(axes[row, 0], w, f"Re={Re}: vorticity, t=30", pts, cfg["alpha"], vmax=8.0, cmap="RdBu_r")
        draw(axes[row, 1], vmag, f"Re={Re}: |velocity|, t=30", pts, cfg["alpha"], vmax=2.0, cmap="viridis")
        re_list.append(Re)
        rms_list.append(rms(w))
        max_list.append(float(np.max(np.abs(w))))
    fig.suptitle(f"{name} Re sweep: does lowering Re from {name}'s usual ~40-61k\n"
                 f"toward the cylinder's clean Re=100 baseline clean up the wake?", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outdir / "re_sweep_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir / 're_sweep_comparison.png'} ({len(available)} Re values)")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(re_list, rms_list, "o-", color="C0")
    ax[0].set_xscale("log"); ax[0].set_xlabel("Re"); ax[0].set_ylabel(r"domain-RMS $|\omega|$, t=30")
    ax[0].set_title(f"{name}: domain-RMS vorticity vs. Re"); ax[0].grid(alpha=0.3)
    ax[1].plot(re_list, max_list, "s-", color="C3")
    ax[1].set_xscale("log"); ax[1].set_xlabel("Re"); ax[1].set_ylabel(r"max $|\omega|$, t=30")
    ax[1].set_title(f"{name}: peak vorticity vs. Re"); ax[1].grid(alpha=0.3)
    fig.suptitle(f"{name} Re sweep: quantitative trend of broadband noise vs. Re", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(outdir / "domain_rms_vs_Re.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir / 'domain_rms_vs_Re.png'}")
    print(f"\n{name} Re, RMS|omega|, max|omega|:")
    for Re, r, m in zip(re_list, rms_list, max_list):
        print(f"  Re={Re}: RMS={r:.4f}  max={m:.4f}")


if __name__ == "__main__":
    main()
