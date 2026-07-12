"""
gen_ngrid3_dx001_flow_figs.py

flow_evolution.png (same style as 2-c++included/gen_flowfield_figs_v2.py)
for the two SD7003-only investigation runs into a mentor question about
weird-looking vorticity fields in 3-small_dt/:

  - ngrid=3 (multi-domain far field, dx=0.02, dt=0.005) -> SD7003/3-ngrid=3/
  - dx=0.01 (finer single grid, dt=0.005)                -> SD7003/3-dx0.01/

Both runs are Python (py/ibpm.py, native FFTW3) vs. C++ (build/ibpm), same
alpha=4.60 deg, Re=61100 flowfield case as the rest of this suite.

Usage: python3 SURF_test/gen_ngrid3_dx001_flow_figs.py
"""
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
SURF = REPO / "SURF_test"
sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

ALPHA = 4.60
DAT = SURF / "airfoils" / "LSAT-SD7003" / "sd7003.dat.txt"
VMAX = 8.0
LEVELS = np.linspace(-VMAX, VMAX, 41)


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


def load_omega(path):
    return State(filename=str(path)).omega._data[0].copy()


def draw_field(ax, field, title, X, Y, pts, alpha_deg, ok=True):
    if ok:
        ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=LEVELS, cmap="RdBu_r", extend="both")
    else:
        ax.set_facecolor("0.85")
        ax.text(0.5, 0.5, "diverged\n(NaN)", ha="center", va="center", transform=ax.transAxes,
                fontsize=9, color="0.3")
    th = -np.deg2rad(alpha_deg)
    c, sn = np.cos(th), np.sin(th)
    xc, yc = pts[:, 0] - 0.25, pts[:, 1]
    xr, yr = xc * c - yc * sn + 0.25, xc * sn + yc * c
    ax.fill(xr, yr, color="0.15", zorder=5)
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def make_flow_evolution(outdir, py_dir, cpp_dir, dt, steps, nx, ny, xoffset, yoffset,
                         subtitle, filename="flow_evolution.png"):
    outdir.mkdir(parents=True, exist_ok=True)
    length = 6.0
    dx = length / nx
    xs = xoffset + np.arange(1, nx) * dx
    ys = yoffset + np.arange(1, ny) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = load_dat_pts(DAT)

    valid_steps = []
    for s in steps:
        pyf = py_dir / f"flow{s:05d}.bin"
        cppf = cpp_dir / f"flow{s:05d}.bin"
        if pyf.exists() and cppf.exists():
            valid_steps.append(s)
    if not valid_steps:
        print(f"{outdir}: no snapshots found, skipping")
        return

    ncols = len(valid_steps)
    fig, axes = plt.subplots(2, ncols, figsize=(3.1 * ncols, 6.6))
    axes = np.atleast_2d(axes)
    for col, s in enumerate(valid_steps):
        py_field = load_omega(py_dir / f"flow{s:05d}.bin")
        cpp_field = load_omega(cpp_dir / f"flow{s:05d}.bin")
        py_ok = np.isfinite(py_field).all()
        cpp_ok = np.isfinite(cpp_field).all()
        draw_field(axes[0, col], py_field, f"py/ibpm.py, t={s * dt:g}", X, Y, pts, ALPHA, ok=py_ok)
        draw_field(axes[1, col], cpp_field, f"C++ build/ibpm, t={s * dt:g}", X, Y, pts, ALPHA, ok=cpp_ok)
        if col > 0:
            axes[0, col].set_yticklabels([]); axes[1, col].set_yticklabels([])
        axes[1, col].set_xlabel("x")
    axes[0, 0].set_ylabel("y (py/ibpm.py)")
    axes[1, 0].set_ylabel("y (C++)")
    fig.suptitle(f"SD7003: vorticity field evolution, py/ibpm.py vs. C++ build/ibpm, "
                 f"$\\alpha$={ALPHA}°, Re=61100, {subtitle}", fontsize=12)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-VMAX, VMAX), cmap="RdBu_r"),
                 ax=axes, shrink=0.7, label="vorticity $\\omega$ (clipped)")
    fig.savefig(outdir / filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir / filename} ({ncols} snapshot columns: {valid_steps})")


def main():
    # ---- ngrid=3 (dx=0.02, dt=0.005, nsteps=6000; ny=152/yoffset=-1.52 to satisfy ny%4==0) ----
    make_flow_evolution(
        outdir=SURF / "airfoils" / "LSAT-SD7003" / "3-ngrid=3",
        py_dir=SURF / "airfoils" / "LSAT-SD7003" / "_run_data_ngrid3" / "flowfield",
        cpp_dir=SURF / "airfoils" / "LSAT-SD7003" / "_run_data_ngrid3_cpp" / "flowfield",
        dt=0.005, steps=[0, 1000, 2000, 3000, 4000, 5000, 6000],
        nx=300, ny=152, xoffset=-2.0, yoffset=-1.52,
        subtitle="ngrid=3 (multi-domain), dx=0.02, dt=0.005 -- both implementations diverge to NaN "
                  "after t~20-21 (same time at dt=0.01 too -- not a CFL/dt artifact)",
    )

    # ---- dx=0.01 (dt=0.005, nsteps=6000, nx=600, ny=300) ----
    make_flow_evolution(
        outdir=SURF / "airfoils" / "LSAT-SD7003" / "3-dx0.01",
        py_dir=SURF / "airfoils" / "LSAT-SD7003" / "_run_data_dx001" / "flowfield",
        cpp_dir=SURF / "airfoils" / "LSAT-SD7003" / "_run_data_dx001_cpp" / "flowfield",
        dt=0.005, steps=[0, 1000, 2000, 3000, 4000, 5000, 6000],
        nx=600, ny=300, xoffset=-2.0, yoffset=-1.5,
        subtitle="dx=0.01, dt=0.005 (halved from dt=0.01 baseline; finest grid needs this per "
                  "SD7003/README.md's documented CFL fix) -- stable for the full t=0..30 window",
    )


if __name__ == "__main__":
    main()
