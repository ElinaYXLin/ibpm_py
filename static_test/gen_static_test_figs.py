"""
gen_static_test_figs.py

Figures for the static-algorithm reproducibility check (static_test/README.md):
py_static/ibpm.py and build_static/ibpm (cpp_static), run 5x each on the
standard NACA0012 flow-evolution case (Re=500, alpha=5deg, dx=0.02, dt=0.01,
nsteps=3000, restart=500 -> 7 snapshots t=0,5,...,30 -- same convention as
SURF_test/low_re/run_naca0012.py), after switching the DST planner flag from
FFTW_EXHAUSTIVE to FFTW_ESTIMATE|FFTW_UNALIGNED so the same sine-transform
algorithm is used every run, in both languages.

Three figures:
- flow_evolution_py_vs_cpp.png: the standard 2-row vorticity field comparison
  (py_run1 top, cpp_run1 bottom), same style as the rest of this suite.
- reproducibility_diff.png: (run5 - run1) vorticity difference, one row per
  implementation. Every run was already confirmed bit-identical to run1
  (max|diff| == 0.0 exactly for all 5 runs, both implementations, all 7
  snapshots) -- these panels are expected to render as flat/blank fields at
  fixed color scale, which is itself the result being shown.
- py_vs_cpp_diff.png: (py_run1 - cpp_run1) vorticity difference at each
  snapshot, on a roundoff-scale color axis (~1e-12), showing the tiny
  floating-point-level disagreement between the two independent
  implementations running the same fixed algorithm.

Usage: python3 static_test/gen_static_test_figs.py
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
pkg = types.ModuleType("py_static")
pkg.__path__ = [str(REPO / "py_static")]
sys.modules["py_static"] = pkg
from py_static.state import State  # noqa: E402

STATIC_TEST = REPO / "static_test"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
NX, NY = 300, 150
DX = DOMAIN["length"] / NX
xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")
STEPS = [0, 500, 1000, 1500, 2000, 2500, 3000]
VMAX = 8.0
DAT_PATH = REPO / "SURF_test" / "low_re" / "NACA0012" / "naca0012.dat.txt"


def load_dat_pts(path):
    lines = pathlib.Path(path).read_text().splitlines()
    pts = []
    for l in lines[1:]:
        l = l.strip()
        if l:
            x, y = l.split()
            pts.append((float(x), float(y)))
    return np.array(pts)


def load_omega(run_dir, step):
    f = run_dir / f"flow{step:05d}.bin"
    return State(filename=str(f)).omega._data[0].copy()


def draw_field(ax, field, title, pts, vmax, cmap="RdBu_r"):
    ax.contourf(X, Y, np.clip(field, -vmax, vmax), levels=41, cmap=cmap, extend="both")
    ax.fill(pts[:, 0], pts[:, 1], color="0.15", zorder=5)
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def main():
    pts = load_dat_pts(DAT_PATH)
    py1 = {s: load_omega(STATIC_TEST / "py_run1" / "flowfield", s) for s in STEPS}
    cpp1 = {s: load_omega(STATIC_TEST / "cpp_run1" / "flowfield", s) for s in STEPS}
    py5 = {s: load_omega(STATIC_TEST / "py_run5" / "flowfield", s) for s in STEPS}
    cpp5 = {s: load_omega(STATIC_TEST / "cpp_run5" / "flowfield", s) for s in STEPS}

    # ---- Figure A: standard py vs cpp flow evolution ----
    ncols = len(STEPS)
    fig, axes = plt.subplots(2, ncols, figsize=(3.1 * ncols, 6.6))
    for col, s in enumerate(STEPS):
        draw_field(axes[0, col], py1[s], f"py_static run1, t={s*0.01:g}", pts, VMAX)
        draw_field(axes[1, col], cpp1[s], f"cpp_static run1, t={s*0.01:g}", pts, VMAX)
        if col > 0:
            axes[0, col].set_yticklabels([]); axes[1, col].set_yticklabels([])
        axes[1, col].set_xlabel("x")
    axes[0, 0].set_ylabel("y (py_static)"); axes[1, 0].set_ylabel("y (cpp_static)")
    fig.suptitle("NACA0012, Re=500, $\\alpha$=5° (dx=0.02): py_static vs cpp_static "
                 "vorticity evolution, fixed-algorithm DST (FFTW_ESTIMATE|FFTW_UNALIGNED)", fontsize=11)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-VMAX, VMAX), cmap="RdBu_r"),
                 ax=axes, shrink=0.7, label="vorticity $\\omega$ (clipped)")
    out = STATIC_TEST / "flow_evolution_py_vs_cpp.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # ---- Figure B: reproducibility diff (run5 - run1), one row per implementation ----
    fig, axes = plt.subplots(2, ncols, figsize=(3.1 * ncols, 6.6))
    rvmax = 1e-10  # fixed tiny scale -- these diffs are exactly 0.0, so panels render flat
    for col, s in enumerate(STEPS):
        dpy = py5[s] - py1[s]
        dcpp = cpp5[s] - cpp1[s]
        draw_field(axes[0, col], dpy, f"py: run5-run1, t={s*0.01:g}\nmax|d|={np.max(np.abs(dpy)):.1e}",
                   pts, rvmax, cmap="PuOr")
        draw_field(axes[1, col], dcpp, f"cpp: run5-run1, t={s*0.01:g}\nmax|d|={np.max(np.abs(dcpp)):.1e}",
                   pts, rvmax, cmap="PuOr")
        if col > 0:
            axes[0, col].set_yticklabels([]); axes[1, col].set_yticklabels([])
        axes[1, col].set_xlabel("x")
    axes[0, 0].set_ylabel("y (py_static)"); axes[1, 0].set_ylabel("y (cpp_static)")
    fig.suptitle("Run-to-run reproducibility: vorticity difference, run5 minus run1\n"
                 "(all 5 runs are bit-identical in both implementations -- panels are exactly zero, "
                 f"color axis fixed at ±{rvmax:g} for scale)", fontsize=11)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-rvmax, rvmax), cmap="PuOr"),
                 ax=axes, shrink=0.7, label="$\\Delta\\omega$ (run5 - run1)")
    out = STATIC_TEST / "reproducibility_diff.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # ---- Figure C: py vs cpp diff at roundoff scale ----
    fig, axes = plt.subplots(1, ncols, figsize=(3.1 * ncols, 3.6))
    dvmax = 2e-12
    for col, s in enumerate(STEPS):
        d = py1[s] - cpp1[s]
        draw_field(axes[col], d, f"t={s*0.01:g}\nmax|d|={np.max(np.abs(d)):.2e}", pts, dvmax, cmap="PRGn")
        if col > 0:
            axes[col].set_yticklabels([])
        axes[col].set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.suptitle("py_static run1 minus cpp_static run1: vorticity difference "
                 "(floating-point-roundoff scale)", fontsize=11)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-dvmax, dvmax), cmap="PRGn"),
                 ax=axes, shrink=0.7, label="$\\Delta\\omega$ (py - cpp)")
    out = STATIC_TEST / "py_vs_cpp_diff.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
