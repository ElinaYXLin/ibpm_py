"""
gen_static_test_figs_sd8000.py

Figures for the SD8000 coarse-grid (dx=0.04, Re=60800, "hardest" -- see
static_test/SD8000/README.md) reproducibility check, mirroring
static_test/gen_static_test_figs.py's three-figure structure but adapted
for this case's numerical instability: py_static blows up to NaN at step
2827 (t=28.27) in every one of its 5 runs (bit-identically -- see
README.md), while cpp_static stays bounded through t=30 in every run.

- flow_evolution_py_vs_cpp.png: 2-row vorticity field evolution, py_run1
  (top) vs. cpp_run1 (bottom); py's t=30 panel is replaced with a "blown
  up (NaN)" placeholder.
- reproducibility_diff.png: (run5 - run1) vorticity difference, one row
  per implementation. All panels are exactly zero (confirmed by byte
  comparison), including the NaN panel (NaN in run5 == NaN in run1,
  bit-identically, at every poisoned array element).
- py_vs_cpp_diff.png: (py_run1 - cpp_run1) vorticity difference at each
  pre-blow-up snapshot, log-magnitude annotated per panel -- shows the
  difference growing from ~1e-12 (t=2.5, roundoff floor) to ~350 (t=25,
  order of the vorticity field itself) as the two implementations'
  trajectories chaotically diverge, consistent with
  SURF_test/gen_port_fidelity_diagnostic.py's Panel C finding.

Usage: python3 static_test/SD8000/gen_static_test_figs_sd8000.py
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

BASE = REPO / "static_test" / "SD8000"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
NX, NY = 150, 75
DX = DOMAIN["length"] / NX
xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")
STEPS = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]
VMAX = 8.0
DAT_PATH = REPO / "SURF_test" / "airfoils" / "LSAT-SD8000" / "sd8000.dat.txt"


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
    return State(filename=str(run_dir / f"flow{step:05d}.bin")).omega._data[0].copy()


def draw_field(ax, field, title, pts, vmax, cmap="RdBu_r"):
    ax.contourf(X, Y, np.clip(field, -vmax, vmax), levels=41, cmap=cmap, extend="both")
    ax.fill(pts[:, 0], pts[:, 1], color="0.15", zorder=5)
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def draw_nan_panel(ax, title, pts):
    ax.set_facecolor("0.85")
    for y0 in np.linspace(-1.5, 1.5, 12):
        ax.plot([-2, 4], [y0, y0 + 3], color="0.6", lw=0.5, zorder=1)
    ax.fill(pts[:, 0], pts[:, 1], color="0.15", zorder=5)
    ax.text(1, 0, "blown up\n(NaN)", ha="center", va="center", fontsize=10,
            color="firebrick", zorder=6, weight="bold")
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def main():
    pts = load_dat_pts(DAT_PATH)
    py1dir = BASE / "py_run1" / "flowfield"
    cpp1dir = BASE / "cpp_run1" / "flowfield"
    py5dir = BASE / "py_run5" / "flowfield"
    cpp5dir = BASE / "cpp_run5" / "flowfield"

    py1 = {s: load_omega(py1dir, s) for s in STEPS}
    cpp1 = {s: load_omega(cpp1dir, s) for s in STEPS}
    py5 = {s: load_omega(py5dir, s) for s in STEPS}
    cpp5 = {s: load_omega(cpp5dir, s) for s in STEPS}

    # ---- Figure A: standard py vs cpp flow evolution ----
    ncols = len(STEPS)
    fig, axes = plt.subplots(2, ncols, figsize=(2.6 * ncols, 5.8))
    for col, s in enumerate(STEPS):
        if np.isnan(py1[s]).any():
            draw_nan_panel(axes[0, col], f"py_static run1, t={s*0.01:g}", pts)
        else:
            draw_field(axes[0, col], py1[s], f"py_static run1, t={s*0.01:g}", pts, VMAX)
        draw_field(axes[1, col], cpp1[s], f"cpp_static run1, t={s*0.01:g}", pts, VMAX)
        if col > 0:
            axes[0, col].set_yticklabels([]); axes[1, col].set_yticklabels([])
        axes[1, col].set_xlabel("x")
    axes[0, 0].set_ylabel("y (py_static)"); axes[1, 0].set_ylabel("y (cpp_static)")
    fig.suptitle("SD8000 coarse grid (dx=0.04), Re=60800, $\\alpha$=-0.81°: py_static vs cpp_static\n"
                 "vorticity evolution -- py_static blows up to NaN at step 2827 in every run "
                 "(fixed-algorithm DST, FFTW_ESTIMATE|FFTW_UNALIGNED)", fontsize=10)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-VMAX, VMAX), cmap="RdBu_r"),
                 ax=axes, shrink=0.7, label="vorticity $\\omega$ (clipped)")
    out = BASE / "flow_evolution_py_vs_cpp.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # ---- Figure B: reproducibility diff (run5 - run1) ----
    fig, axes = plt.subplots(2, ncols, figsize=(2.6 * ncols, 5.8))
    rvmax = 1e-10
    for col, s in enumerate(STEPS):
        dpy = py5[s] - py1[s]
        dcpp = cpp5[s] - cpp1[s]
        if np.isnan(py1[s]).any():
            same = np.array_equal(np.isnan(py5[s]), np.isnan(py1[s])) and \
                   np.array_equal(py5[s][~np.isnan(py5[s])], py1[s][~np.isnan(py1[s])])
            draw_nan_panel(axes[0, col], f"py: run5-run1, t={s*0.01:g}\n"
                           f"{'both NaN, bit-identical' if same else 'MISMATCH'}", pts)
        else:
            draw_field(axes[0, col], dpy, f"py: run5-run1, t={s*0.01:g}\nmax|d|={np.max(np.abs(dpy)):.1e}",
                       pts, rvmax, cmap="PuOr")
        draw_field(axes[1, col], dcpp, f"cpp: run5-run1, t={s*0.01:g}\nmax|d|={np.max(np.abs(dcpp)):.1e}",
                   pts, rvmax, cmap="PuOr")
        if col > 0:
            axes[0, col].set_yticklabels([]); axes[1, col].set_yticklabels([])
        axes[1, col].set_xlabel("x")
    axes[0, 0].set_ylabel("y (py_static)"); axes[1, 0].set_ylabel("y (cpp_static)")
    fig.suptitle("Run-to-run reproducibility: vorticity difference, run5 minus run1\n"
                 "(all 5 runs are bit-identical in both implementations, even through py_static's "
                 f"NaN blow-up -- color axis fixed at ±{rvmax:g} for scale)", fontsize=10)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-rvmax, rvmax), cmap="PuOr"),
                 ax=axes, shrink=0.7, label="$\\Delta\\omega$ (run5 - run1)")
    out = BASE / "reproducibility_diff.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # ---- Figure C: py vs cpp diff, adaptive scale per panel, pre-blow-up steps only ----
    pre_steps = [s for s in STEPS if not np.isnan(py1[s]).any()]
    fig, axes = plt.subplots(1, len(pre_steps), figsize=(2.6 * len(pre_steps), 3.2))
    axes = np.atleast_1d(axes)
    for col, s in enumerate(pre_steps):
        d = py1[s] - cpp1[s]
        dmax = max(np.max(np.abs(d)), 1e-15)
        draw_field(axes[col], d, f"t={s*0.01:g}\nmax|d|={dmax:.2e}", pts, dmax, cmap="PRGn")
        if col > 0:
            axes[col].set_yticklabels([])
        axes[col].set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.suptitle("py_static run1 minus cpp_static run1: vorticity difference\n"
                 "(chaotic amplification of roundoff -- color axis rescaled per panel, "
                 "see title for actual magnitude; py_static blows up at t=28.27)", fontsize=10)
    out = BASE / "py_vs_cpp_diff.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
