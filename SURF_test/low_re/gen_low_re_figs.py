"""
gen_low_re_figs.py

Flow-evolution figures for the two genuinely-low-Re (hundreds) airfoils:
- NACA0012 (Re=500, alpha=5 deg) -- "easy": thin, symmetric, no camber,
  minimal laminar-separation-bubble risk. Freshly run for this folder
  (SURF_test/low_re/run_naca0012.py).
- SD7003 (Re=500, alpha=4.60 deg) -- "hard": the same cambered,
  laminar-separation-bubble-prone airfoil that originally prompted the
  mentor's question, now shown at genuinely low Re instead of its usual
  ~61,100. Reuses SD7003's own Re=500 run from
  ../airfoils/LSAT-SD7003/4-Re_sweep/ rather than duplicating it.

Both: py/ibpm.py vs. C++ build/ibpm, same flow_evolution.png convention
(7 snapshots t=0,5,...,30) as the rest of this suite.

Usage: python3 SURF_test/low_re/gen_low_re_figs.py
Output: SURF_test/low_re/NACA0012/flow_evolution.png,
        SURF_test/low_re/SD7003/flow_evolution.png
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

DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
NX, NY = 300, 150
DX = DOMAIN["length"] / NX
xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")
VMAX = 8.0
STEPS = [0, 500, 1000, 1500, 2000, 2500, 3000]


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
    # NOTE(fix): py/ibpm.py's `-alpha` only tilts the free-stream (BaseFlow);
    # it is never applied to the geometry (confirmed: `geom.load(...)` in
    # ibpm.py never sees alpha). The solved vorticity field `field` is
    # therefore in the frame where the body sits at its raw, UNROTATED
    # orientation. This function used to rotate only the drawn outline by
    # -alpha_deg without rotating `field` to match -- inconsistent, and
    # equivalent to plotting the body about one grid cell away from where
    # it truly sits (at this dx, alpha) in the field being drawn. Plotting
    # the raw points directly (no rotation) matches what was actually
    # solved. `alpha_deg` is kept as a parameter for call-site compatibility
    # but is no longer used for a rotation.
    ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=41, cmap="RdBu_r", extend="both")
    ax.fill(pts[:, 0], pts[:, 1], color="0.15", zorder=5)
    ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def make_flow_evolution(label, alpha, dat_path, py_dir, cpp_dir, filename_stem, out_path, re_val):
    pts = load_dat_pts(dat_path)
    steps = [s for s in STEPS if (py_dir / f"{filename_stem}{s:05d}.bin").exists()
             and (cpp_dir / f"{filename_stem}{s:05d}.bin").exists()]
    if not steps:
        print(f"{label}: no matching snapshots found, skipping")
        return
    ncols = len(steps)
    fig, axes = plt.subplots(2, ncols, figsize=(3.1 * ncols, 6.6))
    axes = np.atleast_2d(axes)
    for col, s in enumerate(steps):
        w_py = State(filename=str(py_dir / f"{filename_stem}{s:05d}.bin")).omega._data[0].copy()
        w_cpp = State(filename=str(cpp_dir / f"{filename_stem}{s:05d}.bin")).omega._data[0].copy()
        draw(axes[0, col], w_py, f"py/ibpm.py, t={s * 0.01:g}", pts, alpha)
        draw(axes[1, col], w_cpp, f"C++ build/ibpm, t={s * 0.01:g}", pts, alpha)
        if col > 0:
            axes[0, col].set_yticklabels([]); axes[1, col].set_yticklabels([])
        axes[1, col].set_xlabel("x")
    axes[0, 0].set_ylabel("y (py/ibpm.py)"); axes[1, 0].set_ylabel("y (C++)")
    fig.suptitle(f"{label}: vorticity field evolution, py/ibpm.py vs. C++ build/ibpm, "
                 f"Re={re_val}, $\\alpha$={alpha}° (dx=0.02)", fontsize=12)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-VMAX, VMAX), cmap="RdBu_r"),
                 ax=axes, shrink=0.7, label="vorticity $\\omega$ (clipped)")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path} ({ncols} snapshot columns)")


def main():
    # ---- NACA0012 (easy), Re=500, alpha=5 -- freshly run for this folder ----
    make_flow_evolution(
        label="NACA0012", alpha=5.0,
        dat_path=REPO / "SURF_test" / "low_re" / "NACA0012" / "naca0012.dat.txt",
        py_dir=REPO / "SURF_test" / "low_re" / "NACA0012" / "_run_data" / "flowfield",
        cpp_dir=REPO / "SURF_test" / "low_re" / "NACA0012" / "_run_data_cpp" / "flowfield",
        filename_stem="flow", out_path=REPO / "SURF_test" / "low_re" / "NACA0012" / "flow_evolution.png",
        re_val=500,
    )
    # ---- SD7003 (hard), Re=500, alpha=4.60 -- reused from ../airfoils/LSAT-SD7003/4-Re_sweep/ ----
    make_flow_evolution(
        label="SD7003", alpha=4.60,
        dat_path=REPO / "SURF_test" / "airfoils" / "LSAT-SD7003" / "sd7003.dat.txt",
        py_dir=REPO / "SURF_test" / "airfoils" / "LSAT-SD7003" / "4-Re_sweep" / "_run_data" / "Re500",
        cpp_dir=REPO / "SURF_test" / "airfoils" / "LSAT-SD7003" / "4-Re_sweep" / "_run_data_cpp" / "Re500",
        filename_stem="run", out_path=REPO / "SURF_test" / "low_re" / "SD7003" / "flow_evolution.png",
        re_val=500,
    )


if __name__ == "__main__":
    main()
