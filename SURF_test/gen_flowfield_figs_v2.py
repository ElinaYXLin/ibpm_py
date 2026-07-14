"""
gen_flowfield_figs_v2.py

C++-vs-Python version of gen_flowfield_figs.py's vorticity evolution
montage (original output now lives in SD7003/1-orig, SD8000/1-orig): same 7
snapshots (t=0..30), but with a C++ build/ibpm row added underneath the
py/ibpm.py row for a direct visual comparison, both read from their own
_run_data{,_cpp}/flowfield/ restart files (see run_flowfield.py /
run_flowfield_cpp.py).

Usage: python3 SURF_test/gen_flowfield_figs_v2.py
Output: SURF_test/airfoils/LSAT-{SD7003,SD8000}/2-c++included/flow_evolution.png
"""
import sys, types, pathlib
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

CASES = {
    "SD7003": dict(alpha=4.60, Re=61100, dat=REPO / "SURF_test" / "airfoils" / "LSAT-SD7003" / "sd7003.dat.txt"),
    "SD8000": dict(alpha=5.36, Re=60800, dat=REPO / "SURF_test" / "airfoils" / "LSAT-SD8000" / "sd8000.dat.txt"),
}


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
    s = State(filename=str(path))
    return s.omega._data[0].copy()


VMAX = 8.0
LEVELS = np.linspace(-VMAX, VMAX, 41)


def draw_field(ax, field, title, airfoil_pts, alpha_deg):
    # NOTE(fix): py/ibpm.py's `-alpha` only tilts the free-stream (BaseFlow);
    # it is never applied to the geometry, so `field` is in the frame where
    # the body sits at its raw, UNROTATED orientation. This used to rotate
    # only the drawn outline by -alpha_deg without rotating `field` to
    # match -- inconsistent, and equivalent to plotting the body about one
    # grid cell away from where it truly sits (at this dx, alpha). Plotting
    # the raw points directly matches what was actually solved. `alpha_deg`
    # is kept for call-site compatibility but no longer used for a rotation.
    im = ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=LEVELS, cmap="RdBu_r", extend="both")
    ax.fill(airfoil_pts[:, 0], airfoil_pts[:, 1], color="0.15", zorder=5)
    ax.set_xlim(-2, 4)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    return im


for name, cfg in CASES.items():
    outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "2-c++included"
    outdir.mkdir(parents=True, exist_ok=True)
    py_rundir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data" / "flowfield"
    cpp_rundir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_cpp" / "flowfield"
    airfoil_pts = load_dat_pts(cfg["dat"])

    steps = [0, 500, 1000, 1500, 2000, 2500, 3000]
    py_steps = [s for s in steps if (py_rundir / f"flow{s:05d}.bin").exists()]
    cpp_steps = [s for s in steps if (cpp_rundir / f"flow{s:05d}.bin").exists()]
    steps = [s for s in steps if s in py_steps and s in cpp_steps]
    if not steps:
        print(f"{name}: no matching Python+C++ flowfield snapshots found, skipping")
        continue

    ncols = len(steps)
    fig, axes = plt.subplots(2, ncols, figsize=(3.1 * ncols, 6.6))
    for col, s in enumerate(steps):
        py_field = load_omega(py_rundir / f"flow{s:05d}.bin")
        cpp_field = load_omega(cpp_rundir / f"flow{s:05d}.bin")
        draw_field(axes[0, col], py_field, f"py/ibpm.py, t={s * 0.01:g}", airfoil_pts, cfg["alpha"])
        draw_field(axes[1, col], cpp_field, f"C++ build/ibpm, t={s * 0.01:g}", airfoil_pts, cfg["alpha"])
        if col > 0:
            axes[0, col].set_yticklabels([])
            axes[1, col].set_yticklabels([])
        axes[1, col].set_xlabel("x")
    axes[0, 0].set_ylabel("y (py/ibpm.py)")
    axes[1, 0].set_ylabel("y (C++)")

    fig.suptitle(f"{name}: vorticity field evolution, py/ibpm.py vs. C++ build/ibpm, "
                 f"Re={cfg['Re']} $\\alpha$={cfg['alpha']}° (dx=0.02)", fontsize=12)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-VMAX, VMAX), cmap="RdBu_r"),
                 ax=axes, shrink=0.7, label="vorticity $\\omega$ (clipped)")
    fig.savefig(outdir / "flow_evolution.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"{name}: wrote {outdir / 'flow_evolution.png'}")

print("done")
