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
    "SD7003": dict(alpha=4.60, Re=61100, dat=REPO / "SURF_test" / "SD7003" / "sd7003.dat.txt"),
    "SD8000": dict(alpha=5.36, Re=60800, dat=REPO / "SURF_test" / "SD8000" / "sd8000.dat.txt"),
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
    im = ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=LEVELS, cmap="RdBu_r", extend="both")
    # rotate airfoil outline into the lab frame the same way the solver rotates the freestream
    # (body fixed at alpha=0, flow rotated by +alpha -> equivalent lab-frame body rotation is -alpha)
    th = -np.deg2rad(alpha_deg)
    c, sn = np.cos(th), np.sin(th)
    xc, yc = airfoil_pts[:, 0] - 0.25, airfoil_pts[:, 1]
    xr = xc * c - yc * sn + 0.25
    yr = xc * sn + yc * c
    ax.fill(xr, yr, color="0.15", zorder=5)
    ax.set_xlim(-2, 4)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x")
    return im


for name, cfg in CASES.items():
    outdir = REPO / "results" / name
    rundir = outdir / "_run_data" / "flowfield"
    bins = sorted(rundir.glob("flow?????.bin"))
    airfoil_pts = load_dat_pts(cfg["dat"])

    # Evolution montage
    steps = [0, 500, 1000, 1500, 2000, 2500, 3000]
    steps = [s for s in steps if (rundir / f"flow{s:05d}.bin").exists()]
    ncols = 3
    nrows = (len(steps) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.3 * nrows))
    axes = np.atleast_1d(axes).flatten()
    for ax, s in zip(axes, steps):
        field = load_omega(rundir / f"flow{s:05d}.bin")
        draw_field(ax, field, f"t = {s * 0.01:g}", airfoil_pts, cfg["alpha"])
    for ax in axes[len(steps):]:
        ax.axis("off")
    fig.suptitle(f"{name}: vorticity field evolution, Re={cfg['Re']} "
                 f"$\\alpha$={cfg['alpha']}° (py/ibpm.py, dx=0.02)", fontsize=11)
    im = axes[0].collections[0] if axes[0].collections else None
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-VMAX, VMAX), cmap="RdBu_r"),
                 ax=axes[:len(steps)], shrink=0.8, label="vorticity $\\omega$ (clipped)")
    fig.savefig(outdir / "flow_evolution.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"{name}: wrote {outdir / 'flow_evolution.png'}")

print("done")
