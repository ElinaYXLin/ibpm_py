"""
gen_cylinder_re_sweep_figs.py

Figures for run_cylinder_re_sweep.py's Re sweep (E1+E2/E3 combined):
- re_sweep_comparison.png: grid of (Re rows) x (vorticity | velocity
  magnitude columns) at t=30, the fully-evolved snapshot of each short
  run. Answers the core visual question directly: does raising Re from
  the clean Re=100 baseline toward SD7003/SD8000's Re~40-61k make the
  cylinder wake go from clean to speckled, and does it look different in
  velocity vs. vorticity?
- domain_rms_vs_Re.png: quantitative companion, domain-RMS |omega| vs Re,
  pinning down the transition numerically (same style as
  airfoils/SD7003/3-ngrid=3/instability_diagnostic.png panel A/C).

Usage: python3 SURF_test/vortall/gen_cylinder_re_sweep_figs.py
Output: SURF_test/vortall/2-Re_sweep/re_sweep_comparison.png,
        SURF_test/vortall/2-Re_sweep/domain_rms_vs_Re.png
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

OUTDIR = REPO / "SURF_test" / "vortall" / "2-Re_sweep"
RUNDIR = OUTDIR / "_run_data_cpp"
RE_VALUES = [100, 500, 1000, 3000, 10000]
# Re>=1000 reran at dt=0.005/nsteps=6000 (see run_cylinder_re_sweep.py's
# DT_OVERRIDE) after diverging at dt=0.02; final step differs accordingly.
FINAL_STEP_OVERRIDE = {1000: 6000, 3000: 6000, 10000: 6000}
DEFAULT_FINAL_STEP = 1500  # t=30 at dt=0.02
DOMAIN = dict(length=9, xoffset=-1, yoffset=-2)
NX, NY = 450, 200
DX = DOMAIN["length"] / NX
xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")


def load_state(path):
    return State(filename=str(path))


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


def draw(ax, field, title, vmax, cmap):
    ax.contourf(X, Y, np.clip(field, -vmax if cmap == "RdBu_r" else 0, vmax),
                levels=41, cmap=cmap, extend="both")
    circ = plt.Circle((0, 0), 0.5, color="0.15", zorder=5)
    ax.add_patch(circ)
    ax.set_xlim(-1, 8); ax.set_ylim(-2, 2); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def main():
    available = [Re for Re in RE_VALUES
                 if (RUNDIR / f"Re{Re}" / f"cyl{FINAL_STEP_OVERRIDE.get(Re, DEFAULT_FINAL_STEP):05d}.bin").exists()]
    if not available:
        print("no completed Re-sweep snapshots found yet")
        return

    fig, axes = plt.subplots(len(available), 2, figsize=(9, 3.4 * len(available)))
    axes = np.atleast_2d(axes)
    re_list, rms_list, max_list = [], [], []
    for row, Re in enumerate(available):
        final_step = FINAL_STEP_OVERRIDE.get(Re, DEFAULT_FINAL_STEP)
        s = load_state(RUNDIR / f"Re{Re}" / f"cyl{final_step:05d}.bin")
        w = vorticity(s)
        vmag = velocity_mag(s)
        draw(axes[row, 0], w, f"Re={Re}: vorticity, t=30", vmax=8.0, cmap="RdBu_r")
        draw(axes[row, 1], vmag, f"Re={Re}: |velocity|, t=30", vmax=2.0, cmap="viridis")
        re_list.append(Re)
        rms_list.append(rms(w))
        max_list.append(float(np.max(np.abs(w))))
    axes[0, 0].set_title(axes[0, 0].get_title() + "\n(vorticity)")
    fig.suptitle("Cylinder Re sweep: does raising Re from the clean Re=100 baseline\n"
                 "toward the airfoils' Re~40-61k make the wake go from clean to speckled?", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUTDIR / "re_sweep_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 're_sweep_comparison.png'} ({len(available)} Re values)")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(re_list, rms_list, "o-", color="C0")
    ax[0].set_xscale("log"); ax[0].set_xlabel("Re"); ax[0].set_ylabel(r"domain-RMS $|\omega|$, t=30")
    ax[0].set_title("Domain-RMS vorticity vs. Re"); ax[0].grid(alpha=0.3)
    ax[1].plot(re_list, max_list, "s-", color="C3")
    ax[1].set_xscale("log"); ax[1].set_xlabel("Re"); ax[1].set_ylabel(r"max $|\omega|$, t=30")
    ax[1].set_title("Peak vorticity vs. Re"); ax[1].grid(alpha=0.3)
    fig.suptitle("Cylinder Re sweep: quantitative onset of broadband noise", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUTDIR / "domain_rms_vs_Re.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 'domain_rms_vs_Re.png'}")
    print("\nRe, RMS|omega|, max|omega|:")
    for Re, r, m in zip(re_list, rms_list, max_list):
        print(f"  Re={Re}: RMS={r:.4f}  max={m:.4f}")


if __name__ == "__main__":
    main()
