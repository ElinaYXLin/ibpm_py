"""
gen_cylinder_re_sweep_figs.py

Figures for run_cylinder_re_sweep.py / run_cylinder_re_sweep_py.py's Re
sweep: grid of (Re rows) x (py vorticity | C++ vorticity columns) at
t=30 -- Python-vs-C++ fidelity check at every swept Re. Plus quantitative
companions: domain-RMS/max vorticity vs Re for both implementations
overlaid, and a py-vs-cpp agreement table.

Usage: python3 SURF_test/vortall/gen_cylinder_re_sweep_figs.py
Output: SURF_test/vortall/2-Re_sweep/re_sweep_comparison.png,
        SURF_test/vortall/2-Re_sweep/domain_rms_vs_Re.png,
        SURF_test/vortall/2-Re_sweep/fidelity_summary.txt
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

OUTDIR = REPO / "SURF_test" / "vortall" / "2-Re_sweep"
PY_DIR = OUTDIR / "_run_data"
CPP_DIR = OUTDIR / "_run_data_cpp"
RE_VALUES = [100, 500, 1000, 3000, 10000]
FINAL_STEP_OVERRIDE = {1000: 6000, 3000: 6000, 10000: 6000}
DEFAULT_FINAL_STEP = 1500
DOMAIN = dict(length=9, xoffset=-1, yoffset=-2)
NX, NY = 450, 200
DX = DOMAIN["length"] / NX
xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")
PY_COLOR, CPP_COLOR = "C0", "C3"


def vorticity(path):
    return State(filename=str(path)).omega._data[0].copy()


def rms(f):
    return float(np.sqrt(np.nanmean(f.astype(np.float64) ** 2)))


def draw(ax, field, title, vmax=8.0):
    ax.contourf(X, Y, np.clip(field, -vmax, vmax), levels=41, cmap="RdBu_r", extend="both")
    circ = plt.Circle((0, 0), 0.5, color="0.15", zorder=5)
    ax.add_patch(circ)
    ax.set_xlim(-1, 8); ax.set_ylim(-2, 2); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)


def main():
    available = []
    for Re in RE_VALUES:
        step = FINAL_STEP_OVERRIDE.get(Re, DEFAULT_FINAL_STEP)
        if (PY_DIR / f"Re{Re}" / f"cyl{step:05d}.bin").exists() and (CPP_DIR / f"Re{Re}" / f"cyl{step:05d}.bin").exists():
            available.append(Re)
    if not available:
        print("no matching Python+C++ Re-sweep snapshots found yet")
        return

    fig, axes = plt.subplots(len(available), 2, figsize=(9, 3.4 * len(available)))
    axes = np.atleast_2d(axes)
    re_list, py_rms_list, cpp_rms_list, py_max_list, cpp_max_list, rel_diff_list = [], [], [], [], [], []
    for row, Re in enumerate(available):
        step = FINAL_STEP_OVERRIDE.get(Re, DEFAULT_FINAL_STEP)
        w_py = vorticity(PY_DIR / f"Re{Re}" / f"cyl{step:05d}.bin")
        w_cpp = vorticity(CPP_DIR / f"Re{Re}" / f"cyl{step:05d}.bin")
        draw(axes[row, 0], w_py, f"Re={Re}: py/ibpm.py, t=30")
        draw(axes[row, 1], w_cpp, f"Re={Re}: C++ build/ibpm, t=30")
        re_list.append(Re)
        py_rms_list.append(rms(w_py)); cpp_rms_list.append(rms(w_cpp))
        py_max_list.append(float(np.max(np.abs(w_py)))); cpp_max_list.append(float(np.max(np.abs(w_cpp))))
        rel_diff_list.append(rms(w_py - w_cpp) / max(rms(w_cpp), 1e-12))
    fig.suptitle("Cylinder Re sweep: py/ibpm.py vs. C++ build/ibpm fidelity, t=30\n"
                 "(does raising Re from the clean Re=100 baseline toward the airfoils' Re~40-61k\n"
                 "make the wake go from clean to speckled, in BOTH implementations?)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(OUTDIR / "re_sweep_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 're_sweep_comparison.png'} ({len(available)} Re values)")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(re_list, py_rms_list, "o-", color=PY_COLOR, label="py/ibpm.py")
    ax[0].plot(re_list, cpp_rms_list, "^--", color=CPP_COLOR, label="C++ build/ibpm")
    ax[0].set_xscale("log"); ax[0].set_xlabel("Re"); ax[0].set_ylabel(r"domain-RMS $|\omega|$, t=30")
    ax[0].set_title("Domain-RMS vorticity vs. Re"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(re_list, py_max_list, "o-", color=PY_COLOR, label="py/ibpm.py")
    ax[1].plot(re_list, cpp_max_list, "^--", color=CPP_COLOR, label="C++ build/ibpm")
    ax[1].set_xscale("log"); ax[1].set_xlabel("Re"); ax[1].set_ylabel(r"max $|\omega|$, t=30")
    ax[1].set_title("Peak vorticity vs. Re"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.suptitle("Cylinder Re sweep: quantitative trend vs. Re, both implementations overlaid", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUTDIR / "domain_rms_vs_Re.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 'domain_rms_vs_Re.png'}")

    with open(OUTDIR / "fidelity_summary.txt", "w") as f:
        f.write("Cylinder Re sweep: py/ibpm.py vs. C++ build/ibpm fidelity at t=30\n\n")
        f.write(f"{'Re':>8} {'RMS_py':>10} {'RMS_cpp':>10} {'max_py':>10} {'max_cpp':>10} {'RMS(py-cpp)/RMS_cpp':>22}\n")
        for Re, rp, rc, mp, mc, rd in zip(re_list, py_rms_list, cpp_rms_list, py_max_list, cpp_max_list, rel_diff_list):
            f.write(f"{Re:8d} {rp:10.4f} {rc:10.4f} {mp:10.3f} {mc:10.3f} {rd:22.4%}\n")
    print(f"wrote {OUTDIR / 'fidelity_summary.txt'}")
    print("\nRe, RMS_py, RMS_cpp, relative diff:")
    for Re, rp, rc, rd in zip(re_list, py_rms_list, cpp_rms_list, rel_diff_list):
        print(f"  Re={Re}: RMS_py={rp:.4f}  RMS_cpp={rc:.4f}  rel_diff={rd:.2%}")


if __name__ == "__main__":
    main()
