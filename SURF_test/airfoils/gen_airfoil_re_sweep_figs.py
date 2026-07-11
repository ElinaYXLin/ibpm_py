"""
gen_airfoil_re_sweep_figs.py

Figures for run_airfoil_re_sweep.py / run_airfoil_re_sweep_py.py's Re
sweep: grid of (Re rows) x (py vorticity | C++ vorticity columns) at
t=30 -- the core Python-vs-C++ fidelity check for every swept Re, not
just the pre-existing baseline. Plus quantitative companions: domain-RMS/
max vorticity vs Re for both implementations overlaid, and a
py-vs-cpp agreement table (RMS of the pointwise difference, relative to
the field's own RMS).

Usage: python3 SURF_test/airfoils/gen_airfoil_re_sweep_figs.py <SD7003|SD8000>
Output: SURF_test/airfoils/<name>/4-Re_sweep/re_sweep_comparison.png,
        SURF_test/airfoils/<name>/4-Re_sweep/domain_rms_vs_Re.png,
        SURF_test/airfoils/<name>/4-Re_sweep/fidelity_summary.txt
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
PY_COLOR, CPP_COLOR = "C0", "C3"


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


def vorticity(path):
    return State(filename=str(path)).omega._data[0].copy()


def rms(f):
    return float(np.sqrt(np.nanmean(f.astype(np.float64) ** 2)))


def draw(ax, field, title, pts, alpha_deg, vmax=8.0):
    ax.contourf(X, Y, np.clip(field, -vmax, vmax), levels=41, cmap="RdBu_r", extend="both")
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
    py_dir = outdir / "_run_data"
    cpp_dir = outdir / "_run_data_cpp"
    pts = load_dat_pts(cfg["dat"])

    available = [Re for Re in RE_VALUES
                 if (py_dir / f"Re{Re}" / f"run{FINAL_STEP:05d}.bin").exists()
                 and (cpp_dir / f"Re{Re}" / f"run{FINAL_STEP:05d}.bin").exists()]
    if not available:
        print(f"{name}: no matching Python+C++ Re-sweep snapshots found")
        return

    fig, axes = plt.subplots(len(available), 2, figsize=(9, 2.6 * len(available)))
    axes = np.atleast_2d(axes)
    re_list, py_rms_list, cpp_rms_list, py_max_list, cpp_max_list, rel_diff_list = [], [], [], [], [], []
    for row, Re in enumerate(available):
        w_py = vorticity(py_dir / f"Re{Re}" / f"run{FINAL_STEP:05d}.bin")
        w_cpp = vorticity(cpp_dir / f"Re{Re}" / f"run{FINAL_STEP:05d}.bin")
        draw(axes[row, 0], w_py, f"Re={Re}: py/ibpm.py, t=30", pts, cfg["alpha"])
        draw(axes[row, 1], w_cpp, f"Re={Re}: C++ build/ibpm, t=30", pts, cfg["alpha"])
        re_list.append(Re)
        py_rms_list.append(rms(w_py))
        cpp_rms_list.append(rms(w_cpp))
        py_max_list.append(float(np.max(np.abs(w_py))))
        cpp_max_list.append(float(np.max(np.abs(w_cpp))))
        diff_rms = rms(w_py - w_cpp)
        rel_diff_list.append(diff_rms / max(rms(w_cpp), 1e-12))
    axes[0, 0].set_ylabel("y")
    fig.suptitle(f"{name} Re sweep: py/ibpm.py vs. C++ build/ibpm fidelity, t=30\n"
                 f"(does lowering Re from ~40-61k toward the cylinder's clean Re=100 baseline "
                 f"clean up the wake, in BOTH implementations?)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / "re_sweep_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir / 're_sweep_comparison.png'} ({len(available)} Re values)")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(re_list, py_rms_list, "o-", color=PY_COLOR, label="py/ibpm.py")
    ax[0].plot(re_list, cpp_rms_list, "^--", color=CPP_COLOR, label="C++ build/ibpm")
    ax[0].set_xscale("log"); ax[0].set_xlabel("Re"); ax[0].set_ylabel(r"domain-RMS $|\omega|$, t=30")
    ax[0].set_title(f"{name}: domain-RMS vorticity vs. Re"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(re_list, py_max_list, "o-", color=PY_COLOR, label="py/ibpm.py")
    ax[1].plot(re_list, cpp_max_list, "^--", color=CPP_COLOR, label="C++ build/ibpm")
    ax[1].set_xscale("log"); ax[1].set_xlabel("Re"); ax[1].set_ylabel(r"max $|\omega|$, t=30")
    ax[1].set_title(f"{name}: peak vorticity vs. Re"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.suptitle(f"{name} Re sweep: quantitative trend vs. Re, both implementations overlaid", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(outdir / "domain_rms_vs_Re.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir / 'domain_rms_vs_Re.png'}")

    with open(outdir / "fidelity_summary.txt", "w") as f:
        f.write(f"{name} Re sweep: py/ibpm.py vs. C++ build/ibpm fidelity at t=30\n\n")
        f.write(f"{'Re':>8} {'RMS_py':>10} {'RMS_cpp':>10} {'max_py':>10} {'max_cpp':>10} {'RMS(py-cpp)/RMS_cpp':>22}\n")
        for Re, rp, rc, mp, mc, rd in zip(re_list, py_rms_list, cpp_rms_list, py_max_list, cpp_max_list, rel_diff_list):
            f.write(f"{Re:8d} {rp:10.4f} {rc:10.4f} {mp:10.3f} {mc:10.3f} {rd:22.4%}\n")
    print(f"wrote {outdir / 'fidelity_summary.txt'}")
    print(f"\n{name} Re, RMS_py, RMS_cpp, relative diff:")
    for Re, rp, rc, rd in zip(re_list, py_rms_list, cpp_rms_list, rel_diff_list):
        print(f"  Re={Re}: RMS_py={rp:.4f}  RMS_cpp={rc:.4f}  rel_diff={rd:.2%}")


if __name__ == "__main__":
    main()
