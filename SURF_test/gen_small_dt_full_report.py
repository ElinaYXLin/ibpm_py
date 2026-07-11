"""
gen_small_dt_full_report.py

Produces the SAME diagram types that live in SD7003/SD8000's
2-c++included/ -- polar_comparison.png, drag_polar.png,
grid_convergence.png, flow_evolution.png -- but from the dt=0.001 sweep
run by run_small_dt_full_sweep.py, into {SD7003,SD8000}/3-small_dt/.

No error-bar whiskers here either (matching 2-c++included/'s convention --
see SD7003/README.md's "About the error bars" section for why); std. dev.
is still written out in full in the accompanying summary.txt.

Usage:  python3 SURF_test/gen_small_dt_full_report.py
"""
import json
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
SURF = REPO / "SURF_test"
sys.path.insert(0, str(SURF))
from parse_uiuc import parse_blocks, nearest_block  # noqa: E402

pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

results = json.loads((SURF / "small_dt_full_results.json").read_text())
DT = 0.001
NSTEPS = 6000

CASES = {"SD7003": dict(Re=61100), "SD8000": dict(Re=60800)}
PY_COLOR, CPP_COLOR = "C0", "C3"


def series(rows, sort_key):
    rows = sorted(rows, key=lambda r: r[sort_key])
    x = np.array([r[sort_key] for r in rows])
    cl = np.array([r["cl_mean"] for r in rows])
    cd = np.array([r["cd_mean"] for r in rows])
    cl_std = np.array([r["cl_std"] for r in rows])
    cd_std = np.array([r["cd_std"] for r in rows])
    return x, cl, cd, cl_std, cd_std


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


def main():
    for name, cfg in CASES.items():
        outdir = SURF / "airfoils" / name / "3-small_dt"
        outdir.mkdir(parents=True, exist_ok=True)

        drg_blocks = parse_blocks(SURF / "airfoils" / name / f"{name}.DRG.txt", "drg")
        exp = nearest_block(drg_blocks, cfg["Re"])
        exp_alpha, exp_cl, exp_cd = np.array(exp["alpha"]), np.array(exp["Cl"]), np.array(exp["Cd"])

        py_alpha, py_cl, py_cd, py_cl_std, py_cd_std = series(results["polar"][name]["py"], "alpha")
        cpp_alpha, cpp_cl, cpp_cd, cpp_cl_std, cpp_cd_std = series(results["polar"][name]["cpp"], "alpha")

        # ---- polar_comparison.png ----
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].plot(exp_alpha, exp_cl, "ko-", label=f"UIUC LSAT experiment (Re={exp['Re']:.0f})", ms=5)
        axes[0].plot(py_alpha, py_cl, marker="s", ls="--", color=PY_COLOR,
                     label=f"py/ibpm.py, native FFTW3, dt={DT}", ms=5)
        axes[0].plot(cpp_alpha, cpp_cl, marker="^", ls="--", color=CPP_COLOR,
                     label=f"C++ build/ibpm, dt={DT}", ms=5)
        axes[0].set_xlabel(r"$\alpha$ (deg)"); axes[0].set_ylabel("$C_l$")
        axes[0].set_title(f"{name}: lift coefficient"); axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

        axes[1].plot(exp_alpha, exp_cd, "ko-", label=f"UIUC LSAT experiment (Re={exp['Re']:.0f})", ms=5)
        axes[1].plot(py_alpha, py_cd, marker="s", ls="--", color=PY_COLOR,
                     label=f"py/ibpm.py, native FFTW3, dt={DT}", ms=5)
        axes[1].plot(cpp_alpha, cpp_cd, marker="^", ls="--", color=CPP_COLOR,
                     label=f"C++ build/ibpm, dt={DT}", ms=5)
        axes[1].set_xlabel(r"$\alpha$ (deg)"); axes[1].set_ylabel("$C_d$")
        axes[1].set_title(f"{name}: drag coefficient"); axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3)
        fig.suptitle(f"{name}: py/ibpm.py (native FFTW3) vs. C++ build/ibpm vs. UIUC LSAT, "
                     f"dt={DT} (10x smaller than 2-c++included/'s dt=0.01)\n"
                     f"averaged over last 60% of a {NSTEPS}-step (t={NSTEPS*DT:g}) run -- shorter "
                     f"window than 2-c++included/'s t=30, see README")
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        fig.savefig(outdir / "polar_comparison.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # ---- drag_polar.png ----
        fig, ax = plt.subplots(figsize=(6, 5.5))
        ax.plot(exp_cd, exp_cl, "ko-", label="UIUC LSAT experiment", ms=5)
        ax.plot(py_cd, py_cl, marker="s", ls="--", color=PY_COLOR, label="py/ibpm.py (native FFTW3)", ms=5)
        ax.plot(cpp_cd, cpp_cl, marker="^", ls="--", color=CPP_COLOR, label="C++ build/ibpm", ms=5)
        ax.set_xlabel("$C_d$"); ax.set_ylabel("$C_l$")
        ax.set_title(f"{name}: drag polar, Re$\\approx${cfg['Re']}, dt={DT}")
        ax.legend(); ax.grid(alpha=0.3)
        fig.savefig(outdir / "drag_polar.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # ---- grid_convergence.png ----
        py_dx, py_conv_cl, py_conv_cd, py_conv_cl_std, py_conv_cd_std = \
            series(results["convergence"][name]["py"], "dx")
        cpp_dx, cpp_conv_cl, cpp_conv_cd, cpp_conv_cl_std, cpp_conv_cd_std = \
            series(results["convergence"][name]["cpp"], "dx")
        conv_alpha = results["convergence"][name]["py"][0]["alpha"]
        exp_row = np.argmin(np.abs(exp_alpha - conv_alpha))
        exp_cl_at_alpha, exp_cd_at_alpha = exp_cl[exp_row], exp_cd[exp_row]
        exp_alpha_matched = exp_alpha[exp_row]

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].plot(py_dx, py_conv_cl, marker="o", ls="-", color=PY_COLOR, label="py/ibpm.py (native FFTW3)")
        axes[0].plot(cpp_dx, cpp_conv_cl, marker="^", ls="-", color=CPP_COLOR, label="C++ build/ibpm")
        axes[0].axhline(exp_cl_at_alpha, color="k", ls="--", label=f"UIUC exp. ($\\alpha$={exp_alpha_matched:.2f}°)")
        axes[0].set_xlabel("grid spacing $dx$"); axes[0].set_ylabel("$C_l$"); axes[0].invert_xaxis()
        axes[0].set_title(f"{name}: $C_l$ grid convergence, $\\alpha$={conv_alpha}°, dt={DT}")
        axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

        axes[1].plot(py_dx, py_conv_cd, marker="o", ls="-", color=PY_COLOR, label="py/ibpm.py (native FFTW3)")
        axes[1].plot(cpp_dx, cpp_conv_cd, marker="^", ls="-", color=CPP_COLOR, label="C++ build/ibpm")
        axes[1].axhline(exp_cd_at_alpha, color="k", ls="--", label=f"UIUC exp. ($\\alpha$={exp_alpha_matched:.2f}°)")
        axes[1].set_xlabel("grid spacing $dx$"); axes[1].set_ylabel("$C_d$"); axes[1].invert_xaxis()
        axes[1].set_title(f"{name}: $C_d$ grid convergence, $\\alpha$={conv_alpha}°, dt={DT}")
        axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3)
        fig.suptitle(f"{name}: grid convergence (dx=0.04,0.02,0.01), py/ibpm.py (native FFTW3) vs. "
                     f"C++ build/ibpm, dt={DT}")
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        fig.savefig(outdir / "grid_convergence.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # ---- summary.txt ----
        with open(outdir / "summary.txt", "w") as f:
            f.write(f"{name} vs. UIUC LSAT experiment (Re_exp={exp['Re']:.0f}, Re_sim={cfg['Re']}), "
                    f"dt={DT} (10x smaller than 2-c++included/'s dt=0.01)\n")
            f.write(f"averaged over last 60% of {NSTEPS} steps (t={NSTEPS*DT:g})\n\n")
            f.write("Polar sweep (dx=0.02, nx=300, ny=150):\n")
            f.write(f"{'alpha':>8} {'Cl_exp':>8} {'Cl_py':>14} {'Cl_cpp':>14} "
                    f"{'Cd_exp':>8} {'Cd_py':>15} {'Cd_cpp':>15}\n")
            for i, a in enumerate(py_alpha):
                row = np.argmin(np.abs(exp_alpha - a))
                j = np.argmin(np.abs(cpp_alpha - a))
                f.write(f"{a:8.2f} {exp_cl[row]:8.3f} {py_cl[i]:7.3f}±{py_cl_std[i]:.3f} "
                        f"{cpp_cl[j]:7.3f}±{cpp_cl_std[j]:.3f} "
                        f"{exp_cd[row]:8.4f} {py_cd[i]:8.4f}±{py_cd_std[i]:.4f} "
                        f"{cpp_cd[j]:8.4f}±{cpp_cd_std[j]:.4f}\n")
            f.write(f"\nGrid convergence at alpha={conv_alpha}° "
                    f"(Cl_exp={exp_cl_at_alpha:.3f}, Cd_exp={exp_cd_at_alpha:.4f}):\n")
            conv_py_sorted = sorted(results["convergence"][name]["py"], key=lambda r: -r["dx"])
            conv_cpp_sorted = sorted(results["convergence"][name]["cpp"], key=lambda r: -r["dx"])
            for cpy, ccpp in zip(conv_py_sorted, conv_cpp_sorted):
                f.write(f"  dx={cpy['dx']:.3f} ({cpy['tag']:6s})  py:  Cl={cpy['cl_mean']:+.4f}"
                        f"±{cpy['cl_std']:.4f}  Cd={cpy['cd_mean']:+.4f}±{cpy['cd_std']:.4f}  "
                        f"({cpy['elapsed']:.0f}s)\n")
                f.write(f"  dx={ccpp['dx']:.3f} ({ccpp['tag']:6s})  cpp: Cl={ccpp['cl_mean']:+.4f}"
                        f"±{ccpp['cl_std']:.4f}  Cd={ccpp['cd_mean']:+.4f}±{ccpp['cd_std']:.4f}  "
                        f"({ccpp['elapsed']:.0f}s)\n")
        print(f"{name}: wrote polar/drag_polar/grid_convergence/summary to {outdir}")

    # ---- flow_evolution.png (separate loop: needs the .dat coordinate files) ----
    DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
    NX, NY = 300, 150
    DX = DOMAIN["length"] / NX
    xs = DOMAIN["xoffset"] + np.arange(1, NX) * DX
    ys = DOMAIN["yoffset"] + np.arange(1, NY) * DX
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    VMAX = 8.0
    LEVELS = np.linspace(-VMAX, VMAX, 41)
    FLOW_CFG = {
        "SD7003": dict(alpha=4.60, dat=SURF / "airfoils" / "SD7003" / "sd7003.dat.txt"),
        "SD8000": dict(alpha=5.36, dat=SURF / "airfoils" / "SD8000" / "sd8000.dat.txt"),
    }

    def draw_field(ax, field, title, pts, alpha_deg):
        ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=LEVELS, cmap="RdBu_r", extend="both")
        th = -np.deg2rad(alpha_deg)
        c, sn = np.cos(th), np.sin(th)
        xc, yc = pts[:, 0] - 0.25, pts[:, 1]
        xr, yr = xc * c - yc * sn + 0.25, xc * sn + yc * c
        ax.fill(xr, yr, color="0.15", zorder=5)
        ax.set_xlim(-2, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
        ax.set_title(title, fontsize=9)

    for name, cfg in FLOW_CFG.items():
        outdir = SURF / "airfoils" / name / "3-small_dt"
        py_dir = SURF / "airfoils" / name / "_run_data_smalldt_full" / "flowfield"
        cpp_dir = SURF / "airfoils" / name / "_run_data_smalldt_full_cpp" / "flowfield"
        pts = load_dat_pts(cfg["dat"])
        steps = [0, 1000, 2000, 3000, 4000, 5000, 6000]
        steps = [s for s in steps if (py_dir / f"run{s:05d}.bin").exists()
                 and (cpp_dir / f"run{s:05d}.bin").exists()]
        if not steps:
            print(f"{name}: no flowfield snapshots yet, skipping flow_evolution.png")
            continue
        ncols = len(steps)
        fig, axes = plt.subplots(2, ncols, figsize=(3.1 * ncols, 6.6))
        axes = np.atleast_2d(axes)
        for col, s in enumerate(steps):
            draw_field(axes[0, col], load_omega(py_dir / f"run{s:05d}.bin"),
                       f"py/ibpm.py, t={s*DT:g}", pts, cfg["alpha"])
            draw_field(axes[1, col], load_omega(cpp_dir / f"run{s:05d}.bin"),
                       f"C++ build/ibpm, t={s*DT:g}", pts, cfg["alpha"])
            if col > 0:
                axes[0, col].set_yticklabels([]); axes[1, col].set_yticklabels([])
            axes[1, col].set_xlabel("x")
        axes[0, 0].set_ylabel("y (py/ibpm.py)"); axes[1, 0].set_ylabel("y (C++)")
        fig.suptitle(f"{name}: vorticity field evolution, py/ibpm.py (native FFTW3) vs. C++ build/ibpm, "
                     f"$\\alpha$={cfg['alpha']}°, dt={DT} (dx=0.02)", fontsize=12)
        fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-VMAX, VMAX), cmap="RdBu_r"),
                     ax=axes, shrink=0.7, label="vorticity $\\omega$ (clipped)")
        fig.savefig(outdir / "flow_evolution.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"{name}: wrote flow_evolution.png ({len(steps)} snapshot pairs)")

    print("done")


if __name__ == "__main__":
    main()
