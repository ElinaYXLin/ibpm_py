"""
gen_airfoil_report_v2.py

Regenerates the SD7003/SD8000 comparison figures with THREE series instead
of two: py/ibpm.py, this repo's C++ reference build (build/ibpm), and the
UIUC LSAT wind-tunnel experiment -- see run_all_airfoils.py (Python) and
run_all_airfoils_cpp.py (C++) for how each dataset was produced.

Difference from the original gen_airfoil_report.py (kept for the record --
its output now lives in SD7003/1-orig, SD8000/1-orig): NO ERROR BARS are drawn on
these figures. The instantaneous Cd/Cl trace this solver produces is
genuinely unsteady at this Re/resolution (vortex shedding, not noise), so
"mean +/- std dev of the trace over the averaging window" is a real,
meaningful quantity -- but a whisker on a polar/convergence plot invites
reading it as "measurement uncertainty" the way an experimentalist's error
bar would be read, which it is not: it's a bandwidth of an oscillation, not
a confidence interval on the mean. To avoid that misreading, the mean is
plotted as a plain marker with no whisker; the underlying std dev is still
written out in full in summary.txt for anyone who wants it.

Usage:
    python3 SURF_test/gen_airfoil_report_v2.py

Output (per case):
    SURF_test/airfoils/LSAT-{SD7003,SD8000}/2-c++included/polar_comparison.png
    SURF_test/airfoils/LSAT-{SD7003,SD8000}/2-c++included/drag_polar.png
    SURF_test/airfoils/LSAT-{SD7003,SD8000}/2-c++included/grid_convergence.png
    SURF_test/airfoils/LSAT-{SD7003,SD8000}/2-c++included/summary.txt
"""
import sys, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
SURF = REPO / "SURF_test"
sys.path.insert(0, str(SURF))
from parse_uiuc import parse_blocks, nearest_block

py_results = json.loads((SURF / "batch_results.json").read_text())
cpp_results = json.loads((SURF / "batch_results_cpp.json").read_text())

CASES = {
    "SD7003": dict(Re=61100),
    "SD8000": dict(Re=60800),
}

PY_COLOR = "C0"
CPP_COLOR = "C3"


def series(results, name, key):
    """Return (alpha_or_dx array, cl array, cd array) for results[key][name],
    sorted by the sweep variable (alpha for polar, dx for convergence)."""
    rows = results[key][name]
    sort_key = "alpha" if key == "polar" else "dx"
    rows = sorted(rows, key=lambda r: r[sort_key])
    x = np.array([r[sort_key] for r in rows])
    cl = np.array([r["cl_mean"] for r in rows])
    cd = np.array([r["cd_mean"] for r in rows])
    cl_std = np.array([r["cl_std"] for r in rows])
    cd_std = np.array([r["cd_std"] for r in rows])
    return x, cl, cd, cl_std, cd_std


for name, cfg in CASES.items():
    outdir = SURF / "airfoils" / f"LSAT-{name}" / "2-c++included"
    outdir.mkdir(parents=True, exist_ok=True)

    drg_blocks = parse_blocks(SURF / "airfoils" / f"LSAT-{name}" / f"{name}.DRG.txt", "drg")
    exp = nearest_block(drg_blocks, cfg["Re"])
    exp_alpha = np.array(exp["alpha"])
    exp_cl = np.array(exp["Cl"])
    exp_cd = np.array(exp["Cd"])

    py_alpha, py_cl, py_cd, py_cl_std, py_cd_std = series(py_results, name, "polar")
    cpp_alpha, cpp_cl, cpp_cd, cpp_cl_std, cpp_cd_std = series(cpp_results, name, "polar")

    # ---- Figure 1: Cl vs alpha, Cd vs alpha ----
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(exp_alpha, exp_cl, "ko-", label=f"UIUC LSAT experiment (Re={exp['Re']:.0f})", ms=5)
    axes[0].plot(py_alpha, py_cl, marker="s", ls="--", color=PY_COLOR,
                 label=f"py/ibpm.py (Re={cfg['Re']}, dx=0.02)", ms=5)
    axes[0].plot(cpp_alpha, cpp_cl, marker="^", ls="--", color=CPP_COLOR,
                 label=f"C++ build/ibpm (Re={cfg['Re']}, dx=0.02)", ms=5)
    axes[0].set_xlabel(r"$\alpha$ (deg)")
    axes[0].set_ylabel("$C_l$")
    axes[0].set_title(f"{name}: lift coefficient")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].plot(exp_alpha, exp_cd, "ko-", label=f"UIUC LSAT experiment (Re={exp['Re']:.0f})", ms=5)
    axes[1].plot(py_alpha, py_cd, marker="s", ls="--", color=PY_COLOR,
                 label=f"py/ibpm.py (Re={cfg['Re']}, dx=0.02)", ms=5)
    axes[1].plot(cpp_alpha, cpp_cd, marker="^", ls="--", color=CPP_COLOR,
                 label=f"C++ build/ibpm (Re={cfg['Re']}, dx=0.02)", ms=5)
    axes[1].set_xlabel(r"$\alpha$ (deg)")
    axes[1].set_ylabel("$C_d$")
    axes[1].set_title(f"{name}: drag coefficient")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"{name}: py/ibpm.py vs. C++ build/ibpm vs. UIUC LSAT wind-tunnel data\n"
                 f"(markers = time-averaged mean only; see summary.txt for the underlying "
                 f"unsteady-trace std. dev.)")
    fig.tight_layout()
    fig.savefig(outdir / "polar_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 2: Cl-Cd polar (drag polar) ----
    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.plot(exp_cd, exp_cl, "ko-", label="UIUC LSAT experiment", ms=5)
    ax.plot(py_cd, py_cl, marker="s", ls="--", color=PY_COLOR, label="py/ibpm.py", ms=5)
    ax.plot(cpp_cd, cpp_cl, marker="^", ls="--", color=CPP_COLOR, label="C++ build/ibpm", ms=5)
    ax.set_xlabel("$C_d$")
    ax.set_ylabel("$C_l$")
    ax.set_title(f"{name}: drag polar, Re$\\approx${cfg['Re']}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(outdir / "drag_polar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 3: grid convergence ----
    py_dx, py_conv_cl, py_conv_cd, py_conv_cl_std, py_conv_cd_std = series(py_results, name, "convergence")
    cpp_dx, cpp_conv_cl, cpp_conv_cd, cpp_conv_cl_std, cpp_conv_cd_std = series(cpp_results, name, "convergence")
    conv_alpha = py_results["convergence"][name][0]["alpha"]
    exp_row = np.argmin(np.abs(exp_alpha - conv_alpha))
    exp_cl_at_alpha = exp_cl[exp_row]
    exp_cd_at_alpha = exp_cd[exp_row]
    exp_alpha_matched = exp_alpha[exp_row]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(py_dx, py_conv_cl, marker="o", ls="-", color=PY_COLOR, label="py/ibpm.py")
    axes[0].plot(cpp_dx, cpp_conv_cl, marker="^", ls="-", color=CPP_COLOR, label="C++ build/ibpm")
    axes[0].axhline(exp_cl_at_alpha, color="k", ls="--",
                     label=f"UIUC exp. ($\\alpha$={exp_alpha_matched:.2f}°)")
    axes[0].set_xlabel("grid spacing $dx$")
    axes[0].set_ylabel("$C_l$")
    axes[0].invert_xaxis()
    axes[0].set_title(f"{name}: $C_l$ grid convergence, $\\alpha$={conv_alpha}°")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].plot(py_dx, py_conv_cd, marker="o", ls="-", color=PY_COLOR, label="py/ibpm.py")
    axes[1].plot(cpp_dx, cpp_conv_cd, marker="^", ls="-", color=CPP_COLOR, label="C++ build/ibpm")
    axes[1].axhline(exp_cd_at_alpha, color="k", ls="--",
                     label=f"UIUC exp. ($\\alpha$={exp_alpha_matched:.2f}°)")
    axes[1].set_xlabel("grid spacing $dx$")
    axes[1].set_ylabel("$C_d$")
    axes[1].invert_xaxis()
    axes[1].set_title(f"{name}: $C_d$ grid convergence, $\\alpha$={conv_alpha}°")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"{name}: grid convergence study (dx = 0.04, 0.02, 0.01), py/ibpm.py vs. C++ "
                 f"build/ibpm\n(markers = time-averaged mean only; see summary.txt for the "
                 f"underlying unsteady-trace std. dev.)")
    fig.tight_layout()
    fig.savefig(outdir / "grid_convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- summary text ----
    with open(outdir / "summary.txt", "w") as f:
        f.write(f"{name} vs. UIUC LSAT experiment (Re_exp={exp['Re']:.0f}, Re_sim={cfg['Re']})\n")
        f.write("py/ibpm.py vs. C++ build/ibpm vs. UIUC LSAT experiment\n\n")
        f.write(
            "NOTE on the +/- figures below: these are the std. dev. of the instantaneous\n"
            "Cl/Cd trace over the last 60% of each run (this flow genuinely vortex-sheds at\n"
            "this Re/resolution, so this is a real oscillation amplitude, not measurement\n"
            "noise). The comparison figures in this folder plot the mean only, with no\n"
            "error-bar whiskers, to avoid that oscillation-amplitude number being misread as\n"
            "a confidence interval the way an experimental error bar would be -- see this\n"
            "script's module docstring for the full explanation.\n\n"
        )
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
                f"(experiment closest alpha={exp_alpha_matched:.2f}°, "
                f"Cl_exp={exp_cl_at_alpha:.3f}, Cd_exp={exp_cd_at_alpha:.4f}):\n")
        conv_py_rows = sorted(py_results["convergence"][name], key=lambda r: -r["dx"])
        conv_cpp_rows = sorted(cpp_results["convergence"][name], key=lambda r: -r["dx"])
        for cpy, ccpp in zip(conv_py_rows, conv_cpp_rows):
            f.write(f"  dx={cpy['dx']:.3f} ({cpy['tag']:6s})  "
                    f"py:  Cl={cpy['cl_mean']:+.4f}±{cpy['cl_std']:.4f}  "
                    f"Cd={cpy['cd_mean']:+.4f}±{cpy['cd_std']:.4f}  ({cpy['elapsed']:.0f}s)\n")
            f.write(f"  dx={ccpp['dx']:.3f} ({ccpp['tag']:6s})  "
                    f"cpp: Cl={ccpp['cl_mean']:+.4f}±{ccpp['cl_std']:.4f}  "
                    f"Cd={ccpp['cd_mean']:+.4f}±{ccpp['cd_std']:.4f}  ({ccpp['elapsed']:.0f}s)\n")

    print(f"{name}: wrote figures/summary to {outdir}")

print("done")
