import sys, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
SURF = REPO / "SURF_test"
sys.path.insert(0, str(SURF))
from parse_uiuc import parse_blocks, nearest_block
results = json.loads((SURF / "batch_results.json").read_text())

CASES = {
    "SD7003": dict(Re=61100),
    "SD8000": dict(Re=60800),
}

for name, cfg in CASES.items():
    outdir = REPO / "results" / name
    outdir.mkdir(exist_ok=True)

    drg_blocks = parse_blocks(SURF / name / f"{name}.DRG.txt", "drg")
    exp = nearest_block(drg_blocks, cfg["Re"])
    exp_alpha = np.array(exp["alpha"])
    exp_cl = np.array(exp["Cl"])
    exp_cd = np.array(exp["Cd"])

    polar = results["polar"][name]
    sim_alpha = np.array([p["alpha"] for p in polar])
    sim_cl = np.array([p["cl_mean"] for p in polar])
    sim_cl_std = np.array([p["cl_std"] for p in polar])
    sim_cd = np.array([p["cd_mean"] for p in polar])
    sim_cd_std = np.array([p["cd_std"] for p in polar])

    # ---- Figure 1: Cl vs alpha, Cd vs alpha ----
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(exp_alpha, exp_cl, "ko-", label=f"UIUC LSAT experiment (Re={exp['Re']:.0f})", ms=5)
    axes[0].errorbar(sim_alpha, sim_cl, yerr=sim_cl_std, fmt="C0s--",
                      label=f"py/ibpm.py (Re={cfg['Re']}, dx=0.02)", capsize=3, ms=5)
    axes[0].set_xlabel(r"$\alpha$ (deg)")
    axes[0].set_ylabel("$C_l$")
    axes[0].set_title(f"{name}: lift coefficient")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].plot(exp_alpha, exp_cd, "ko-", label=f"UIUC LSAT experiment (Re={exp['Re']:.0f})", ms=5)
    axes[1].errorbar(sim_alpha, sim_cd, yerr=sim_cd_std, fmt="C1s--",
                      label=f"py/ibpm.py (Re={cfg['Re']}, dx=0.02)", capsize=3, ms=5)
    axes[1].set_xlabel(r"$\alpha$ (deg)")
    axes[1].set_ylabel("$C_d$")
    axes[1].set_title(f"{name}: drag coefficient")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"{name}: py/ibpm.py polar vs. UIUC LSAT wind-tunnel data "
                 f"(error bars = ±1 std. dev. of unsteady force trace)")
    fig.tight_layout()
    fig.savefig(outdir / "polar_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 2: Cl-Cd polar (drag polar) ----
    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.plot(exp_cd, exp_cl, "ko-", label="UIUC LSAT experiment", ms=5)
    ax.errorbar(sim_cd, sim_cl, xerr=sim_cd_std, yerr=sim_cl_std, fmt="C0s--",
                label="py/ibpm.py", capsize=3, ms=5)
    ax.set_xlabel("$C_d$")
    ax.set_ylabel("$C_l$")
    ax.set_title(f"{name}: drag polar, Re$\\approx${cfg['Re']}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(outdir / "drag_polar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 3: grid convergence ----
    conv = results["convergence"][name]
    conv_dx = np.array([c["dx"] for c in conv])
    conv_cl = np.array([c["cl_mean"] for c in conv])
    conv_cl_std = np.array([c["cl_std"] for c in conv])
    conv_cd = np.array([c["cd_mean"] for c in conv])
    conv_cd_std = np.array([c["cd_std"] for c in conv])
    conv_alpha = conv[0]["alpha"]
    exp_row = np.argmin(np.abs(exp_alpha - conv_alpha))
    exp_cl_at_alpha = exp_cl[exp_row]
    exp_cd_at_alpha = exp_cd[exp_row]
    exp_alpha_matched = exp_alpha[exp_row]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].errorbar(conv_dx, conv_cl, yerr=conv_cl_std, fmt="C0o-", capsize=3,
                      label="py/ibpm.py")
    axes[0].axhline(exp_cl_at_alpha, color="k", ls="--",
                     label=f"UIUC exp. ($\\alpha$={exp_alpha_matched:.2f}°)")
    axes[0].set_xlabel("grid spacing $dx$")
    axes[0].set_ylabel("$C_l$")
    axes[0].invert_xaxis()
    axes[0].set_title(f"{name}: $C_l$ grid convergence, $\\alpha$={conv_alpha}°")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].errorbar(conv_dx, conv_cd, yerr=conv_cd_std, fmt="C1o-", capsize=3,
                      label="py/ibpm.py")
    axes[1].axhline(exp_cd_at_alpha, color="k", ls="--",
                     label=f"UIUC exp. ($\\alpha$={exp_alpha_matched:.2f}°)")
    axes[1].set_xlabel("grid spacing $dx$")
    axes[1].set_ylabel("$C_d$")
    axes[1].invert_xaxis()
    axes[1].set_title(f"{name}: $C_d$ grid convergence, $\\alpha$={conv_alpha}°")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"{name}: grid convergence study (dx = 0.04, 0.02, 0.01; error bars = "
                 f"±1 std. dev. of unsteady force trace)")
    fig.tight_layout()
    fig.savefig(outdir / "grid_convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- summary text ----
    with open(outdir / "summary.txt", "w") as f:
        f.write(f"{name} vs. UIUC LSAT experiment (Re_exp={exp['Re']:.0f}, Re_sim={cfg['Re']})\n\n")
        f.write("Polar sweep (dx=0.02, nx=300, ny=150):\n")
        f.write(f"{'alpha':>8} {'Cl_exp':>8} {'Cl_sim':>10} {'Cd_exp':>8} {'Cd_sim':>10}\n")
        for a, cl_s, cl_sd, cd_s, cd_sd in zip(sim_alpha, sim_cl, sim_cl_std, sim_cd, sim_cd_std):
            row = np.argmin(np.abs(exp_alpha - a))
            f.write(f"{a:8.2f} {exp_cl[row]:8.3f} {cl_s:6.3f}±{cl_sd:.3f} "
                    f"{exp_cd[row]:8.4f} {cd_s:7.4f}±{cd_sd:.4f}\n")
        f.write(f"\nGrid convergence at alpha={conv_alpha}° "
                f"(experiment closest alpha={exp_alpha_matched:.2f}°, "
                f"Cl_exp={exp_cl_at_alpha:.3f}, Cd_exp={exp_cd_at_alpha:.4f}):\n")
        for c in conv:
            f.write(f"  {c['tag']:8s} dx={c['dx']:.3f} nx={c['nx']:4d} npts={c['npts']:4d}  "
                    f"Cl={c['cl_mean']:+.4f}±{c['cl_std']:.4f}  Cd={c['cd_mean']:+.4f}±{c['cd_std']:.4f}  "
                    f"({c['elapsed']:.0f}s)\n")

    print(f"{name}: wrote figures/summary to {outdir}")

print("done")
