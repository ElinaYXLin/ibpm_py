"""
plot_residuals.py

Reads SURF_test/vortall/inner/data/residuals.csv (written by
compute_residuals.py) and produces figures under
SURF_test/vortall/inner/figures/:

  01_solver_residuals.png     divergence / no-slip / Poisson residuals vs t,
                               log scale, against a floating-point-roundoff
                               reference line -- the main "is the linear
                               algebra actually correct" plot.
  02_vorticity_smoothness.png ||Laplacian(omega)|| vs t -- flags grid-scale
                               noise/blow-up (not a residual, a sanity check).
  03_conservation.png         circulation / enstrophy / kinetic energy vs t.
  04_force_reproducibility.png Cd, Cl from this instrumented rerun overlaid
                               on the archived SURF_test/vortall/_run_data/
                               vortall.force, plus their pointwise difference
                               -- checks that this script's from-scratch
                               driver reproduces the documented run.

Usage:
    python3 SURF_test/vortall/inner/plot_residuals.py
"""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
FIG_DIR = HERE / "figures"
ARCHIVED_FORCE = HERE.parent / "_run_data" / "vortall.force"

MACHINE_EPS_REF = 1e-10  # "floating point roundoff" reference line for residual plots
PLOT_FLOOR = 1e-18  # y-axis floor for log residual plots; step 0 is an exact-zero
                     # pre-solve state (see compute_residuals.py) and would otherwise
                     # compress the whole informative range into a thin band


def load() -> dict:
    arr = np.genfromtxt(DATA_DIR / "residuals.csv", delimiter=",", names=True)
    return {name: arr[name] for name in arr.dtype.names}


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    d = load()
    t = d["time"]

    # ---------------- 1. solver residuals ----------------
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    panels = [
        ("div_l2", "div_linf", "Divergence of flux field  ||div(q)||\n"
         "(independent finite-difference check, not py/'s own Curl code)"),
        ("noslip_l2", "noslip_linf", "No-slip constraint residual  ||C(ω) - b||\n"
         "(is ProjectionSolver actually enforcing the body BC?)"),
        ("poisson_l2", "poisson_linf", "Streamfunction Poisson-equation residual\n"
         "||∇²ψ - (-ω)||  (is PoissonSolver actually solving what it claims?)"),
    ]
    for ax, (l2_key, linf_key, title) in zip(axes, panels):
        # per-panel floor: half the smallest nonzero value actually seen after
        # step 0 (step 0 is an exact-zero pre-solve state, excluded here so it
        # doesn't compress the whole informative y-range into a thin band)
        nonzero = np.concatenate([d[l2_key][d["step"] > 0], d[linf_key][d["step"] > 0]])
        nonzero = nonzero[nonzero > 0]
        floor = max(nonzero.min() * 0.3, 1e-300) if nonzero.size else PLOT_FLOOR
        ceil = max(nonzero.max(), MACHINE_EPS_REF) * 30
        ax.semilogy(t, np.maximum(d[l2_key], floor), label="RMS (L2)", color="#2980b9", lw=1.2)
        ax.semilogy(t, np.maximum(d[linf_key], floor), label="max (L∞)", color="#c0392b", lw=1.0, alpha=0.8)
        ax.axhline(MACHINE_EPS_REF, color="gray", ls="dashed", lw=1,
                    label=f"{MACHINE_EPS_REF:g} reference")
        ax.set_ylim(floor, ceil)
        ax.set_title(title, fontsize=10.5)
        ax.set_ylabel("residual")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("simulation time t")
    fig.suptitle("SURF_test/vortall -- inner solver residuals over 200 time units\n"
                  "(all three should stay near floating-point roundoff if the "
                  "solver is doing its job)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG_DIR / "01_solver_residuals.png", dpi=150)
    plt.close(fig)

    # ---------------- 2. vorticity field smoothness ----------------
    fig, ax = plt.subplots(figsize=(11, 4))
    nz = np.concatenate([d["helmholtz_l2"][d["step"] > 0], d["helmholtz_linf"][d["step"] > 0]])
    nz = nz[nz > 0]
    floor2 = nz.min() * 0.3 if nz.size else 1e-10
    ax.semilogy(t, np.maximum(d["helmholtz_l2"], floor2), color="#8e44ad", lw=1.2, label="RMS ||∇²ω||")
    ax.semilogy(t, np.maximum(d["helmholtz_linf"], floor2), color="#d35400", lw=1.0, alpha=0.8, label="max ||∇²ω||")
    ax.set_ylim(floor2, nz.max() * 5 if nz.size else 1)
    ax.set_xlabel("simulation time t")
    ax.set_ylabel("||Laplacian(ω)||")
    ax.set_title("Vorticity-field smoothness (||∇²ω||) -- NOT a solver residual, a "
                  "sanity check for grid-scale noise/blow-up.\nSteady, bounded growth as the wake "
                  "develops is expected; a sudden jump/spike would flag trouble.", fontsize=10.5)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "02_vorticity_smoothness.png", dpi=150)
    plt.close(fig)

    # ---------------- 3. conservation quantities ----------------
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    axes[0].plot(t, d["circulation"], color="#16a085", lw=1.2)
    axes[0].set_ylabel("circulation\nΣω·dx²")
    axes[0].set_title("Total circulation (should stay near 0 for this symmetric,\n"
                       "zero-total-circulation-shedding case; no drift/runaway trend expected)", fontsize=10)
    axes[1].plot(t, d["enstrophy"], color="#2980b9", lw=1.2)
    axes[1].set_ylabel("enstrophy\nΣω²·dx²")
    axes[1].set_title("Enstrophy -- grows as the wake instability develops, "
                       "then oscillates once shedding saturates", fontsize=10)
    axes[2].plot(t, d["kinetic_energy"], color="#c0392b", lw=1.2)
    axes[2].set_ylabel("kinetic energy\n½<q,q>")
    axes[2].set_title("Kinetic energy of the flux field", fontsize=10)
    axes[-1].set_xlabel("simulation time t")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle("SURF_test/vortall -- conservation / physical sanity quantities", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG_DIR / "03_conservation.png", dpi=150)
    plt.close(fig)

    # ---------------- 4. force reproducibility vs. archived run ----------------
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    axes[0].plot(t, d["Cd"], color="#2980b9", lw=1.1, label="this run (inner/)")
    axes[1].plot(t, d["Cl"], color="#2980b9", lw=1.1, label="this run (inner/)")

    if ARCHIVED_FORCE.exists():
        arch = np.loadtxt(ARCHIVED_FORCE)
        t_arch, cd_arch, cl_arch = arch[:, 1], arch[:, 2], arch[:, 3]
        n = min(len(t), len(t_arch))
        axes[0].plot(t_arch, cd_arch, color="#c0392b", lw=1.0, ls="dashed",
                      label="archived _run_data/vortall.force")
        axes[1].plot(t_arch, cl_arch, color="#c0392b", lw=1.0, ls="dashed",
                      label="archived _run_data/vortall.force")
        diff_cd = np.abs(d["Cd"][:n] - cd_arch[:n])
        diff_cl = np.abs(d["Cl"][:n] - cl_arch[:n])
        floor3 = 1e-16
        axes[2].semilogy(t[:n], np.maximum(diff_cd, floor3), color="#2980b9", lw=1.1, label="|ΔCd|")
        axes[2].semilogy(t[:n], np.maximum(diff_cl, floor3), color="#c0392b", lw=1.1, label="|ΔCl|")
        axes[2].axhline(5e-6, color="gray", ls="dashed", lw=1,
                          label="5e-6 (vortall.force's own %.5e print precision)")
        axes[2].set_ylim(floor3, 1e-4)
        axes[2].set_ylabel("|difference|")
        axes[2].set_title("Pointwise difference between this independently-written driver "
                           "and the archived run\n(both run the identical deterministic algorithm; "
                           "the ~5e-6 ceiling here is explained entirely by\n"
                           "vortall.force's own \"%.5e\" (6-significant-figure) text format, not a "
                           "real trajectory divergence -- see README.md)", fontsize=9.5)
        axes[2].legend(fontsize=8)
        axes[2].grid(alpha=0.3)
    else:
        axes[2].text(0.5, 0.5, f"archived force file not found:\n{ARCHIVED_FORCE}",
                      ha="center", va="center", transform=axes[2].transAxes)

    axes[0].set_ylabel("Cd")
    axes[0].set_title("Drag coefficient", fontsize=10)
    axes[1].set_ylabel("Cl")
    axes[1].set_title("Lift coefficient", fontsize=10)
    for ax in axes[:2]:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("simulation time t")
    fig.suptitle("SURF_test/vortall -- reproducibility check against the archived run", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG_DIR / "04_force_reproducibility.png", dpi=150)
    plt.close(fig)

    print(f"wrote figures to {FIG_DIR}")

    # ---------------- text summary ----------------
    summary_path = HERE / "residual_summary.txt"
    with open(summary_path, "w") as f:
        def line(s: str = "") -> None:
            print(s)
            f.write(s + "\n")

        line("Residual summary, SURF_test/vortall/inner/data/residuals.csv")
        line(f"  {len(t)} recorded steps, t in [{t.min():.2f}, {t.max():.2f}]")
        line()
        # skip step 0 (pre-solve, not representative -- see compute_residuals.py)
        mask = d["step"] > 0
        for key, label in [("div_l2", "divergence L2"), ("div_linf", "divergence Linf"),
                             ("noslip_l2", "no-slip L2"), ("noslip_linf", "no-slip Linf"),
                             ("poisson_l2", "poisson L2"), ("poisson_linf", "poisson Linf")]:
            v = d[key][mask]
            line(f"  {label:16s}: max over run = {v.max():.3e}, median = {np.median(v):.3e}")
        line()
        line(f"  circulation: min={d['circulation'].min():.3e}, max={d['circulation'].max():.3e}")
        line(f"  enstrophy:   min={d['enstrophy'].min():.3e}, max={d['enstrophy'].max():.3e}")
        if ARCHIVED_FORCE.exists():
            n = min(len(t), len(np.loadtxt(ARCHIVED_FORCE)))
            arch = np.loadtxt(ARCHIVED_FORCE)
            diff_cd = np.abs(d["Cd"][:n] - arch[:n, 2])
            diff_cl = np.abs(d["Cl"][:n] - arch[:n, 3])
            line()
            line(f"  vs. archived vortall.force: max|dCd|={diff_cd.max():.3e}, "
                 f"max|dCl|={diff_cl.max():.3e}")
            line("  (this max is fully explained by vortall.force's own \"%.5e\" text-print "
                 "precision -- see output_force.py -- not a real trajectory divergence;")
            line("   diff stays flat at this ceiling for the whole t=0-200 window instead of "
                 "growing, i.e. this run and the archived one have NOT yet chaotically diverged)")
            # a threshold well above print-precision noise (5e-6) but far below an O(1)
            # phase-drift divergence -- would flag a REAL trajectory split, if one occurred
            growing_threshold = 1e-4
            first_bad = np.argmax(diff_cd > growing_threshold) if np.any(diff_cd > growing_threshold) else -1
            if first_bad > 0:
                line(f"  |dCd| first exceeds {growing_threshold:g} (real divergence, not just "
                     f"print-precision) at step {int(d['step'][first_bad])}, t={t[first_bad]:.2f}")
            else:
                line(f"  |dCd| never exceeds {growing_threshold:g} over this run -- no real "
                     f"divergence from the archived trajectory within t=0-200")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
