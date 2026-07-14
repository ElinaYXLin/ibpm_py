"""
gen_chaos_sensitivity_figs.py

Reads this directory's _run_data/ (produced by run_chaos_sensitivity.py)
plus the existing SURF_test/airfoils/LSAT-{SD7003,SD8000}/4-Re_sweep/Re200
restart data (no rerun needed for that part -- it already exists) and
writes every figure/number in this directory's README.md. Nothing here is
hand-drawn; every figure traces back to an actual py/ibpm.py or
build/ibpm run.

Outputs:
  re200_comparison.png              -- Cl(t), Cd(t), phase-space, Re=200, both airfoils
  phase_space_coarse_ext4000.png    -- phase-space, extended coarse-grid runs, both airfoils
  blowup_histogram.png              -- blow-up step histogram, 16-run perturbation ensemble
  SD8000_ext4000_py_snapshots.png, SD8000_ext4000_cpp_snapshots.png,
  SD7003_ext4000_py_snapshots.png   -- vorticity snapshots leading up to each blow-up

Usage: python3 SURF_test/airfoils/7-chaos_sensitivity/gen_chaos_sensitivity_figs.py
"""
from __future__ import annotations

import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
OUTDIR = REPO / "SURF_test" / "airfoils" / "7-chaos_sensitivity"
RUN_DATA = OUTDIR / "_run_data"
AIRFOILS_DIR = REPO / "SURF_test" / "airfoils"

sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

BLOWUP_THRESHOLD = 20.0  # |Cl| or |Cd| exceeding this = "blown up" (normal osc. amplitude is ~1-5)

# Grid for the coarse (dx=0.04) cases: nx=150, ny=75, length=6, xoffset=-2, yoffset=-1.5
NX, NY = 150, 75
LENGTH, XOFFSET, YOFFSET = 6.0, -2.0, -1.5
DX = LENGTH / NX
XS = XOFFSET + np.arange(1, NX) * DX
YS = YOFFSET + np.arange(1, NY) * DX
XCOARSE, YCOARSE = np.meshgrid(XS, YS, indexing="ij")


def load_force(p):
    d = np.loadtxt(p)
    if d.ndim == 1:
        d = d[None, :]
    return d[:, 1], d[:, 2], d[:, 3]  # t, Cd, Cl


def find_blowup_step(force_path, threshold=BLOWUP_THRESHOLD):
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    step, cd, cl = d[:, 0], d[:, 2], d[:, 3]
    mask = (np.abs(cl) > threshold) | (np.abs(cd) > threshold)
    if not mask.any():
        return None
    return int(step[np.argmax(mask)])


# ======================================================================
# 1. Re=200 comparison (existing data, no rerun) + phase space
# ======================================================================
def gen_re200_comparison():
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for row, name in enumerate(["SD7003", "SD8000"]):
        py_dir = AIRFOILS_DIR / f"LSAT-{name}" / "4-Re_sweep" / "_run_data" / "Re200"
        cpp_dir = AIRFOILS_DIR / f"LSAT-{name}" / "4-Re_sweep" / "_run_data_cpp" / "Re200"
        t_py, cd_py, cl_py = load_force(py_dir / "run.force")
        t_cpp, cd_cpp, cl_cpp = load_force(cpp_dir / "run.force")

        ax = axes[row, 0]
        ax.plot(t_py, cl_py, color="C0", lw=0.8, label="py Cl")
        ax.plot(t_cpp, cl_cpp, color="C3", lw=0.8, ls="--", label="cpp Cl")
        ax.set_title(f"{name} Re=200: Cl(t)"); ax.set_xlabel("t"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axes[row, 1]
        ax.plot(t_py, cd_py, color="C0", lw=0.8, label="py Cd")
        ax.plot(t_cpp, cd_cpp, color="C3", lw=0.8, ls="--", label="cpp Cd")
        ax.set_title(f"{name} Re=200: Cd(t)"); ax.set_xlabel("t"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axes[row, 2]
        ax.plot(cd_py, cl_py, color="C0", lw=0.8, label="py")
        ax.plot(cd_cpp, cl_cpp, color="C3", lw=0.8, ls="--", label="cpp")
        ax.set_xlabel("Cd"); ax.set_ylabel("Cl")
        ax.set_title(f"{name} Re=200: phase space Cl vs Cd"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

        n = min(len(cl_py), len(cl_cpp))
        maxdiff_cl = float(np.abs(cl_py[:n] - cl_cpp[:n]).max())
        maxdiff_cd = float(np.abs(cd_py[:n] - cd_cpp[:n]).max())
        print(f"{name} Re=200: max|Cl diff|={maxdiff_cl:.3e}  max|Cd diff|={maxdiff_cd:.3e}  (n={n} steps)")

    fig.suptitle("Re=200 py vs cpp: force traces + phase space "
                 "(no experimental reference needed -- both-implementation agreement check)")
    fig.tight_layout()
    outp = OUTDIR / "re200_comparison.png"
    fig.savefig(outp, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outp}")


# ======================================================================
# 2. Extended-to-4000-steps: blow-up detection + phase space
# ======================================================================
def gen_extended_runs():
    blowup_steps = {}
    for name in ["SD7003", "SD8000"]:
        for impl in ["py", "cpp"]:
            fp = RUN_DATA / f"{name}_ext4000_{impl}" / "run.force"
            step = find_blowup_step(fp)
            blowup_steps[(name, impl)] = step
            if step is not None:
                print(f"{name} {impl}: blows up at step {step}")
            else:
                print(f"{name} {impl}: no blow-up through 4000 steps")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, name in zip(axes, ["SD7003", "SD8000"]):
        for impl, color in [("py", "C0"), ("cpp", "C3")]:
            d = np.loadtxt(RUN_DATA / f"{name}_ext4000_{impl}" / "run.force")
            cd, cl = d[:, 2], d[:, 3]
            mask = (np.abs(cd) < 10) & (np.abs(cl) < 10)
            ax.plot(cd[mask], cl[mask], color=color, lw=0.5, alpha=0.8, label=impl)
        ax.set_xlabel("Cd"); ax.set_ylabel("Cl")
        ax.set_title(f"{name} coarse grid (dx=0.04), extended to 4000 steps\n"
                     "phase space Cl vs Cd (|Cd|,|Cl|<10 shown)")
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.suptitle("Phase-space trajectories: broadband/chaotic, no clean limit cycle "
                 "(contrast with Re=200 spiral-to-point above)")
    fig.tight_layout()
    outp = OUTDIR / "phase_space_coarse_ext4000.png"
    fig.savefig(outp, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outp}")
    return blowup_steps


# ======================================================================
# 3. Vorticity snapshots leading up to each blow-up
# ======================================================================
def gen_blowup_snapshots(blowup_steps):
    cases = [
        ("SD8000", "py", [2800, 2825, 2850, 2875, 2900]),
        ("SD8000", "cpp", [2875, 2900, 2925, 2950, 2975]),
        ("SD7003", "py", [3850, 3875, 3900, 3925, 3975]),
    ]
    for name, impl, steps in cases:
        blowup_step = blowup_steps.get((name, impl))
        if blowup_step is None:
            print(f"{name} {impl}: no recorded blow-up, skipping snapshot figure")
            continue
        rundir = RUN_DATA / f"{name}_ext4000_{impl}"
        fig, axes = plt.subplots(1, len(steps), figsize=(3.2 * len(steps), 4))
        for ax, s in zip(axes, steps):
            fp = rundir / f"run{s:05d}.bin"
            if not fp.exists():
                ax.set_visible(False)
                continue
            st = State(filename=str(fp))
            w = st.omega._data[0]
            ax.contourf(XCOARSE, YCOARSE, np.clip(w, -8.0, 8.0), levels=41, cmap="RdBu_r", extend="both")
            tag = "STABLE" if s < blowup_step else "POST-BLOWUP"
            ax.set_title(f"t={s*0.01:.2f} ({tag})", fontsize=9)
            ax.set_aspect("equal")
        fig.suptitle(f"{name} {impl} (blow-up step {blowup_step}): vorticity snapshots "
                     "leading up to blow-up (every 25 steps)")
        fig.tight_layout()
        outp = OUTDIR / f"{name}_ext4000_{impl}_snapshots.png"
        fig.savefig(outp, dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {outp}")


# ======================================================================
# 4. Perturbation ensemble: blow-up step histogram
# ======================================================================
def gen_blowup_histogram():
    rel_perturbations = [-1e-5, -5e-6, -2e-6, -1e-6, -5e-7, -2e-7, -1e-7, -5e-8,
                          5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5]
    blowup_steps = []
    n_total = len(rel_perturbations)
    for i in range(n_total):
        fp = RUN_DATA / f"SD8000_perturb_{i:02d}" / "run.force"
        step = find_blowup_step(fp)
        if step is not None:
            blowup_steps.append(step)
            print(f"  perturb {i:02d} (rel={rel_perturbations[i]:+.1e}): blows up at step {step}")
        else:
            print(f"  perturb {i:02d} (rel={rel_perturbations[i]:+.1e}): no blow-up through 4000 steps")

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(blowup_steps, bins=12, color="C3", edgecolor="black", alpha=0.85)
    ax.axvline(3000, color="k", ls="--", lw=1.2,
               label="original run length (3000 steps)\n(where the first-documented blow-up was seen)")
    ax.set_xlabel("timestep at which |Cl| or |Cd| first exceeds 20")
    ax.set_ylabel("count (of 16 runs)")
    ax.set_title(f"SD8000 coarse grid (dx=0.04), C++ only: blow-up step vs. tiny Re perturbation\n"
                 f"Re perturbed by only 5e-8 to 1e-5 (relative) around 60800 -- {len(blowup_steps)}/{n_total} "
                 f"runs blow up (spread: {max(blowup_steps)-min(blowup_steps)} steps), "
                 f"{n_total-len(blowup_steps)}/{n_total} stay bounded through 4000 steps")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    outp = OUTDIR / "blowup_histogram.png"
    fig.savefig(outp, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outp}")


def main():
    gen_re200_comparison()
    blowup_steps = gen_extended_runs()
    gen_blowup_snapshots(blowup_steps)
    gen_blowup_histogram()
    print("done")


if __name__ == "__main__":
    main()
