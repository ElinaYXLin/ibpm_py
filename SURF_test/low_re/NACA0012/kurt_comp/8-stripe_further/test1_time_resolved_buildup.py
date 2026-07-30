"""
test1_time_resolved_buildup.py

Follow-up test #1 from ../7-stripe_investigation/README.md's Proposals
section: Test 8's "seeded immediately by projection, but full multi-chord
reach builds up over many iterations via the elliptic solve" claim
(H6c) was inferred from exactly TWO points -- step 1 and the converged
step 3000. The baseline dx=0.02 run already saved a restart snapshot
every 250 steps; this plots upstream enstrophy and reach L_up vs. step
number through the whole run, turning that inferred claim into an
observed curve. Zero new runs.

The existing 250-step cadence turned out to be too coarse: enstrophy and
L_up are ALREADY at ~100% of the converged value by the very first
snapshot (step 250), so the coarse data alone only proves buildup
happens somewhere in [step 1, step 250] -- it doesn't resolve HOW. This
script therefore also launches one cheap new run (same case exactly:
dx=0.02, alpha=0, Re=1000, dt=0.01, zero IC -- deterministic, so it's
the identical trajectory, just saved every 5 steps instead of every 250
for the first 300 steps) to actually resolve the buildup curve. Cheap:
dx=0.02, 300x150 grid, 300 steps.

Usage:
  python3 test1_time_resolved_buildup.py run       # launch the fine-cadence rerun
  python3 test1_time_resolved_buildup.py analyze   # analyze + plot (default)
Output: figures/test1_time_resolved_buildup.png,
        data/test1_time_resolved_buildup.csv
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

DX = 0.02
BUFFER_DX = 2.0
YLIM = (-0.5, 0.5)
NOISE_FLOOR = 1e-3
STEPS = list(range(0, 3001, 250))
RUN_DIR = c.KURT1 / "runs" / "dx0.020" / "steady_py_a00"

FINE_NSTEPS = 300
FINE_RESTART = 5
FINE_STEPS = list(range(0, FINE_NSTEPS + 1, FINE_RESTART))
FINE_RUN_DIR = c.RUNS / "fine_cadence_a00"


def do_run():
    if c.is_done(FINE_RUN_DIR, FINE_NSTEPS):
        print("fine-cadence rerun: already done")
        return
    print("fine-cadence rerun: launching (300 steps, saved every 5)...", flush=True)
    ok, elapsed = c.run_case("py", c.BASE_GEOM_DX002, FINE_RUN_DIR, 300, 150, 0.01, FINE_NSTEPS,
                               alpha=0.0, restart=FINE_RESTART)
    print(f"  {'OK' if ok else 'FAILED'} in {elapsed:.1f}s")


def main():
    X, Y = c.grid_xy(DX)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    rows = []
    for step in STEPS:
        om = c.load_omega(RUN_DIR, step)
        m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
        prof = c.upstream_profile(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
        L_up, _, _ = c.reach_L_up(prof["xs"], prof["max_abs"], x_le, NOISE_FLOOR)
        rows.append(dict(step=step, enstrophy=m["enstrophy"], peak=m["peak"],
                          L_up_chord=L_up, L_up_cells=L_up / DX))
        print(f"step={step}: enstrophy={m['enstrophy']:.6f}, peak={m['peak']:.4f}, "
              f"L_up={L_up:.4f}c ({L_up/DX:.2f} cells)")

    # fine-cadence rerun (same trajectory, deterministic zero IC, just saved more often)
    fine_rows = []
    if c.is_done(FINE_RUN_DIR, FINE_NSTEPS):
        for step in FINE_STEPS:
            om = c.load_omega(FINE_RUN_DIR, step)
            m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            prof = c.upstream_profile(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            L_up, _, _ = c.reach_L_up(prof["xs"], prof["max_abs"], x_le, NOISE_FLOOR)
            fine_rows.append(dict(step=step, enstrophy=m["enstrophy"], peak=m["peak"],
                                   L_up_chord=L_up, L_up_cells=L_up / DX))
        print(f"\nfine-cadence rerun: {len(fine_rows)} snapshots loaded (0-{FINE_NSTEPS}, every {FINE_RESTART})")
    else:
        print("\nfine-cadence rerun not found -- run `python3 test1_time_resolved_buildup.py run` first")

    with open(c.DATA / "test1_time_resolved_buildup.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["step", "enstrophy", "peak", "L_up_chord", "L_up_cells"])
        w.writeheader(); w.writerows(rows)
    if fine_rows:
        with open(c.DATA / "test1_time_resolved_buildup_fine.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["step", "enstrophy", "peak", "L_up_chord", "L_up_cells"])
            w.writeheader(); w.writerows(fine_rows)
    print("wrote test1_time_resolved_buildup.csv" + (" and _fine.csv" if fine_rows else ""))

    steps = [r["step"] for r in rows]
    ens = [r["enstrophy"] for r in rows]
    lup = [r["L_up_chord"] for r in rows]
    converged_ens = ens[-1]
    converged_lup = lup[-1]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for col, (key, label, color) in enumerate([("enstrophy", "upstream enstrophy", "#c0392b"),
                                                  ("L_up_chord", "reach L_up (chord)", "#2980b9")]):
        converged = ens[-1] if key == "enstrophy" else lup[-1]
        # top row: full 0-3000 range
        ax = axes[0, col]
        ax.plot(steps, [r[key] for r in rows], "o-", color=color, lw=1.4, ms=5, alpha=0.5,
                label="coarse (every 250 steps, existing)")
        if fine_rows:
            ax.plot([r["step"] for r in fine_rows], [r[key] for r in fine_rows], "o-",
                     color=color, lw=1.8, ms=4, label="fine (every 5 steps, new)")
        ax.axhline(converged, color="gray", ls="--", lw=1, alpha=0.7, label="converged (step 3000)")
        ax.set_xlabel("timestep"); ax.set_ylabel(label)
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
        ax.set_title(f"{label} vs. timestep, full range", fontsize=10)
        # bottom row: zoomed to the fine-cadence window (0-300) where all the action is
        ax = axes[1, col]
        if fine_rows:
            fsteps = [r["step"] for r in fine_rows]
            ax.plot(fsteps, [r[key] for r in fine_rows], "o-", color=color, lw=1.8, ms=5)
            ax.axhline(converged, color="gray", ls="--", lw=1, alpha=0.7, label="converged (step 3000)")
            ax.axhline(0.9 * converged, color="gray", ls=":", lw=1, alpha=0.6, label="90% of converged")
            ax.set_xlim(0, FINE_NSTEPS)
        ax.set_xlabel("timestep (t = step x 0.01)"); ax.set_ylabel(label)
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
        ax.set_title(f"{label}, zoomed to 0-{FINE_NSTEPS} (fine cadence)", fontsize=10)

    fig.suptitle("Test 1: is the far-upstream reach built up gradually over many\n"
                 "iterations (H6c), or already fully developed early on?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = c.FIGS / "test1_time_resolved_buildup.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    # quantify: step at which enstrophy/L_up first reach 90% of converged value
    def step_to_90pct(step_list, vals, converged):
        thresh = 0.9 * converged
        for s, v in zip(step_list, vals):
            if v >= thresh:
                return s
        return None
    print(f"\n[coarse] enstrophy reaches 90% of converged by step "
          f"{step_to_90pct(steps, ens, converged_ens)}")
    print(f"[coarse] L_up reaches 90% of converged by step "
          f"{step_to_90pct(steps, lup, converged_lup)}")
    if fine_rows:
        fsteps = [r["step"] for r in fine_rows]
        fens = [r["enstrophy"] for r in fine_rows]
        flup = [r["L_up_chord"] for r in fine_rows]
        print(f"[fine] enstrophy reaches 90% of converged by step "
              f"{step_to_90pct(fsteps, fens, converged_ens)}")
        print(f"[fine] L_up reaches 90% of converged by step "
              f"{step_to_90pct(fsteps, flup, converged_lup)}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        do_run()
    else:
        main()
