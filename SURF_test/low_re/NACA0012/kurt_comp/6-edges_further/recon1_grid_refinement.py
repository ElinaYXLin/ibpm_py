"""
recon1_grid_refinement.py

Reconciliation thread #1 from 5-leading_edge/README.md's "Open thread": at
Re=1000 (5-leading_edge's Group 2, y=0 lineout metric), the LE peak SHRINKS
with grid refinement (22.2->8.5->2.95). At Re=500 (the older
../2-leading_edge_investigation, 2-D field-max metric), the LE peak GROWS
(64.7->68.3->74.7). Two things differ at once: Reynolds number AND how
"peak" is measured. This computes BOTH metrics at BOTH Re from data that
ALREADY EXISTS on disk (zero new runs) -- if that alone explains the
disagreement, no new simulation is needed at all.

Uses py_static's State loader directly on ../2-leading_edge_investigation's
existing snapshots (produced by the OLDER `py`/`src` ports, never
modified here) -- confirmed binary-compatible (py_static reads them and
reproduces the exact reported peak, 64.711 vs the old investigation's
documented 64.7).

Usage: python3 recon1_grid_refinement.py
Output: figures/recon1_grid_refinement.png, data/recon1_grid_refinement.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)  # 2-D window for the field-max metric


def field_max_metric(X, Y, om):
    """2-D field maximum near the LE, matching the OLDER investigation's
    method (a raw max|omega| over a window around the LE, not a 1-D cut)."""
    _, _, sub, _, _ = c.window(X, Y, om, LE_XLIM, LE_YLIM)
    return float(np.abs(sub).max())


def lineout_metric(X, Y, om):
    """y=0 lineout maximum, matching 5-leading_edge's Group 2 method."""
    ys = Y[0, :]
    iy0 = c.nearest_index(ys, 0.0)
    xs = X[:, 0]
    m = (xs >= LE_XLIM[0]) & (xs <= LE_XLIM[1])
    return float(np.abs(om[:, iy0][m]).max())


def main():
    rows = []

    # ---------------- Re=500 (existing ../2-leading_edge_investigation data) ----------------
    re500_cases = [
        (0.02, c.PRIOR / "_run_data" / "le_uniform_baseline_snap", 3000),
        (0.01, c.PRIOR / "_run_data" / "gridconv_dx0.01_snap", 6000),
        (0.005, c.PRIOR / "_run_data" / "gridconv_dx0.005_snap", 12000),
    ]
    for dx, run_dir, step in re500_cases:
        X, Y = c.grid_xy(dx)
        om = c.load_omega(run_dir, step)
        rows.append(dict(re=500, dx=dx, field_max=field_max_metric(X, Y, om),
                          lineout_max=lineout_metric(X, Y, om)))

    # ---------------- Re=1000 (existing 5-leading_edge / 1-paper_based data) ----------------
    re1000_cases = [
        (0.02, c.KURT / "1-paper_based" / "runs" / "dx0.020" / "steady_py_a00", 3000),
        (0.01, c.KURT5 / "runs" / "grid_refine" / "dx0.0100", 6000),
        (0.005, c.KURT5 / "runs" / "grid_refine" / "dx0.0050", 12000),
    ]
    for dx, run_dir, step in re1000_cases:
        X, Y = c.grid_xy(dx)
        om = c.load_omega(run_dir, step)
        rows.append(dict(re=1000, dx=dx, field_max=field_max_metric(X, Y, om),
                          lineout_max=lineout_metric(X, Y, om)))

    with open(c.DATA / "recon1_grid_refinement.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["re", "dx", "field_max", "lineout_max"])
        w.writeheader(); w.writerows(rows)

    print("Reconciliation 1: LE peak vs dx, cross-matched Re x metric (zero new runs)")
    print(f"{'Re':>6} {'dx':>7} {'field_max (2-D)':>16} {'lineout_max (1-D)':>18}")
    for r in rows:
        print(f"{r['re']:>6} {r['dx']:>7} {r['field_max']:>16.3f} {r['lineout_max']:>18.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = {500: "#c0392b", 1000: "#2980b9"}
    for ax, metric, title in [(axes[0], "field_max", "2-D field-max metric\n(the OLDER investigation's method)"),
                               (axes[1], "lineout_max", "y=0 lineout-max metric\n(5-leading_edge's method)")]:
        for re in (500, 1000):
            sub = [r for r in rows if r["re"] == re]
            sub.sort(key=lambda r: -r["dx"])
            dxs = [r["dx"] for r in sub]
            vals = [r[metric] for r in sub]
            ax.plot(dxs, vals, "o-", color=colors[re], label=f"Re={re}", lw=1.8, ms=7)
        ax.set_xscale("log"); ax.invert_xaxis()
        ax.set_xlabel("dx (log scale, refining -->)"); ax.set_ylabel("LE peak |omega|")
        ax.set_title(title, fontsize=10.5)
        ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.suptitle("Reconciliation 1: does the grid-refinement DIRECTION disagreement (Group 2 vs.\n"
                 "the Re=500 investigation) come from Reynolds number, or from how 'peak' is measured?",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "recon1_grid_refinement.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nwrote {out.name}")


if __name__ == "__main__":
    main()
