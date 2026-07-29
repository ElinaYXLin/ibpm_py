"""
recon1_le_te_field_max.py

Extension of recon1_grid_refinement.py: that script's left panel plots
the 2-D field-max metric (established, per Group G/Test A-F, as the
reliable metric -- unlike the y=0 lineout, which Groups B/C/recon1/recon2
all found unreliable) vs. dx, for the LEADING edge only, at both Re=500
and Re=1000. This reproduces the exact same experiment (same runs, same
dx sweep, same field-max metric, same Re=500/1000 overlay) but puts LE
and TE side by side, so both edges' grid-refinement trend under the
correct metric can be read directly, at both Reynolds numbers.

Zero new runs -- reuses exactly the same on-disk snapshots as
recon1_grid_refinement.py (../2-leading_edge_investigation for Re=500,
../5-leading_edge / ../1-paper_based for Re=1000); the TE window is just
a second crop of fields already being loaded, not a new simulation.

Usage: python3 recon1_le_te_field_max.py
Output: figures/recon1_le_te_field_max.png, data/recon1_le_te_field_max.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
TE_XLIM, TE_YLIM = (0.65, 1.35), (-0.25, 0.25)


def field_max_metric(X, Y, om, xlim, ylim):
    _, _, sub, _, _ = c.window(X, Y, om, xlim, ylim)
    return float(np.abs(sub).max())


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
        rows.append(dict(re=500, dx=dx,
                          le_field_max=field_max_metric(X, Y, om, LE_XLIM, LE_YLIM),
                          te_field_max=field_max_metric(X, Y, om, TE_XLIM, TE_YLIM)))

    # ---------------- Re=1000 (existing 5-leading_edge / 1-paper_based data) ----------------
    re1000_cases = [
        (0.02, c.KURT / "1-paper_based" / "runs" / "dx0.020" / "steady_py_a00", 3000),
        (0.01, c.KURT5 / "runs" / "grid_refine" / "dx0.0100", 6000),
        (0.005, c.KURT5 / "runs" / "grid_refine" / "dx0.0050", 12000),
    ]
    for dx, run_dir, step in re1000_cases:
        X, Y = c.grid_xy(dx)
        om = c.load_omega(run_dir, step)
        rows.append(dict(re=1000, dx=dx,
                          le_field_max=field_max_metric(X, Y, om, LE_XLIM, LE_YLIM),
                          te_field_max=field_max_metric(X, Y, om, TE_XLIM, TE_YLIM)))

    with open(c.DATA / "recon1_le_te_field_max.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["re", "dx", "le_field_max", "te_field_max"])
        w.writeheader(); w.writerows(rows)

    print("Recon1 LE/TE (2-D field-max metric only), both Re, vs dx (zero new runs)")
    print(f"{'Re':>6} {'dx':>7} {'LE field_max':>14} {'TE field_max':>14}")
    for r in rows:
        print(f"{r['re']:>6} {r['dx']:>7} {r['le_field_max']:>14.3f} {r['te_field_max']:>14.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = {500: "#c0392b", 1000: "#2980b9"}
    for ax, key, title in [(axes[0], "le_field_max", "Leading edge"),
                            (axes[1], "te_field_max", "Trailing edge")]:
        for re in (500, 1000):
            sub = [r for r in rows if r["re"] == re]
            sub.sort(key=lambda r: -r["dx"])
            dxs = [r["dx"] for r in sub]
            vals = [r[key] for r in sub]
            ax.plot(dxs, vals, "o-", color=colors[re], label=f"Re={re}", lw=1.8, ms=7)
        ax.set_xscale("log"); ax.invert_xaxis()
        ax.set_xlabel("dx (log scale, refining -->)"); ax.set_ylabel("peak |omega| (2-D field-max)")
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.suptitle("Recon1 (LE/TE side by side): 2-D field-max peak vs. grid refinement,\n"
                 "Re=500 and Re=1000 -- the metric established as reliable (Groups B/C/G, recon1/recon2)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "recon1_le_te_field_max.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nwrote {out.name}")


if __name__ == "__main__":
    main()
