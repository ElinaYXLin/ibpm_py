"""
gen_figs_6.py

Figures for recon2 and Groups B-E (recon1/A/F have their own gen code in
recon1_grid_refinement.py / testA_metric_and_geom_audit.py /
testF_conditioning.py). Run analyze_6.py first.

Usage: python3 gen_figs_6.py
Output: figures/*.png
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

C_FIELD = "#2980b9"
C_LINE = "#c0392b"


def load(name):
    p = c.DATA / f"{name}.csv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f))


def fig_recon2():
    rows = load("recon2_boundary_density")
    if not rows:
        return
    # augment with the 4 "existing data" points computed inline (see
    # analyze_6.analyze_recon2's printed context) -- hardcode here since
    # they come from separate, already-existing runs across two folders
    table = {
        ("LE-only", 500): 89.88, ("LE-only", 1000): 108.37,
        ("LE+TE", 500): 90.41, ("LE+TE", 1000): 109.87,
    }
    baseline = {500: 64.71, 1000: 71.66}
    fig, ax = plt.subplots(figsize=(8, 5.5))
    x = np.array([500, 1000])
    ax.plot(x, [baseline[re] for re in x], "k--", marker="s", label="baseline (no densification)", lw=1.5)
    ax.plot(x, [table[("LE-only", re)] for re in x], "-o", color="#2980b9", label="LE-only dense", lw=2, ms=8)
    ax.plot(x, [table[("LE+TE", re)] for re in x], "-o", color="#c0392b", label="LE+TE dense", lw=2, ms=8)
    ax.set_xticks(x); ax.set_xlabel("Reynolds number"); ax.set_ylabel("LE peak |omega| (2-D field-max metric)")
    ax.set_title("Reconciliation 2: densifying boundary points, crossed with Re\n"
                 "(both LE-only and LE+TE make it WORSE, at both Re -- not the 'big win'\n"
                 "previously reported using the lineout metric)", fontsize=10.5)
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    out = c.FIGS / "recon2_boundary_density.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_B1():
    rows = [r for r in load("testB1_phase_sweep") if r["impl"] == "py"]
    if not rows:
        return
    labels = [f"({r['shift_x']},{r['shift_y']})" for r in rows]
    field = [float(r["field_max"]) for r in rows]
    line = [float(r["lineout_max"]) for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, vals, title, color in [(axes[0], field, "2-D field-max metric", C_FIELD),
                                     (axes[1], line, "y=0 lineout-max metric", C_LINE)]:
        ax.bar(labels, vals, color=color)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("grid shift (fraction of dx in x, y)")
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle("Test B1: LE peak vs. grid sub-cell phase, NACA0012 held fixed\n"
                 "(same shape, same background grid -- only WHERE the grid sits relative to the body changes)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "testB1_phase_sweep.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_B2():
    rows = [r for r in load("testB2_phase_equalized") if r["impl"] == "py"]
    orig = {"naca0006": (106.74, 8.30), "naca0012": (71.66, 22.18), "naca0018": (64.74, 6.83)}
    if not rows:
        return
    names = [r["shape"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    x = np.arange(len(names))
    for ax, key, idx, title in [(axes[0], "field_max", 0, "2-D field-max metric"),
                                  (axes[1], "lineout_max", 1, "y=0 lineout-max metric")]:
        eq_vals = [float(r[key]) for r in rows]
        orig_vals = [orig[n][idx] for n in names]
        ax.bar(x - 0.18, orig_vals, width=0.36, color="#95a5a6", label="native phase")
        ax.bar(x + 0.18, eq_vals, width=0.36, color="#16a085", label="phase-equalized")
        ax.set_xticks(x); ax.set_xticklabels(names)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8)
    fig.suptitle("Test B2: does equalizing LE sub-cell phase change the\n"
                 "NACA0006/0012/0018 trend?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "testB2_phase_equalized.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_C1():
    rows = [r for r in load("testC1_shape_refinement") if r["impl"] == "py"]
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"naca0006": "#2980b9", "naca0018": "#c0392b"}
    # 0012's existing Group-2 data (field_max recomputed via recon1's method)
    naca0012_dx = [0.02, 0.01, 0.005]
    naca0012_field = [71.664, 86.561, 97.107]  # from recon1_grid_refinement.csv
    for ax, key, title in [(axes[0], "field_max", "2-D field-max metric"),
                             (axes[1], "lineout_max", "y=0 lineout-max metric")]:
        for shape in ("naca0006", "naca0018"):
            sub = [r for r in rows if r["shape"] == shape]
            sub.sort(key=lambda r: -float(r["dx"]))
            dxs = [float(r["dx"]) for r in sub]
            vals = [float(r[key]) for r in sub]
            ax.plot(dxs, vals, "o-", color=colors[shape], label=shape, lw=1.8, ms=7)
        if key == "field_max":
            ax.plot(naca0012_dx, naca0012_field, "o-", color="#16a085", label="naca0012 (existing)", lw=1.8, ms=7)
        ax.set_xscale("log"); ax.invert_xaxis()
        ax.set_xlabel("dx (log, refining -->)"); ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
    fig.suptitle("Test C1: per-shape grid refinement -- does the converged (fine-grid)\n"
                 "peak recover a monotonic sharper-is-worse trend?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "testC1_shape_refinement.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_C2():
    rows = [r for r in load("testC2_thickness_family") if r["impl"] == "py"]
    if not rows:
        return
    rows.sort(key=lambda r: float(r["r_le_over_dx"]))
    r_le = [float(r["r_le_over_dx"]) for r in rows]
    field = [float(r["field_max"]) for r in rows]
    line = [float(r["lineout_max"]) for r in rows]
    names = [r["shape"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, vals, title, color in [(axes[0], field, "2-D field-max metric", C_FIELD),
                                     (axes[1], line, "y=0 lineout-max metric", C_LINE)]:
        ax.plot(r_le, vals, "o-", color=color, lw=1.8, ms=8)
        for n, x, y in zip(names, r_le, vals):
            ax.annotate(n.replace("naca", ""), (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax.axvline(1.0, color="0.6", lw=1, ls="--")
        ax.set_xlabel("r_LE / dx"); ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    fig.suptitle("Test C2: continuous thickness family -- LE peak vs. nose radius / grid spacing", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    out = c.FIGS / "testC2_thickness_family.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_D1():
    rows = [r for r in load("testD1_point_density") if r["impl"] == "py"]
    if not rows:
        return
    rows.sort(key=lambda r: float(r["density_factor"]))
    diverged = [r["density_factor"] for r in rows if r["field_max"].lower() == "nan"]
    rows = [r for r in rows if r["field_max"].lower() != "nan"]
    factors = [float(r["density_factor"]) for r in rows]
    field = [float(r["field_max"]) for r in rows]
    line = [float(r["lineout_max"]) for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, vals, title, color in [(axes[0], field, "2-D field-max metric", C_FIELD),
                                     (axes[1], line, "y=0 lineout-max metric", C_LINE)]:
        ax.plot(factors, vals, "o-", color=color, lw=1.8, ms=8)
        ax.set_xscale("log")
        ax.set_xlabel("LE point-density factor (dx/ds, log scale)"); ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3, which="both")
    note = f"factor(s) {', '.join(diverged)}x: run diverged to NaN (numerical breakdown)" if diverged else ""
    fig.suptitle("Test D1: LE point density vs. peak (NACA0012, LE only, background dx=0.02 fixed)\n" + note,
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    out = c.FIGS / "testD1_point_density.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_E():
    rows1 = [r for r in load("testE1_decouple") if r["impl"] == "py"]
    rows2 = [r for r in load("testE2_common_TE") if r["impl"] == "py"]
    if not rows1 and not rows2:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    if rows1:
        ax = axes[0]
        names = [r["case"].replace("naca0012_dx0.0200_", "") for r in rows1]
        le = [float(r["le_field_max"]) for r in rows1]
        te = [float(r["te_field_max"]) for r in rows1]
        x = np.arange(len(names))
        ax.bar(x - 0.18, le, width=0.36, color="#2980b9", label="LE field_max")
        ax.bar(x + 0.18, te, width=0.36, color="#c0392b", label="TE field_max")
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
        ax.set_title("E1: decoupled front-only / TE-only variants", fontsize=10)
        ax.legend(fontsize=8)
    if rows2:
        ax = axes[1]
        names = [r["base_shape"] for r in rows2]
        le = [float(r["le_field_max"]) for r in rows2]
        te = [float(r["te_field_max"]) for r in rows2]
        x = np.arange(len(names))
        ax.bar(x - 0.18, le, width=0.36, color="#2980b9", label="LE field_max")
        ax.bar(x + 0.18, te, width=0.36, color="#c0392b", label="TE field_max")
        ax.set_xticks(x); ax.set_xticklabels([f"{n}+roundTE" for n in names])
        ax.set_title("E2: common-TE resweep (removes TE variation)", fontsize=10)
        ax.legend(fontsize=8)
    fig.suptitle("Group E: isolating the trailing-edge coupling confound", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    out = c.FIGS / "testE_decouple_and_commonTE.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    fig_recon2()
    fig_B1()
    fig_B2()
    fig_C1()
    fig_C2()
    fig_D1()
    fig_E()
