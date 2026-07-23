"""
testA_metric_and_geom_audit.py

Group A: rule out metric fragility in Test 3b's non-monotonic thickness
trend (NACA0006/0012/0018 LE peaks: 8.3/22.2/6.8) before attributing it to
physics or geometry. Zero new runs -- reuses ../5-leading_edge's existing
Group 3b fields and geometry files.

A1 -- recompute the LE quantity 4 ways: 2-D window max (already used),
integrated enstrophy over the window, RMS of the window, and the spatial
extent (area) where |omega| exceeds a threshold. If the non-monotonicity
survives all four, it's real; if only the point-peak is non-monotonic
while integral metrics are monotonic, the "trend" is measurement noise
from point-sampling a grid-scale oscillation.

A2 -- geometry-only audit: each airfoil's exact LE sub-cell offset (LE x
position modulo dx) and residual boundary asymmetry about y=0. Feeds
Group B (phase) and Group F.

Usage: python3 testA_metric_and_geom_audit.py
Output: data/testA1_metric_robustness.csv, data/testA2_geometry_audit.csv,
        figures/testA1_metric_robustness.png
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

DX = 0.02
LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
V_THRESHOLD = 3.0  # |omega| threshold for the "extent" metric

SHAPES = {
    "naca0006": c.KURT5 / "geom" / "naca0006_dx0.0200.geom",
    "naca0012": c.BASE_GEOM_DX002,
    "naca0018": c.KURT5 / "geom" / "naca0018_dx0.0200.geom",
}
RUN_DIRS = {
    "naca0006": c.KURT5 / "runs" / "shape_spacing" / "naca0006",
    "naca0012": c.KURT5 / "runs" / "shape_spacing" / "naca0012_baseline",
    "naca0018": c.KURT5 / "runs" / "shape_spacing" / "naca0018",
}
NSTEPS = 3000


def a1():
    X, Y = c.grid_xy(DX)
    dx = DX
    cell_area = dx * dx
    rows = []
    for name, run_dir in RUN_DIRS.items():
        om = c.load_omega(run_dir, NSTEPS)
        Xw, Yw, sub, _, _ = c.window(X, Y, om, LE_XLIM, LE_YLIM)
        point_max = float(np.abs(sub).max())
        enstrophy = float(np.sum(sub ** 2) * cell_area)
        rms = float(np.sqrt(np.mean(sub ** 2)))
        extent = float(np.sum(np.abs(sub) > V_THRESHOLD) * cell_area)
        rows.append(dict(shape=name, point_max=point_max, enstrophy=enstrophy,
                          rms=rms, extent_above_thresh=extent))

    with open(c.DATA / "testA1_metric_robustness.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["shape", "point_max", "enstrophy", "rms", "extent_above_thresh"])
        w.writeheader(); w.writerows(rows)

    print("Test A1: LE quantity under 4 metrics (NACA0006/0012/0018)")
    for r in rows:
        print(f"  {r['shape']}: point_max={r['point_max']:.2f}, enstrophy={r['enstrophy']:.3f}, "
              f"rms={r['rms']:.3f}, extent(|w|>{V_THRESHOLD})={r['extent_above_thresh']:.4f}")

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    names = list(RUN_DIRS.keys())
    for ax, key, label in zip(axes,
                               ["point_max", "enstrophy", "rms", "extent_above_thresh"],
                               ["point max |omega|\n(original metric)", "enstrophy (sum(w^2)*dA)",
                                "RMS(omega) in window", f"area with |omega|>{V_THRESHOLD}"]):
        vals = [next(r[key] for r in rows if r["shape"] == n) for n in names]
        ax.bar(names, vals, color=["#2980b9", "#c0392b", "#16a085"])
        ax.set_title(label, fontsize=9)
        ax.tick_params(axis="x", rotation=20)
    fig.suptitle("Test A1: is the NACA0006/0012/0018 non-monotonicity robust across "
                 "4 different LE quantity metrics?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    out = c.FIGS / "testA1_metric_robustness.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out.name}")


def a2():
    rows = []
    for name, geom_path in SHAPES.items():
        g = c.load_geom_points(geom_path)
        x_le, y_le = g["x"][g["i_le"]], g["y"][g["i_le"]]
        # sub-cell phase: where does the LE sit within its grid cell, as a
        # fraction of dx, relative to the grid's own cell edges (grid_xy's
        # xoffset=-2.0, so cell edges are at -2.0 + k*dx)
        phase_x = ((x_le - (-2.0)) % DX) / DX
        phase_y = ((y_le - (-1.5)) % DX) / DX
        # residual asymmetry: mean |y| of points near x=0 above vs below,
        # as a check for accidental up/down bias in the resampled boundary
        near_le = np.abs(g["x"]) < 0.05
        y_near_le = g["y"][near_le]
        asym = float(np.mean(y_near_le)) if len(y_near_le) else float("nan")
        rows.append(dict(shape=name, le_x=x_le, le_y=y_le, phase_x=phase_x, phase_y=phase_y,
                          near_le_y_asymmetry=asym, n_points=g["n"]))

    with open(c.DATA / "testA2_geometry_audit.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["shape", "le_x", "le_y", "phase_x", "phase_y",
                                           "near_le_y_asymmetry", "n_points"])
        w.writeheader(); w.writerows(rows)

    print("\nTest A2: LE sub-cell phase audit (phase = fraction of one dx cell, 0=cell edge, 0.5=cell center)")
    for r in rows:
        print(f"  {r['shape']}: LE=({r['le_x']:.5f},{r['le_y']:.5f}), "
              f"phase_x={r['phase_x']:.3f}, phase_y={r['phase_y']:.3f}, "
              f"near-LE y-asymmetry={r['near_le_y_asymmetry']:.5f}")
    return rows


if __name__ == "__main__":
    a1()
    a2()
