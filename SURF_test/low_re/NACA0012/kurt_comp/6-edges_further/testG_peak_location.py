"""
testG_peak_location.py

Group G: where, physically, is the reported "peak |omega|" actually
sitting -- inside the solid body (not physically meaningful; this
solver's Cartesian grid is not body-fitted, so vorticity is computed
at every grid node including ones geometrically inside the airfoil,
see README's Group G section), on the body surface / in the
boundary-layer artifact region the whole 5-leading_edge/6-edges_further
investigation is about, or out in an ordinary part of the wake?

Zero new runs -- reuses existing dx=0.02 (../1-paper_based) and dx=0.005
(../5-leading_edge/runs/grid_refine) py_static fields and geometry files.

For each case, both the y=0 lineout metric (5-leading_edge's original,
Reconciliation 1 found unreliable) and the 2-D field-max metric (the one
Reconciliation 1/2 and Group C treat as the real signal) are recomputed,
their (x,y) locations marked directly on the field, and each location is
classified inside/outside the body via a point-in-polygon test against
the actual boundary geometry -- so the two metrics' disagreement, and
whether either one is (wrongly) picking up interior grid noise, is
visible directly rather than inferred.

Usage: python3 testG_peak_location.py
Output: figures/testG_peak_location.png (per dx: overview + LE zoom + TE
        zoom, both metrics marked), data/testG_peak_location.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

import common as c

LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
TE_XLIM, TE_YLIM = (0.65, 1.35), (-0.25, 0.25)
OVERVIEW_XLIM, OVERVIEW_YLIM = (-0.3, 1.3), (-0.35, 0.35)
V = 8.0  # colormap clip, matches 5-leading_edge's test2a convention

CASES = [
    dict(dx=0.02, impl="py", label="dx=0.02",
         run_dir=c.KURT / "1-paper_based" / "runs" / "dx0.020" / "steady_py_a00", step=3000,
         geom=c.BASE_GEOM_DX002),
    dict(dx=0.005, impl="py", label="dx=0.005",
         run_dir=c.KURT5 / "runs" / "grid_refine" / "dx0.0050", step=12000,
         geom=c.REPO / "SURF_test" / "geom" / "naca0012_dx0.0050.geom"),
]


def lineout_peak_loc(X, Y, om, xlim, ylim):
    """y=0 lineout metric (5-leading_edge's original): max |omega| along
    the single grid row nearest y=0, restricted to xlim. Returns (x, y,
    value)."""
    ys = Y[0, :]
    iy0 = c.nearest_index(ys, 0.0)
    xs = X[:, 0]
    m = (xs >= xlim[0]) & (xs <= xlim[1])
    lineout = om[:, iy0][m]
    xs_m = xs[m]
    k = int(np.argmax(np.abs(lineout)))
    return float(xs_m[k]), float(ys[iy0]), float(lineout[k])


def field_max_loc(X, Y, om, xlim, ylim):
    """2-D field-max metric (Reconciliation 1/2, Group C): max |omega|
    anywhere in the (xlim, ylim) window. Returns (x, y, value)."""
    Xw, Yw, fw, ix, iy = c.window(X, Y, om, xlim, ylim)
    idx = np.unravel_index(np.argmax(np.abs(fw)), fw.shape)
    return float(Xw[idx]), float(Yw[idx]), float(fw[idx])


def inside_body(x, y, poly_x, poly_y):
    verts = np.column_stack([poly_x, poly_y])
    return bool(MplPath(verts).contains_point((x, y)))


def draw_panel(ax, X, Y, om, xlim, ylim, poly_x, poly_y, markers, title):
    Xw, Yw, fw, _, _ = c.window(X, Y, om, xlim, ylim)
    ax.pcolormesh(Xw, Yw, np.clip(fw[:-1, :-1], -V, V), shading="flat", cmap="jet", vmin=-V, vmax=V)
    # closed body outline, drawn on top of the field
    ax.plot(np.append(poly_x, poly_x[0]), np.append(poly_y, poly_y[0]),
             color="white", lw=2.2, zorder=5)
    ax.plot(np.append(poly_x, poly_x[0]), np.append(poly_y, poly_y[0]),
             color="black", lw=0.9, zorder=6)
    for (mx, my, mlabel, mcolor, mmarker) in markers:
        ax.scatter([mx], [my], s=90, facecolor=mcolor, edgecolor="black",
                   marker=mmarker, linewidth=1.0, zorder=6, label=mlabel)
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=8)


def main():
    rows = []
    n_cases = len(CASES)
    fig, axes = plt.subplots(n_cases, 3, figsize=(13, 4.3 * n_cases), squeeze=False)

    for row_i, case in enumerate(CASES):
        X, Y = c.grid_xy(case["dx"])
        om = c.load_omega(case["run_dir"], case["step"])
        geom_pts = c.load_geom_points(case["geom"])
        poly_x, poly_y = geom_pts["x"], geom_pts["y"]

        case_markers = {"LE": [], "TE": []}
        for region, xlim, ylim in [("LE", LE_XLIM, LE_YLIM), ("TE", TE_XLIM, TE_YLIM)]:
            lx, ly, lval = lineout_peak_loc(X, Y, om, xlim, ylim)
            fx, fy, fval = field_max_loc(X, Y, om, xlim, ylim)
            l_inside = inside_body(lx, ly, poly_x, poly_y)
            f_inside = inside_body(fx, fy, poly_x, poly_y)
            rows.append(dict(dx=case["dx"], region=region, metric="lineout",
                              x=lx, y=ly, value=lval, inside_body=l_inside))
            rows.append(dict(dx=case["dx"], region=region, metric="field_max",
                              x=fx, y=fy, value=fval, inside_body=f_inside))
            case_markers[region] = [
                (lx, ly, f"lineout peak ({'IN' if l_inside else 'OUT'})", "white", "o"),
                (fx, fy, f"field-max peak ({'IN' if f_inside else 'OUT'})", "lime", "^"),
            ]
            print(f"dx={case['dx']} {region}: lineout peak at ({lx:.4f},{ly:.4f})="
                  f"{lval:.2f} [{'inside body' if l_inside else 'in fluid'}];  "
                  f"field-max peak at ({fx:.4f},{fy:.4f})={fval:.2f} "
                  f"[{'inside body' if f_inside else 'in fluid'}]")

        all_markers = case_markers["LE"] + case_markers["TE"]
        draw_panel(axes[row_i, 0], X, Y, om, OVERVIEW_XLIM, OVERVIEW_YLIM, poly_x, poly_y,
                   all_markers, f"{case['label']}: overview (entire airfoil)")
        draw_panel(axes[row_i, 1], X, Y, om, LE_XLIM, LE_YLIM, poly_x, poly_y,
                   case_markers["LE"], f"{case['label']}: LE zoom")
        draw_panel(axes[row_i, 2], X, Y, om, TE_XLIM, TE_YLIM, poly_x, poly_y,
                   case_markers["TE"], f"{case['label']}: TE zoom")
        axes[row_i, 1].legend(fontsize=6, loc="upper right")

    fig.suptitle("Test G: where is the reported peak |omega| actually located?\n"
                 "black outline = actual body surface; white circle = y=0 lineout peak; "
                 "lime triangle = 2-D field-max peak (IN=inside body, OUT=in fluid)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = c.FIGS / "testG_peak_location.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\nwrote {out.name}")

    with open(c.DATA / "testG_peak_location.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dx", "region", "metric", "x", "y", "value", "inside_body"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote testG_peak_location.csv")


if __name__ == "__main__":
    main()
