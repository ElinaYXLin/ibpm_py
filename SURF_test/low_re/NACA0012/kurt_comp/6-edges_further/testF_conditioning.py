"""
testF_conditioning.py

Group F: is the LE peak actually tracking projection-matrix conditioning
rather than nose radius directly? A static linear-algebra diagnostic, no
time-stepping at all -- py_static's own CholeskySolver already builds the
dense projection matrix M = C*Ainv*B explicitly (computeMatrixM(), used
internally before factoring it); this script reuses that method directly
on each geometry and reports the condition number of M.

Usage: python3 testF_conditioning.py
Output: data/testF_conditioning.csv, figures/testF_conditioning.png
"""
import csv
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

sys.path.insert(0, str(c.REPO))
from py_static.base_flow import BaseFlow  # noqa: E402
from py_static.cholesky_solver import CholeskySolver  # noqa: E402
from py_static.geometry import Geometry  # noqa: E402
from py_static.grid import Grid  # noqa: E402
from py_static.navier_stokes_model import NavierStokesModel  # noqa: E402

RE = 1000.0
DX = 0.02
NX, NY = 300, 150
DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.5)
BETA = 0.01  # representative beta=(an+bn)*dt; only relative conditioning across
             # geometries matters here, so any fixed, consistent beta is fine

GEOMETRIES = {
    "naca0006": c.KURT5 / "geom" / "naca0006_dx0.0200.geom",
    "naca0012_baseline": c.BASE_GEOM_DX002,
    "naca0018": c.KURT5 / "geom" / "naca0018_dx0.0200.geom",
    "naca0012_LTEdense": c.KURT5 / "geom" / "naca0012_dx0.0200_LTEdense.geom",
    "naca0012_LTEsparse": c.KURT5 / "geom" / "naca0012_dx0.0200_LTEsparse.geom",
    "naca0012_roundTE": c.KURT5 / "geom" / "naca0012_dx0.0200_roundTE.geom",
    "cylinder": c.REPO / "SURF_test" / "vortall" / "3-grid_refine" / "geom" / "cylinder_dx0.0200.geom",
}

# Test 3b's already-measured LE peaks (py_static, 2-D field-max metric would
# be more consistent with this test's framing, but the original reported
# number -- the y=0 lineout peak -- is what's being asked "does this track
# conditioning") -- see 5-leading_edge/data/test3_3b_shape.csv
LE_PEAK = {
    "naca0006": 8.302, "naca0012_baseline": 22.176, "naca0018": 6.833,
    "naca0012_LTEdense": 5.802, "naca0012_LTEsparse": 30.992,
    "naca0012_roundTE": 1.213, "cylinder": 0.000,
}


def condition_number(geom_path):
    grid = Grid(NX, NY, 1, DOMAIN["length"], DOMAIN["xoffset"], DOMAIN["yoffset"], 0.0, 0.0)
    geom = Geometry(str(geom_path))
    q_potential = BaseFlow(grid, 1.0, 0.0)
    model = NavierStokesModel(grid, geom, RE, q_potential)
    model.init()
    solver = CholeskySolver(grid, model, BETA)
    n = model.getNumPoints()
    matrixM = np.zeros((2 * n, 2 * n), dtype=np.float64)
    solver.computeMatrixM(matrixM)
    cond = np.linalg.cond(matrixM)
    eigvals = np.linalg.eigvalsh(0.5 * (matrixM + matrixM.T))
    return float(cond), float(eigvals.min()), float(eigvals.max()), n


def main():
    rows = []
    for name, geom_path in GEOMETRIES.items():
        cond, emin, emax, n = condition_number(geom_path)
        rows.append(dict(geometry=name, n_points=n, cond_number=cond,
                          min_eig=emin, max_eig=emax, le_peak_lineout=LE_PEAK.get(name, float("nan"))))
        print(f"{name}: n={n}, cond(M)={cond:.4e}, min_eig={emin:.4e}, max_eig={emax:.4e}, "
              f"LE_peak(lineout)={LE_PEAK.get(name, float('nan'))}")

    with open(c.DATA / "testF_conditioning.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["geometry", "n_points", "cond_number", "min_eig",
                                           "max_eig", "le_peak_lineout"])
        w.writeheader(); w.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    names = [r["geometry"] for r in rows]
    conds = [r["cond_number"] for r in rows]
    peaks = [r["le_peak_lineout"] for r in rows]

    ax = axes[0]
    ax.bar(names, conds, color="#8e44ad")
    ax.set_yscale("log")
    ax.set_ylabel("condition number of M (log scale)")
    ax.set_title("Projection-matrix conditioning per geometry", fontsize=10)
    ax.tick_params(axis="x", rotation=30)

    ax = axes[1]
    ax.scatter(conds, peaks, color="#c0392b", s=60, zorder=3)
    for n, x, y in zip(names, conds, peaks):
        ax.annotate(n, (x, y), fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_xlabel("condition number of M (log scale)")
    ax.set_ylabel("LE peak |omega| (y=0 lineout metric)")
    ax.set_title("Does the LE peak track conditioning, or r_LE directly?", fontsize=10)
    ax.grid(alpha=0.3)

    fig.suptitle("Test F: projection-matrix (M = C*Ainv*B) conditioning across geometries", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = c.FIGS / "testF_conditioning.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nwrote {out.name}")


if __name__ == "__main__":
    main()
