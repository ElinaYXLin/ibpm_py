"""
test2_density_conditioning_extended.py

Follow-up test #2 from ../7-stripe_investigation/README.md's Proposals
section: ../7-stripe_investigation Test 6 found condition number and
upstream enstrophy monotonic together, but only across 3 geometries
(LTEsparse/baseline/LTEdense). ../6-edges_further Group D's existing
LE-only density sweep (0.5x, 2x, 4x=recon2's LEonly_dense, 8x, and the
diverged 16x) already has known flow output; this recomputes upstream
enstrophy for those geometries and builds each one's projection-matrix
condition number fresh (Group F never computed conditioning for the D1
geometries, only for LTEsparse/baseline/LTEdense/shape family), turning
Test 6's 3-point relationship into a 6-7-point dose-response curve
spanning ~5 orders of magnitude in cond(M). Zero new runs (existing D1
flow output); the new work is building 6 projection matrices, each
cheap (101-112 boundary points, same order as Group F's own baseline
matrix).

Usage: python3 test2_density_conditioning_extended.py
Output: figures/test2_density_conditioning_extended.png,
        data/test2_density_conditioning_extended.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

DX = 0.02
BUFFER_DX = 2.0
YLIM = (-0.5, 0.5)
NSTEPS = 3000
GEOMDIR = c.KURT6 / "geom"

CASES = {
    0.5: dict(geom=GEOMDIR / "naca0012_dx0.0200_LEdensity0.5x.geom",
              run=c.KURT6 / "runs" / "D1_point_density" / "LEdensity0.5x"),
    1.0: dict(geom=c.BASE_GEOM_DX002,
              run=c.KURT5 / "runs" / "shape_spacing" / "naca0012_baseline"),
    2.0: dict(geom=GEOMDIR / "naca0012_dx0.0200_LEdensity2x.geom",
              run=c.KURT6 / "runs" / "D1_point_density" / "LEdensity2x"),
    4.0: dict(geom=GEOMDIR / "naca0012_dx0.0200_LEonly_dense.geom",
              run=c.KURT6 / "runs" / "recon2" / "LEonly_dense_Re1000"),
    8.0: dict(geom=GEOMDIR / "naca0012_dx0.0200_LEdensity8x.geom",
              run=c.KURT6 / "runs" / "D1_point_density" / "LEdensity8x"),
    16.0: dict(geom=GEOMDIR / "naca0012_dx0.0200_LEdensity16x.geom",
               run=c.KURT6 / "runs" / "D1_point_density" / "LEdensity16x"),
}


def main():
    X, Y = c.grid_xy(DX)
    grid = c.Grid(300, 150, 1, c.DOMAIN["length"], c.DOMAIN["xoffset"], c.DOMAIN["yoffset"])

    rows = []
    for factor, case in CASES.items():
        geom = c.Geometry(str(case["geom"]))
        M, n = c.build_projection_matrix(grid, geom, c.RE)
        eigvals = np.linalg.eigvalsh(M)
        # at extreme ill-conditioning (factor=16), roundoff can push the
        # numerically smallest eigenvalue of a theoretically SPD matrix
        # slightly negative -- cond as max/min then comes out negative,
        # which isn't a meaningful "negative condition number", it's a
        # symptom of just how close to singular M already is. Report
        # magnitude and flag it rather than let a negative value silently
        # break the log-scale plot.
        raw_cond = float(eigvals.max() / eigvals.min())
        cond = abs(raw_cond)
        eig_min_negative = bool(eigvals.min() < 0)

        om = c.load_omega(case["run"], NSTEPS)
        diverged = bool(np.any(~np.isfinite(om)))
        if diverged:
            enstrophy = float("nan")
            flag = " (min eigenvalue of M is numerically NEGATIVE -- matrix has lost SPD-ness)" if eig_min_negative else ""
            print(f"factor={factor}: n={n}, |cond|={cond:.3e}{flag}, DIVERGED (NaN/Inf in field)")
        else:
            g = c.load_geom_points(case["geom"])
            x_le = g["x"][g["i_le"]]
            m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            enstrophy = m["enstrophy"]
            print(f"factor={factor}: n={n}, cond={cond:.3e}, upstream_enstrophy={enstrophy:.5f}")

        rows.append(dict(density_factor=factor, n_points=n, cond=cond, raw_cond=raw_cond,
                          eig_min_negative=eig_min_negative, enstrophy=enstrophy, diverged=diverged))

    with open(c.DATA / "test2_density_conditioning_extended.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["density_factor", "n_points", "cond", "raw_cond",
                                            "eig_min_negative", "enstrophy", "diverged"])
        w.writeheader(); w.writerows(rows)
    print("wrote test2_density_conditioning_extended.csv")

    fig, ax = plt.subplots(figsize=(8, 6))
    finite_rows = [r for r in rows if not r["diverged"]]
    ax.plot([r["cond"] for r in finite_rows], [r["enstrophy"] for r in finite_rows],
            "o-", color="#2980b9", lw=1.8, ms=8, zorder=5)
    for r in finite_rows:
        ax.annotate(f"{r['density_factor']:g}x", (r["cond"], r["enstrophy"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=8)
    diverged_rows = [r for r in rows if r["diverged"]]
    for r in diverged_rows:
        ax.axvline(r["cond"], color="#c0392b", ls="--", lw=1.2, alpha=0.7)
        ax.annotate(f"{r['density_factor']:g}x DIVERGED\n(cond={r['cond']:.2e})",
                    (r["cond"], 0), textcoords="offset points", xytext=(-60, 20),
                    fontsize=8, color="#c0392b")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("condition number of projection matrix M")
    ax.set_ylabel("upstream enstrophy")
    ax.set_title("Test 2: conditioning vs. upstream enstrophy, extended to 6 LE-density levels\n"
                 "(Test 6 in ../7-stripe_investigation had only 3 points)", fontsize=11)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = c.FIGS / "test2_density_conditioning_extended.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
