"""
test6_spacing_conditioning.py

Test 6: boundary-point spacing and conditioning, retested on the clean
upstream ground. `../6-edges_further` Group F found the projection
matrix's condition number does NOT track the LE peak -- but that
comparison used the LE peak, which mixes real boundary-layer physics
with numerical artifact (per Group G). Conditioning is a pure
numerical-error quantity, so if it predicts anything, it should predict
a pure numerical-error signal -- upstream enstrophy is exactly that,
with no physics mixed in. Recomputes upstream metrics for the LTEsparse
/ baseline / LTEdense geometry variants (dx=0.02, Re=1000, alpha=0,
steady) and scatters them against Group F's existing condition numbers.
Zero new runs -- reuses ../5-leading_edge/runs/shape_spacing.

Usage: python3 test6_spacing_conditioning.py
Output: figures/test6_spacing_conditioning.png, data/test6_spacing_conditioning.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0
DX = 0.02
NSTEPS = 3000

CASES = {
    "naca0012_LTEsparse": dict(geom=c.KURT5 / "geom" / "naca0012_dx0.0200_LTEsparse.geom",
                                py=c.KURT5 / "runs" / "shape_spacing" / "naca0012_LTEsparse",
                                cpp=c.KURT5 / "runs" / "shape_spacing" / "naca0012_LTEsparse_cpp",
                                cond=4049.656866563488),
    "naca0012_baseline": dict(geom=c.BASE_GEOM_DX002,
                               py=c.KURT5 / "runs" / "shape_spacing" / "naca0012_baseline",
                               cpp=c.KURT5 / "runs" / "shape_spacing" / "naca0012_baseline_cpp",
                               cond=12708.490686345693),
    "naca0012_LTEdense": dict(geom=c.KURT5 / "geom" / "naca0012_dx0.0200_LTEdense.geom",
                               py=c.KURT5 / "runs" / "shape_spacing" / "naca0012_LTEdense",
                               cpp=c.KURT5 / "runs" / "shape_spacing" / "naca0012_LTEdense_cpp",
                               cond=114751269.4259997),
}
CASE_COLOR = {"naca0012_LTEsparse": "#27ae60", "naca0012_baseline": "#1f77b4", "naca0012_LTEdense": "#c0392b"}


def main():
    X, Y = c.grid_xy(DX)
    rows = []
    for name, case in CASES.items():
        g = c.load_geom_points(case["geom"])
        x_le = g["x"][g["i_le"]]
        for impl in ("py", "cpp"):
            om = c.load_omega(case[impl], NSTEPS)
            m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            rows.append(dict(case=name, impl=impl, cond=case["cond"], enstrophy=m["enstrophy"],
                              int_abs=m["int_abs"], peak=m["peak"]))
            print(f"{name} {impl}: cond={case['cond']:.3e} upstream_enstrophy={m['enstrophy']:.6f} "
                  f"upstream_peak={m['peak']:.3f}")

    with open(c.DATA / "test6_spacing_conditioning.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "impl", "cond", "enstrophy", "int_abs", "peak"])
        w.writeheader(); w.writerows(rows)
    print("wrote test6_spacing_conditioning.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    for name in CASES:
        py_r = next(r for r in rows if r["case"] == name and r["impl"] == "py")
        cpp_r = next(r for r in rows if r["case"] == name and r["impl"] == "cpp")
        ax.scatter([py_r["cond"]], [py_r["enstrophy"]], s=90, color=CASE_COLOR[name],
                   marker="o", label=f"{name} (py)", zorder=5)
        ax.scatter([cpp_r["cond"]], [cpp_r["enstrophy"]], s=60, color=CASE_COLOR[name],
                   marker="x", zorder=5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("condition number of projection matrix M (Group F)")
    ax.set_ylabel("upstream enstrophy")
    ax.set_title("Does conditioning predict upstream (pure-artifact) enstrophy?", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7, loc="upper left")
    # explicit minor-tick labels (e.g. 0.2, 0.5, 2) so the axis direction/scale
    # is legible without having to infer it from just the decade labels --
    # matplotlib hides minor tick labels on log axes by default even with a
    # locator/formatter set, so set values+labels directly instead
    ymin, ymax = ax.get_ylim()
    minor_vals = [v for v in (0.2, 0.3, 0.5, 2, 3, 5) if ymin <= v <= ymax]
    ax.set_yticks(minor_vals, minor=True)
    ax.set_yticklabels([f"{v:g}" for v in minor_vals], minor=True, fontsize=7)

    ax = axes[1]
    names = list(CASES)
    xpos = np.arange(len(names))
    conds = [CASES[n]["cond"] for n in names]
    ens = [next(r["enstrophy"] for r in rows if r["case"] == n and r["impl"] == "py") for n in names]
    ax2 = ax.twinx()
    b1 = ax.bar(xpos - 0.2, conds, 0.4, color="#7f8c8d", label="condition number")
    b2 = ax2.bar(xpos + 0.2, ens, 0.4, color="#e67e22", label="upstream enstrophy")
    ax.set_yscale("log"); ax2.set_yscale("log")
    ax.set_xticks(xpos); ax.set_xticklabels(names, fontsize=8, rotation=15)
    ax.set_ylabel("condition number", color="#7f8c8d")
    ax2.set_ylabel("upstream enstrophy", color="#e67e22")
    ax.set_title("Side-by-side: does monotonic conditioning match monotonic enstrophy?", fontsize=10)

    fig.suptitle("Test 6: boundary-point spacing & conditioning vs. upstream (clean) noise\n"
                 "NACA0012 LTEsparse/baseline/LTEdense, Re=1000, alpha=0, steady, dx=0.02", fontsize=11)
    fig.tight_layout(rect=[0.02, 0, 1, 0.88])
    out = c.FIGS / "test6_spacing_conditioning.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
