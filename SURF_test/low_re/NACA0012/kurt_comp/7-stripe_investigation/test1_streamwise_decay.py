"""
test1_streamwise_decay.py

Test 1: streamwise decay profile of the upstream (x < x_LE) striping.
For each x-column upstream of the LE, over all y in a window, compute
max|omega|, RMS, and integral(|omega| dy), and plot vs distance upstream
on a log-y axis, overlaying dx=0.02/0.01/0.005 -- NACA0012, alpha=0,
steady, Re=1000. py_static and cpp_static both plotted (solid/dashed) to
confirm agreement. Zero new runs -- reuses ../1-paper_based (dx=0.02) and
../5-leading_edge/runs/grid_refine (dx=0.01/0.005).

Usage: python3 test1_streamwise_decay.py
Output: figures/test1_streamwise_decay.png, data/test1_streamwise_decay.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0  # exclude 2 cells nearest the LE (genuine wraparound boundary layer)

CASES = [
    dict(dx=0.02, py=c.KURT1 / "runs" / "dx0.020" / "steady_py_a00",
         cpp=c.KURT1 / "runs" / "dx0.020" / "steady_cpp_a00", nsteps=3000),
    dict(dx=0.01, py=c.KURT5 / "runs" / "grid_refine" / "dx0.0100",
         cpp=c.KURT5 / "runs" / "grid_refine" / "dx0.0100_cpp", nsteps=6000),
    dict(dx=0.005, py=c.KURT5 / "runs" / "grid_refine" / "dx0.0050",
         cpp=c.KURT5 / "runs" / "grid_refine" / "dx0.0050_cpp", nsteps=12000),
]
DX_COLOR = {0.02: "#1f77b4", 0.01: "#e67e22", 0.005: "#27ae60"}


def main():
    x_le = 0.0  # NACA0012 raw coords: LE at x~0 (see common.load_geom_points on BASE_GEOM_DX002)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    profiles = {}  # dx -> {impl: profile dict}
    rows = []
    for case in CASES:
        dx = case["dx"]
        X, Y = c.grid_xy(dx)
        profiles[dx] = {}
        for impl in ("py", "cpp"):
            om = c.load_omega(case[impl], case["nsteps"])
            prof = c.upstream_profile(X, Y, om, x_le, dx, buffer_dx=BUFFER_DX, ylim=YLIM)
            profiles[dx][impl] = prof
            for xi, ma, r, ia, isg in zip(prof["xs"], prof["max_abs"], prof["rms"],
                                           prof["int_abs_dy"], prof["int_signed_dy"]):
                rows.append(dict(dx=dx, impl=impl, x=xi, dist_upstream=x_le - xi,
                                  max_abs=ma, rms=r, int_abs_dy=ia, int_signed_dy=isg))

    # py vs cpp agreement check: max relative difference in max_abs profile, per dx
    print("py vs cpp agreement (max relative diff in max|omega| profile):")
    agree_lines = []
    for dx in profiles:
        py_v = profiles[dx]["py"]["max_abs"]
        cpp_v = profiles[dx]["cpp"]["max_abs"]
        reldiff = np.abs(py_v - cpp_v) / np.maximum(np.abs(py_v), 1e-12)
        line = f"dx={dx}: max reldiff={reldiff.max():.2e}, mean reldiff={reldiff.mean():.2e}"
        print(f"  {line}")
        agree_lines.append(line)

    with open(c.DATA / "test1_streamwise_decay.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dx", "impl", "x", "dist_upstream", "max_abs", "rms",
                                           "int_abs_dy", "int_signed_dy"])
        w.writeheader(); w.writerows(rows)
    print("wrote test1_streamwise_decay.csv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, key, title in [(axes[0], "max_abs", "max|omega| over y"),
                             (axes[1], "rms", "RMS(omega) over y"),
                             (axes[2], "int_abs_dy", "integral |omega| dy")]:
        for dx in sorted(profiles, reverse=True):
            for impl, ls in (("py", "-"), ("cpp", "--")):
                prof = profiles[dx][impl]
                order = np.argsort(prof["xs"])
                dist = x_le - prof["xs"][order]
                vals = np.clip(prof[key][order], 1e-6, None)
                lbl = f"dx={dx} {impl}" if impl == "py" else None
                ax.plot(dist, vals, ls, color=DX_COLOR[dx], lw=1.6 if impl == "py" else 1.0,
                         alpha=1.0 if impl == "py" else 0.6, label=lbl)
        ax.set_yscale("log")
        ax.set_xlabel("distance upstream of LE (chord)")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.suptitle("Test 1: streamwise decay of upstream (x<x_LE) vorticity noise\n"
                 "NACA0012, alpha=0, steady, Re=1000 -- solid=py_static, dashed=cpp_static "
                 "(dashed overlaps solid almost exactly -- see agreement box)", fontsize=11)
    fig.text(0.01, 0.90, "py vs cpp agreement (max|omega| profile):\n" + "\n".join(agree_lines),
              fontsize=8, family="monospace", va="top")
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "test1_streamwise_decay.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
