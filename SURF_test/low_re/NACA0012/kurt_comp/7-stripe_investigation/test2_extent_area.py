"""
test2_extent_area.py

Test 2: how much of the domain is contaminated, not just how strong.
Threshold-area curves A(tau) (area where |omega|>tau) for tau in a small
ladder, plus total upstream enstrophy, integral |omega|, and the
signed-vs-absolute integral contrast (a coherent vortical structure would
carry net circulation; a sign-alternating checkerboard cancels -- if
signed << absolute, that's direct evidence of oscillatory noise, not a
structure). NACA0012, alpha=0, steady, Re=1000, dx=0.02/0.01/0.005, both
implementations. Zero new runs.

Usage: python3 test2_extent_area.py
Output: figures/test2_extent_area.png, data/test2_extent_area.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0
THRESHOLDS = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)

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
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    metrics = {}  # dx -> {impl: metrics dict}
    rows = []
    for case in CASES:
        dx = case["dx"]
        X, Y = c.grid_xy(dx)
        metrics[dx] = {}
        for impl in ("py", "cpp"):
            om = c.load_omega(case[impl], case["nsteps"])
            m = c.upstream_scalar_metrics(X, Y, om, x_le, dx, buffer_dx=BUFFER_DX,
                                            ylim=YLIM, thresholds=THRESHOLDS)
            metrics[dx][impl] = m
            row = dict(dx=dx, impl=impl, enstrophy=m["enstrophy"], int_abs=m["int_abs"],
                       int_signed=m["int_signed"], peak=m["peak"],
                       signed_over_abs=m["int_signed"] / m["int_abs"] if m["int_abs"] else 0.0)
            for tau in THRESHOLDS:
                row[f"area_tau{tau}"] = m["areas"][tau]
            rows.append(row)

    print("py vs cpp agreement (enstrophy, relative diff):")
    for dx in metrics:
        py_e = metrics[dx]["py"]["enstrophy"]
        cpp_e = metrics[dx]["cpp"]["enstrophy"]
        print(f"  dx={dx}: py={py_e:.5f} cpp={cpp_e:.5f} reldiff={abs(py_e-cpp_e)/py_e:.2e}")

    with open(c.DATA / "test2_extent_area.csv", "w", newline="") as f:
        fieldnames = ["dx", "impl", "enstrophy", "int_abs", "int_signed", "peak", "signed_over_abs"] + \
                     [f"area_tau{tau}" for tau in THRESHOLDS]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)
    print("wrote test2_extent_area.csv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # panel 1: threshold-area curves
    ax = axes[0]
    for dx in sorted(metrics, reverse=True):
        areas_py = [metrics[dx]["py"]["areas"][tau] for tau in THRESHOLDS]
        areas_cpp = [metrics[dx]["cpp"]["areas"][tau] for tau in THRESHOLDS]
        ax.plot(THRESHOLDS, areas_py, "o-", color=DX_COLOR[dx], label=f"dx={dx}", lw=1.8)
        ax.plot(THRESHOLDS, areas_cpp, "x--", color=DX_COLOR[dx], alpha=0.6, ms=5)
    ax.set_xlabel("threshold tau"); ax.set_ylabel("upstream area where |omega|>tau (chord^2)")
    ax.set_title("Threshold-area A(tau)", fontsize=10)
    ax.set_yscale("log"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    # panel 2: enstrophy and int_abs vs dx
    ax = axes[1]
    dxs = sorted(metrics)
    ens_py = [metrics[dx]["py"]["enstrophy"] for dx in dxs]
    ens_cpp = [metrics[dx]["cpp"]["enstrophy"] for dx in dxs]
    intabs_py = [metrics[dx]["py"]["int_abs"] for dx in dxs]
    ax.plot(dxs, ens_py, "o-", color="#8e44ad", label="enstrophy (py)", lw=1.8)
    ax.plot(dxs, ens_cpp, "x--", color="#8e44ad", alpha=0.6, ms=6, label="enstrophy (cpp)")
    ax.plot(dxs, intabs_py, "o-", color="#16a085", label="integral |omega| (py)", lw=1.8)
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xlabel("dx (log, refining -->)"); ax.set_title("Upstream enstrophy & |omega| integral vs dx", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    # panel 3: signed vs absolute integral (cancellation check)
    ax = axes[2]
    width = 0.25
    xpos = np.arange(len(dxs))
    abs_vals = [metrics[dx]["py"]["int_abs"] for dx in dxs]
    signed_vals = [metrics[dx]["py"]["int_signed"] for dx in dxs]
    ax.bar(xpos - width / 2, abs_vals, width, label="integral |omega| dA", color="#c0392b")
    ax.bar(xpos + width / 2, np.abs(signed_vals), width, label="|integral omega dA| (signed)", color="#2980b9")
    ax.set_xticks(xpos); ax.set_xticklabels([f"dx={dx}" for dx in dxs])
    ax.set_yscale("log")
    ax.set_title("Cancellation check: |signed| << |absolute|\nimplies oscillatory noise, not a coherent vortex", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both", axis="y")

    fig.suptitle("Test 2: extent and area of upstream vorticity contamination\n"
                 "NACA0012, alpha=0, steady, Re=1000 (solid/o=py_static, dashed/x=cpp_static)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = c.FIGS / "test2_extent_area.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
