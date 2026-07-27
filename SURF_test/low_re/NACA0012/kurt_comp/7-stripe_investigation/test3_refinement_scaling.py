"""
test3_refinement_scaling.py

Test 3: refinement scaling, measured in two units. Fits upstream
enstrophy vs dx to a power law (enstrophy ~ dx^p): p>0 and enstrophy->0
means a consistent, converging scheme; a plateau means the artifact is
persistent (much more serious). Then reports L_up (the reach of the
decay profile out to a fixed noise floor, from Test 1's per-column
max|omega| profile) in BOTH chord units and cells (L_up/dx) -- constant
in cells implicates the discrete delta function's stencil support
(3-4 cells, a local mechanism); constant in chord implicates something
non-local. NACA0012, alpha=0, steady, Re=1000, both implementations.
Zero new runs (reuses Test 1's profiles).

Usage: python3 test3_refinement_scaling.py
Output: figures/test3_refinement_scaling.png, data/test3_refinement_scaling.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0
NOISE_FLOOR = 1e-3  # well above float noise (~1e-13*peak), well below any real LE signal

CASES = [
    dict(dx=0.02, py=c.KURT1 / "runs" / "dx0.020" / "steady_py_a00",
         cpp=c.KURT1 / "runs" / "dx0.020" / "steady_cpp_a00", nsteps=3000),
    dict(dx=0.01, py=c.KURT5 / "runs" / "grid_refine" / "dx0.0100",
         cpp=c.KURT5 / "runs" / "grid_refine" / "dx0.0100_cpp", nsteps=6000),
    dict(dx=0.005, py=c.KURT5 / "runs" / "grid_refine" / "dx0.0050",
         cpp=c.KURT5 / "runs" / "grid_refine" / "dx0.0050_cpp", nsteps=12000),
]


def main():
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    rows = []
    for case in CASES:
        dx = case["dx"]
        X, Y = c.grid_xy(dx)
        for impl in ("py", "cpp"):
            om = c.load_omega(case[impl], case["nsteps"])
            m = c.upstream_scalar_metrics(X, Y, om, x_le, dx, buffer_dx=BUFFER_DX, ylim=YLIM)
            prof = c.upstream_profile(X, Y, om, x_le, dx, buffer_dx=BUFFER_DX, ylim=YLIM)
            L_up_chord, _, _ = c.reach_L_up(prof["xs"], prof["max_abs"], x_le, NOISE_FLOOR)
            rows.append(dict(dx=dx, impl=impl, enstrophy=m["enstrophy"],
                              L_up_chord=L_up_chord, L_up_cells=L_up_chord / dx))

    with open(c.DATA / "test3_refinement_scaling.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dx", "impl", "enstrophy", "L_up_chord", "L_up_cells"])
        w.writeheader(); w.writerows(rows)
    print("wrote test3_refinement_scaling.csv")
    print(f"\n(noise floor = {NOISE_FLOOR} used for L_up)")
    for r in rows:
        print(f"  dx={r['dx']} {r['impl']}: enstrophy={r['enstrophy']:.5f}, "
              f"L_up={r['L_up_chord']:.4f}c = {r['L_up_cells']:.2f} cells")

    # power-law fit: enstrophy ~ dx^p, using py data (3 points)
    py_rows = sorted([r for r in rows if r["impl"] == "py"], key=lambda r: r["dx"])
    dxs = np.array([r["dx"] for r in py_rows])
    ens = np.array([r["enstrophy"] for r in py_rows])
    log_dx = np.log(dxs)
    log_ens = np.log(np.maximum(ens, 1e-12))
    p, log_c = np.polyfit(log_dx, log_ens, 1)
    print(f"\npower-law fit (py, 3 points): enstrophy ~ dx^{p:.2f}  (const={np.exp(log_c):.4g})")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    for impl, marker in (("py", "o"), ("cpp", "x")):
        sub = sorted([r for r in rows if r["impl"] == impl], key=lambda r: r["dx"])
        ax.plot([r["dx"] for r in sub], [max(r["enstrophy"], 1e-8) for r in sub],
                marker + "-", label=f"enstrophy ({impl})", lw=1.8 if impl == "py" else 1.0,
                alpha=1.0 if impl == "py" else 0.6, ms=8)
    fit_line = np.exp(log_c) * dxs ** p
    ax.plot(dxs, fit_line, "k--", lw=1.2, label=f"fit: dx^{p:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("dx (log, refining -->)"); ax.set_ylabel("upstream enstrophy")
    ax.set_title("Enstrophy vs dx: power-law fit", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    ax2 = axes[1]
    ax2b = ax2.twinx()
    py_sub = sorted([r for r in rows if r["impl"] == "py"], key=lambda r: r["dx"])
    dxs_p = [r["dx"] for r in py_sub]
    l1 = ax2.plot(dxs_p, [r["L_up_chord"] for r in py_sub], "o-", color="#c0392b",
                   label="L_up (chord units)", lw=1.8, ms=8)
    l2 = ax2b.plot(dxs_p, [r["L_up_cells"] for r in py_sub], "s-", color="#2980b9",
                    label="L_up (cells = L_up/dx)", lw=1.8, ms=8)
    ax2.set_xscale("log"); ax2.invert_xaxis()
    ax2.set_xlabel("dx (log, refining -->)")
    ax2.set_ylabel("L_up (chord units)", color="#c0392b")
    ax2b.set_ylabel("L_up (cells)", color="#2980b9")
    ax2.set_title(f"Reach L_up (|omega|>{NOISE_FLOOR} floor):\nflat in chord -> non-local; "
                   "flat in cells -> local (delta-fn stencil)", fontsize=9)
    lines = l1 + l2
    ax2.legend(lines, [l.get_label() for l in lines], fontsize=8)
    ax2.grid(alpha=0.3, which="both")

    fig.suptitle("Test 3: refinement scaling of upstream contamination -- rate and mechanism", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = c.FIGS / "test3_refinement_scaling.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
