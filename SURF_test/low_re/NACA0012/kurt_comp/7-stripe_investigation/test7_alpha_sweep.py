"""
test7_alpha_sweep.py

Test 7: angle-of-attack sweep. This solver imposes alpha by rotating the
free-stream rather than the body (see ../1-paper_based/README.md's "Wake
vorticity fields" section), so the body-to-grid discretization is
IDENTICAL across the whole alpha=0..40,50,60 sweep that already exists
at dx=0.02 -- an unusually clean controlled experiment: alpha changes
the flow and boundary-force magnitude while holding the geometric
discretization exactly fixed. If upstream noise varies strongly with
alpha, it's driven by force magnitude; if flat, it's purely geometric.

Definitional note: "upstream" is fixed here as x < x_LE in the solver's
own (unrotated-body) frame, same region used throughout this folder --
NOT rotated to track the incoming free-stream direction. At alpha!=0
that region is no longer directly "ahead of" the oncoming flow in the
lab sense; it is still the same fixed geometric region the rest of this
investigation calls "upstream", so the comparison across alpha is
internally consistent even though its physical interpretation shifts.
Only alpha=0,9,12 turned out to have actual precomputed output in
../1-paper_based despite many more `steady_{py,cpp}_aXX` directory names
existing there (empty placeholders, never actually run) -- so this test
is NOT zero-new-runs as originally scoped; it launches the missing
alphas itself (cheap: dx=0.02, 300x150 grid, 3000 steps, same convention
as the existing alpha=0/9/12 cases, seconds each -- nothing like the
faithful2 runs occupying the machine elsewhere). New runs land in
./runs/alpha_sweep/, existing 0/9/12 are reused from ../1-paper_based.

Usage:
  python3 test7_alpha_sweep.py run       # launch missing alphas (blocking)
  python3 test7_alpha_sweep.py analyze   # analyze whatever exists (default)
Output: figures/test7_alpha_sweep.png, data/test7_alpha_sweep.csv
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0
DX = 0.02
NSTEPS = 3000
EXISTING_ALPHAS = [0, 9, 12]
NEW_ALPHAS = [3, 6, 15, 18, 21, 24, 27, 30, 33, 36, 40, 50, 60]
ALPHAS = EXISTING_ALPHAS + NEW_ALPHAS
NX, NY = 300, 150


def run_dir_for(alpha, impl):
    if alpha in EXISTING_ALPHAS:
        return c.KURT1 / "runs" / "dx0.020" / f"steady_{impl}_a{alpha:02d}"
    return c.RUNS / "alpha_sweep" / f"a{alpha:02d}_{impl}"


def do_run():
    for alpha in NEW_ALPHAS:
        for impl in ("py", "cpp"):
            outdir = run_dir_for(alpha, impl)
            if c.is_done(outdir, NSTEPS):
                print(f"alpha={alpha} {impl}: already done")
                continue
            print(f"alpha={alpha} {impl}: launching...", flush=True)
            ok, elapsed = c.run_case(impl, c.BASE_GEOM_DX002, outdir, NX, NY, 0.01, NSTEPS,
                                       alpha=float(alpha), restart=250)
            print(f"  {'OK' if ok else 'FAILED'} in {elapsed:.1f}s")


def main():
    X, Y = c.grid_xy(DX)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    rows = []
    for alpha in ALPHAS:
        for impl in ("py", "cpp"):
            run_dir = run_dir_for(alpha, impl)
            if not c.is_done(run_dir, NSTEPS):
                print(f"alpha={alpha} {impl}: not done, skipping")
                continue
            om = c.load_omega(run_dir, NSTEPS)
            m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            rows.append(dict(alpha=alpha, impl=impl, enstrophy=m["enstrophy"],
                              int_abs=m["int_abs"], peak=m["peak"]))

    with open(c.DATA / "test7_alpha_sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["alpha", "impl", "enstrophy", "int_abs", "peak"])
        w.writeheader(); w.writerows(rows)
    print("wrote test7_alpha_sweep.csv")

    py_rows = sorted([r for r in rows if r["impl"] == "py"], key=lambda r: r["alpha"])
    cpp_rows = sorted([r for r in rows if r["impl"] == "cpp"], key=lambda r: r["alpha"])
    py_ens = np.array([r["enstrophy"] for r in py_rows])
    cpp_ens = np.array([r["enstrophy"] for r in cpp_rows])
    reldiff = np.abs(py_ens - cpp_ens) / np.maximum(py_ens, 1e-12)
    print(f"py vs cpp: max reldiff over sweep = {reldiff.max():.2e}, mean = {reldiff.mean():.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    ax.plot([r["alpha"] for r in py_rows], py_ens, "o-", color="#1f77b4", label="py_static", lw=1.6, ms=4)
    ax.plot([r["alpha"] for r in cpp_rows], cpp_ens, "x--", color="#d62728", label="cpp_static",
             lw=1.0, ms=5, alpha=0.7)
    ax.set_yscale("log")
    ax.set_xlabel("alpha (degrees)"); ax.set_ylabel("upstream enstrophy")
    ax.set_title("Upstream enstrophy vs angle of attack\n(fixed geometric discretization throughout)", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    ax = axes[1]
    py_peak = [r["peak"] for r in py_rows]
    ax.plot([r["alpha"] for r in py_rows], py_peak, "o-", color="#8e44ad", lw=1.6, ms=4)
    ax.set_xlabel("alpha (degrees)"); ax.set_ylabel("upstream peak |omega|")
    ax.set_title("Upstream peak vs angle of attack", fontsize=10)
    ax.grid(alpha=0.3)

    fig.suptitle("Test 7: angle-of-attack sweep -- force-magnitude-driven or purely geometric?\n"
                 f"NACA0012, Re=1000, steady, dx=0.02, alpha=0..40,50,60 (py/cpp max reldiff={reldiff.max():.1e})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "test7_alpha_sweep.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        do_run()
    else:
        main()
