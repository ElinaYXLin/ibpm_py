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

Definitional note, and the lab-frame follow-up: "upstream" is fixed
throughout most of this folder as x < x_LE in the solver's own
(unrotated-body) frame -- NOT rotated to track the incoming free-stream
direction. At alpha!=0 that region is no longer directly "ahead of" the
oncoming flow in the lab sense. This script now computes BOTH versions
at every alpha: the fixed body-frame region (as before) and a
lab-frame region that rotates with alpha to stay aligned with the
actual oncoming flow direction (`common.upstream_mask_2d_rotated`,
using the free-stream unit vector (cos(alpha),sin(alpha)) implied by
py_static/ibpm.py's own drag/lift rotation convention). Comparing the
two answers the natural follow-up question directly: does "how much
upstream noise" depend on which region you call upstream?

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
            ens_rot = c.upstream_enstrophy_rotated(X, Y, om, x_le, alpha, DX,
                                                     buffer_dx=BUFFER_DX, half_width=0.5)
            rows.append(dict(alpha=alpha, impl=impl, enstrophy=m["enstrophy"],
                              enstrophy_rotated=ens_rot, int_abs=m["int_abs"], peak=m["peak"]))

    with open(c.DATA / "test7_alpha_sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["alpha", "impl", "enstrophy", "enstrophy_rotated", "int_abs", "peak"])
        w.writeheader(); w.writerows(rows)
    print("wrote test7_alpha_sweep.csv")

    py_rows_ = [r for r in rows if r["impl"] == "py"]
    cpp_rows_ = [r for r in rows if r["impl"] == "cpp"]
    rot_py = np.array([r["enstrophy_rotated"] for r in sorted(py_rows_, key=lambda r: r["alpha"])])
    rot_cpp = np.array([r["enstrophy_rotated"] for r in sorted(cpp_rows_, key=lambda r: r["alpha"])])
    rot_reldiff = np.abs(rot_py - rot_cpp) / np.maximum(rot_py, 1e-12)
    print(f"py vs cpp (ROTATED frame): max reldiff = {rot_reldiff.max():.2e}, mean = {rot_reldiff.mean():.2e}")

    py_rows = sorted([r for r in rows if r["impl"] == "py"], key=lambda r: r["alpha"])
    cpp_rows = sorted([r for r in rows if r["impl"] == "cpp"], key=lambda r: r["alpha"])
    py_ens = np.array([r["enstrophy"] for r in py_rows])
    cpp_ens = np.array([r["enstrophy"] for r in cpp_rows])
    reldiff = np.abs(py_ens - cpp_ens) / np.maximum(py_ens, 1e-12)
    print(f"py vs cpp: max reldiff over sweep = {reldiff.max():.2e}, mean = {reldiff.mean():.2e}")

    rot_py_ens = np.array([r["enstrophy_rotated"] for r in py_rows])

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    ax = axes[0]
    ax.plot([r["alpha"] for r in py_rows], py_ens, "o-", color="#1f77b4",
             label="body-frame (fixed x<x_LE)", lw=1.8, ms=5)
    ax.plot([r["alpha"] for r in py_rows], rot_py_ens, "s-", color="#c0392b",
             label="lab-frame (rotated, tracks flow)", lw=1.8, ms=5)
    ax.set_yscale("log")
    ax.set_xlabel("alpha (degrees)"); ax.set_ylabel("upstream enstrophy")
    ax.set_title("Upstream enstrophy vs angle of attack:\nbody-frame vs. lab-frame-rotated region", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    ax = axes[1]
    py_peak = [r["peak"] for r in py_rows]
    ax.plot([r["alpha"] for r in py_rows], py_peak, "o-", color="#8e44ad", lw=1.6, ms=4)
    ax.set_xlabel("alpha (degrees)"); ax.set_ylabel("upstream peak |omega| (body-frame)")
    ax.set_title("Upstream peak vs angle of attack", fontsize=10)
    ax.grid(alpha=0.3)

    # third panel: visual explainer -- show both regions overlaid on one
    # example field (a middling alpha where the rotation is clearly visible)
    ax = axes[2]
    example_alpha = 40
    run_dir = run_dir_for(example_alpha, "py")
    om_ex = c.load_omega(run_dir, NSTEPS)
    body_mask = np.zeros_like(om_ex, dtype=bool)
    body_mask[c.upstream_mask(X, x_le, BUFFER_DX, DX), :] = \
        (Y[0, :] >= YLIM[0]) & (Y[0, :] <= YLIM[1])
    rot_mask = c.upstream_mask_2d_rotated(X, Y, x_le, example_alpha, BUFFER_DX, DX, half_width=0.5)
    Xw, Yw, om_w, ix, iy = c.window(X, Y, om_ex, (-2.0, 0.5), (-1.0, 1.0))
    V = 3.0
    ax.pcolormesh(Xw, Yw, np.clip(om_w[:-1, :-1], -V, V), shading="flat", cmap="Greys", vmin=-V, vmax=V)
    body_w = body_mask[np.ix_(ix, iy)]
    rot_w = rot_mask[np.ix_(ix, iy)]
    overlay = np.zeros((*body_w.shape, 4))
    overlay[body_w] = (0.12, 0.47, 0.71, 0.45)   # blue = body-frame region
    overlay[rot_w] = (0.75, 0.22, 0.17, 0.45)    # red = lab-frame region
    both = body_w & rot_w
    overlay[both] = (0.5, 0.35, 0.44, 0.55)      # blended = overlap
    ax.imshow(np.transpose(overlay, (1, 0, 2)), extent=(Xw.min(), Xw.max(), Yw.min(), Yw.max()),
              origin="lower", aspect="auto", zorder=2)
    ax.set_xlim(-2.0, 0.5); ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(f"What the two regions look like\n(example: alpha={example_alpha}deg; "
                 "blue=body-frame, red=lab-frame)", fontsize=10)

    fig.suptitle("Test 7: angle-of-attack sweep -- force-magnitude-driven or purely geometric?\n"
                 f"NACA0012, Re=1000, steady, dx=0.02, alpha=0..40,50,60 "
                 f"(py/cpp max reldiff: body-frame={reldiff.max():.1e}, lab-frame={rot_reldiff.max():.1e})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.86])
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
