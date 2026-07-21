"""
test3_shape_and_spacing.py

Group 3 of the LE/TE striping investigation (see README.md): "is it the
boundary discretization, or the body shape?" All cases: Re=1000, alpha=0,
steady, dx=0.02 (background grid unchanged throughout). Runs BOTH
py_static and cpp_static for every case (per mentor request) so every
plot below shows both implementations. Geometry variants built by
make_geoms.py (run that first).

Test 3a -- boundary-point spacing at fixed dx: LTEdense (ds->dx/4 at LE+TE)
and LTEsparse (ds->4dx at LE+TE), vs. the NACA0012 baseline (ds=dx).
Extends ../../2-leading_edge_investigation's LE-only/denser-only result
(worse at Re=500) to Re=1000, both ends, and adds a sparser arm.

Test 3b -- curvature/bluntness sweep: NACA0006 (sharper) / NACA0012
(baseline) / NACA0018 (blunter) / cylinder (very blunt, constant
curvature), plus a NACA0012-roundTE variant isolating the TE specifically.

Usage:
  python3 test3_shape_and_spacing.py run
  python3 test3_shape_and_spacing.py analyze
Output: runs/shape_spacing/<case>/, figures/test3*.png, data/test3*.csv
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

ALPHA = 0.0
DT, NSTEPS = 0.01, 3000  # matches ../1-paper_based's dx=0.02 steady convention (t=30)
HALF_W_X, HALF_W_Y = 0.35, 0.25  # window half-widths around each shape's own LE/TE

CASES_3A = {
    "naca0012_baseline": c.BASE_GEOM_DX002,
    "naca0012_LTEdense": c.GEOMDIR / "naca0012_dx0.0200_LTEdense.geom",
    "naca0012_LTEsparse": c.GEOMDIR / "naca0012_dx0.0200_LTEsparse.geom",
}
CASES_3B = {
    "naca0006": c.GEOMDIR / "naca0006_dx0.0200.geom",
    "naca0012_baseline": c.BASE_GEOM_DX002,
    "naca0018": c.GEOMDIR / "naca0018_dx0.0200.geom",
    "naca0012_roundTE": c.GEOMDIR / "naca0012_dx0.0200_roundTE.geom",
    "cylinder": c.REPO / "SURF_test" / "vortall" / "3-grid_refine" / "geom" / "cylinder_dx0.0200.geom",
}
ALL_CASES = {**CASES_3A, **CASES_3B}
NX, NY = 300, 150  # dx=0.02, 6c x 3c


IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}


def outdir_for(name, impl):
    # py kept its original (pre-mentor-request) unsuffixed dir name so
    # already-completed py runs aren't disturbed; cpp is new, so it gets
    # an explicit suffix
    return c.RUNS / "shape_spacing" / (name if impl == "py" else f"{name}_cpp")


def do_run(impls=("py", "cpp")):
    for name, geom in ALL_CASES.items():
        for impl in impls:
            outdir = outdir_for(name, impl)
            if c.is_done(outdir, NSTEPS):
                print(f"{name} {impl}: already done")
                continue
            print(f"{name} {impl}: running ({geom.name}) -> {outdir}", flush=True)
            ok, elapsed = c.run_case(impl, geom, outdir, NX, NY, DT, NSTEPS, alpha=ALPHA,
                                      restart=NSTEPS // 2)
            print(f"  {'OK' if ok else 'FAILED'} in {elapsed:.0f}s", flush=True)


def le_te_windows(geom_path):
    g = c.load_geom_points(geom_path)
    x_le, y_le = g["x"][g["i_le"]], g["y"][g["i_le"]]
    x_te, y_te = g["x"][g["i_te"]], g["y"][g["i_te"]]
    le_win = ((x_le - HALF_W_X, x_le + HALF_W_X), (y_le - HALF_W_Y, y_le + HALF_W_Y))
    te_win = ((x_te - HALF_W_X, x_te + HALF_W_X), (y_te - HALF_W_Y, y_te + HALF_W_Y))
    return le_win, te_win, g


def peak_and_amp(X, Y, om, xlim, ylim, ys):
    # lineout at y=0: alpha=0 here, and every shape (airfoils + cylinder) is
    # symmetric about y=0, so this always cuts straight through the LE/TE
    iy0 = c.nearest_index(ys, 0.0)
    xs = X[:, 0]
    m = (xs >= xlim[0]) & (xs <= xlim[1])
    lineout = om[:, iy0][m]
    return float(np.abs(lineout).max()), float(lineout.max() - lineout.min())


def analyze(cases, tag, title):
    X, Y, _ = c.grid_xy(0.02)
    ys = Y[0, :]
    rows = []
    impls_present = set()
    fields_cache = {}
    for name, geom in cases.items():
        for impl in ("py", "cpp"):
            outdir = outdir_for(name, impl)
            if not c.is_done(outdir, NSTEPS):
                print(f"{name} {impl}: run not finished, skipping")
                continue
            om = c.load_omega(outdir, NSTEPS)
            fields_cache[(name, impl)] = om
            impls_present.add(impl)
            le_win, te_win, g = le_te_windows(geom)
            for win, label in [(le_win, "LE"), (te_win, "TE")]:
                xlim, ylim = win
                peak, amp = peak_and_amp(X, Y, om, xlim, ylim, ys)
                rows.append(dict(case=name, impl=impl, region=label, peak_abs_omega=peak,
                                  amplitude=amp, n_points=g["n"], perimeter=g["perimeter"]))

    impls_present = [i for i in ("py", "cpp") if i in impls_present]

    # fields figure: one column per (case, impl), 2 rows (LE, TE)
    n_cols = len(cases) * len(impls_present)
    fig, axes = plt.subplots(2, n_cols, figsize=(3.4 * n_cols, 8), squeeze=False)
    col = 0
    for name, geom in cases.items():
        le_win, te_win, g = le_te_windows(geom)
        for impl in impls_present:
            if (name, impl) not in fields_cache:
                col += 1
                continue
            om = fields_cache[(name, impl)]
            for row, (win, label) in enumerate([(le_win, "LE"), (te_win, "TE")]):
                xlim, ylim = win
                peak, _ = peak_and_amp(X, Y, om, xlim, ylim, ys)
                ax = axes[row, col]
                Xw, Yw, fw, _, _ = c.window(X, Y, om, xlim, ylim)
                V = 8.0
                ax.pcolormesh(Xw, Yw, np.clip(fw[:-1, :-1], -V, V), shading="flat", cmap="jet", vmin=-V, vmax=V)
                ax.set_aspect("equal")
                ax.set_title(f"{name}, {IMPL_LABEL[impl]}\n{label}: peak|omega|={peak:.2f}", fontsize=7.5)
            col += 1
    fig.suptitle(title + " -- py_static vs cpp_static", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = c.FIGS / f"test3_{tag}_fields.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    with open(c.DATA / f"test3_{tag}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "impl", "region", "peak_abs_omega", "amplitude", "n_points", "perimeter"])
        w.writeheader()
        w.writerows(rows)

    # grouped bar chart: LE/TE panels, py+cpp bars side by side per case
    fig, axes = plt.subplots(1, 2, figsize=(6.5 * len(cases) / 3 + 3, 4.8))
    names = list(cases.keys())
    xpos = np.arange(len(names))
    width = 0.8 / max(len(impls_present), 1)
    for ax, region in zip(axes, ("LE", "TE")):
        for k, impl in enumerate(impls_present):
            vals = [next((r["peak_abs_omega"] for r in rows
                          if r["case"] == n and r["region"] == region and r["impl"] == impl), np.nan)
                     for n in names]
            offset = (k - (len(impls_present) - 1) / 2) * width
            ax.bar(xpos + offset, vals, width=width, label=IMPL_LABEL[impl], color=IMPL_COLOR[impl])
        ax.set_xticks(xpos); ax.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel(f"{region} peak |omega|"); ax.legend(fontsize=8); ax.grid(alpha=0.3, axis="y")
        ax.set_title(f"{region}", fontsize=10)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    out = c.FIGS / f"test3_{tag}_bar.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    print(f"\n{tag} summary:")
    for r in rows:
        print(f"  {r['case']:<22} {IMPL_LABEL[r['impl']]:<10} {r['region']}: "
              f"peak|omega|={r['peak_abs_omega']:.3f}, amplitude={r['amplitude']:.3f} (n={r['n_points']})")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        do_run()
    else:
        analyze(CASES_3A, "3a_spacing", "Test 3a: boundary-point spacing at fixed dx=0.02 "
                 "(NACA0012, alpha=0, Re=1000, steady)")
        analyze(CASES_3B, "3b_shape", "Test 3b: curvature/bluntness sweep at fixed dx=0.02 "
                 "(alpha=0, Re=1000, steady)")
