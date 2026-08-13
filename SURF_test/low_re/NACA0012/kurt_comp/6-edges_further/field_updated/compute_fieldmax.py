"""
compute_fieldmax.py

Regenerates every LE/TE "peak |omega|" number 5-leading_edge originally
reported with the unreliable y=0 lineout metric, using the 2-D field-max
metric `../` (6-edges_further) established as reliable instead (see
`../README.md`'s "Bottom line up front" and 5-leading_edge/README.md's
"Correction" section). `../` already recomputed several of these
(Reconciliation 1/2, Group A, Group C1/C2) -- this script's job is to
assemble ALL of them into one self-contained, complete set of tables that
mirror 5-leading_edge's original Group 0b/2a/3a/3b structure exactly, one
table per group, LE+TE, filling in the specific cases `../` never
actually computed with field-max:

  - Test 0b (alpha=0,9,12 wavelength/amplitude table): `../` never
    touched this -- it only ever worked at alpha=0. NEW field-max
    numbers computed here for all 3 angles, both regions.
  - Test 3a's naca0012_LTEsparse case (ds=4dx, LE+TE together): `../`
    Reconciliation 2 covered baseline/LE-only-dense/LE+TE-dense, but
    never LTEsparse specifically. NEW here.
  - Test 3b's naca0012_roundTE and cylinder cases: `../` Group A only
    covered naca0006/0012/0018. NEW here (both LE and TE windows).

Everything else (Test 2a grid refinement; Test 3a baseline/LTEdense;
Test 3b naca0006/0012/0018) is RECOMPUTED here too (not just cited) from
the same on-disk snapshots `../`'s scripts already used, so this folder
is a single complete reference rather than half-original/half-pointer.

**Every table below now also loads cpp_static** for the same case (both
implementations exist on disk for every run already used here) and
reports the relative py/cpp difference on the field-max metric, so the
py/cpp agreement this repo establishes everywhere else is confirmed
here too, not assumed.

Table 1-4: ZERO NEW RUNS (every case reuses a flow snapshot that already
exists on disk, in `../../1-paper_based/`, `../../5-leading_edge/`, or
`../../2-leading_edge_investigation/`). Table 5 (ngrid sweep) launches 8
new runs -- see run_ngrid_sweep.py.

Usage: python3 compute_fieldmax.py
Output: data/*.csv, figures/*.png
"""
import csv
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as c  # noqa: E402  (../common.py -- 6-edges_further's shared helpers)

DATA = HERE / "data"
FIGS = HERE / "figures"
DATA.mkdir(exist_ok=True)
FIGS.mkdir(exist_ok=True)

LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
TE_XLIM, TE_YLIM = (0.65, 1.35), (-0.25, 0.25)
WINDOWS = {"LE": (LE_XLIM, LE_YLIM), "TE": (TE_XLIM, TE_YLIM)}

KURT1 = c.KURT / "1-paper_based"
KURT5 = c.KURT5


def field_max(X, Y, om, region):
    xlim, ylim = WINDOWS[region]
    _, _, sub, _, _ = c.window(X, Y, om, xlim, ylim)
    return float(np.abs(sub).max())


def lineout_max(X, Y, om, region):
    xlim, _ = WINDOWS[region]
    ys = Y[0, :]
    iy0 = c.nearest_index(ys, 0.0)
    xs = X[:, 0]
    m = (xs >= xlim[0]) & (xs <= xlim[1])
    return float(np.abs(om[:, iy0][m]).max())


def cpp_path(py_dir):
    """Map a py_static run directory to its cpp_static counterpart.
    Two naming conventions coexist on disk: 1-paper_based's
    "steady_py_aXX" -> "steady_cpp_aXX" (infix), and 5-leading_edge's
    "<name>" -> "<name>_cpp" (suffix)."""
    py_dir = pathlib.Path(py_dir)
    name = py_dir.name
    if "_py_" in name:
        return py_dir.with_name(name.replace("_py_", "_cpp_"))
    return py_dir.with_name(name + "_cpp")


def load_both_impls(py_dir, step, dx=0.02):
    X, Y = c.grid_xy(dx)
    out = {}
    for impl, run_dir in (("py", py_dir), ("cpp", cpp_path(py_dir))):
        om = c.load_omega(run_dir, step)
        out[impl] = om
    return X, Y, out


def rows_for_case(key_fields, X, Y, om_by_impl):
    """One row per (region, impl), plus the py/cpp relative difference on
    field-max computed once per region and duplicated onto both impl rows
    (so every row is self-describing without a second join)."""
    rows = []
    for region in ("LE", "TE"):
        fm = {impl: field_max(X, Y, om, region) for impl, om in om_by_impl.items()}
        lm = {impl: lineout_max(X, Y, om, region) for impl, om in om_by_impl.items()}
        relerr = abs(fm["py"] - fm["cpp"]) / max(abs(fm["py"]), 1e-300)
        for impl in ("py", "cpp"):
            rows.append(dict(**key_fields, region=region, impl=impl,
                              field_max=fm[impl], lineout_max=lm[impl],
                              fieldmax_py_cpp_relerr=relerr))
    return rows


def write_csv(rows, fname, key_names):
    fieldnames = key_names + ["region", "impl", "field_max", "lineout_max", "fieldmax_py_cpp_relerr"]
    with open(DATA / fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)
    print(f"wrote {fname}")


# ================================================================
# Table 1 -- Test 0b companion: alpha=0,9,12, LE+TE, field-max, py+cpp
# ================================================================
def table1_test0b():
    rows = []
    for alpha in (0, 9, 12):
        py_dir = KURT1 / "runs" / "dx0.020" / f"steady_py_a{alpha:02d}"
        X, Y, om_by_impl = load_both_impls(py_dir, 3000)
        rows += rows_for_case(dict(alpha=alpha), X, Y, om_by_impl)
    write_csv(rows, "table1_test0b_angles.csv", ["alpha"])
    return rows


# ================================================================
# Table 2 -- Test 2a companion: grid refinement, LE+TE, field-max, py+cpp
# ================================================================
def table2_test2a():
    cases = [
        (0.02, KURT1 / "runs" / "dx0.020" / "steady_py_a00", 3000),
        (0.01, KURT5 / "runs" / "grid_refine" / "dx0.0100", 6000),
        (0.005, KURT5 / "runs" / "grid_refine" / "dx0.0050", 12000),
    ]
    rows = []
    for dx, py_dir, step in cases:
        X, Y, om_by_impl = load_both_impls(py_dir, step, dx=dx)
        rows += rows_for_case(dict(dx=dx), X, Y, om_by_impl)
    write_csv(rows, "table2_test2a_grid_refinement.csv", ["dx"])
    return rows


# ================================================================
# Table 3 -- Test 3a companion: LTEsparse / baseline / LTEdense,
#            LE+TE, field-max, py+cpp
# ================================================================
def table3_test3a():
    cases = [
        ("naca0012_LTEsparse (ds=4dx)", KURT5 / "runs" / "shape_spacing" / "naca0012_LTEsparse"),
        ("naca0012_baseline (ds=dx)", KURT5 / "runs" / "shape_spacing" / "naca0012_baseline"),
        ("naca0012_LTEdense (ds=dx/4)", KURT5 / "runs" / "shape_spacing" / "naca0012_LTEdense"),
    ]
    rows = []
    for label, py_dir in cases:
        X, Y, om_by_impl = load_both_impls(py_dir, 3000)
        rows += rows_for_case(dict(case=label), X, Y, om_by_impl)
    write_csv(rows, "table3_test3a_spacing.csv", ["case"])
    return rows


# ================================================================
# Table 4 -- Test 3b companion: naca0006/0012/0018/roundTE/cylinder,
#            LE+TE, field-max, py+cpp
# ================================================================
def table4_test3b():
    cases = [
        ("naca0006", KURT5 / "runs" / "shape_spacing" / "naca0006"),
        ("naca0012 (baseline)", KURT5 / "runs" / "shape_spacing" / "naca0012_baseline"),
        ("naca0018", KURT5 / "runs" / "shape_spacing" / "naca0018"),
        ("naca0012_roundTE", KURT5 / "runs" / "shape_spacing" / "naca0012_roundTE"),
        ("cylinder", KURT5 / "runs" / "shape_spacing" / "cylinder"),
    ]
    rows = []
    for label, py_dir in cases:
        X, Y, om_by_impl = load_both_impls(py_dir, 3000)
        rows += rows_for_case(dict(case=label), X, Y, om_by_impl)
    write_csv(rows, "table4_test3b_shape.csv", ["case"])
    return rows


# ================================================================
# Table 5 -- ngrid sweep (NEW runs, see run_ngrid_sweep.py): default
#            settings (dx=0.02, ds=dx, alpha=0, Re=1000) with ngrid
#            1,2,3,4 as the only variable. LE+TE, field-max, py+cpp.
# ================================================================
NGRID_DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.52)  # ny%4==0 fix, same as elsewhere in this repo
NGRID_NX, NGRID_NY = 300, 152
NGRID_DT, NGRID_NSTEPS = 0.01, 3000


def table5_ngrid_sweep():
    rows = []
    for ngrid in (1, 2, 3, 4):
        py_dir = HERE / "runs" / "ngrid_sweep" / f"ngrid{ngrid}_py"
        X, Y = c.grid_xy(0.02, length=NGRID_DOMAIN["length"], xoffset=NGRID_DOMAIN["xoffset"],
                          yoffset=NGRID_DOMAIN["yoffset"], yheight=NGRID_NY * 0.02)
        om_by_impl = {}
        for impl in ("py", "cpp"):
            run_dir = HERE / "runs" / "ngrid_sweep" / f"ngrid{ngrid}_{impl}"
            om_by_impl[impl] = c.load_omega(run_dir, NGRID_NSTEPS)
        rows += rows_for_case(dict(ngrid=ngrid), X, Y, om_by_impl)
    write_csv(rows, "table5_ngrid_sweep.csv", ["ngrid"])
    return rows


def bar_compare(rows, key, fname, title, group_order=None):
    """Grouped bar chart: field-max vs lineout-max (py only, for
    readability), LE and TE side by side."""
    keys = group_order or sorted(set(r[key] for r in rows), key=str)
    fig, axes = plt.subplots(1, 2, figsize=(max(7, 1.6 * len(keys)), 5))
    for ax, region in zip(axes, ("LE", "TE")):
        fm = [next(r["field_max"] for r in rows if r[key] == k and r["region"] == region and r["impl"] == "py")
              for k in keys]
        lm = [next(r["lineout_max"] for r in rows if r[key] == k and r["region"] == region and r["impl"] == "py")
              for k in keys]
        x = np.arange(len(keys))
        ax.bar(x - 0.2, fm, 0.4, color="#c0392b", label="field-max (reliable)")
        ax.bar(x + 0.2, lm, 0.4, color="#7f8c8d", label="lineout-max (unreliable)")
        ax.set_xticks(x); ax.set_xticklabels([str(k) for k in keys], rotation=30, ha="right", fontsize=8)
        ax.set_title(f"{region} window", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        if ax is axes[0]:
            ax.set_ylabel("peak |omega|")
            ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(FIGS / fname, dpi=140)
    plt.close(fig)
    print(f"wrote {fname}")


def print_pycpp_summary(all_rows, label):
    relerrs = [r["fieldmax_py_cpp_relerr"] for r in all_rows]
    print(f"  {label}: max py/cpp relative diff (field-max) = {max(relerrs):.3e}, "
          f"n={len(relerrs)//2} cases")


def main():
    r1 = table1_test0b()
    r2 = table2_test2a()
    r3 = table3_test3a()
    r4 = table4_test3b()

    bar_compare(r1, "alpha", "table1_test0b_angles.png",
                "Test 0b companion: field-max vs lineout, alpha=0/9/12",
                group_order=[0, 9, 12])
    bar_compare(r2, "dx", "table2_test2a_grid_refinement.png",
                "Test 2a companion: field-max vs lineout, grid refinement",
                group_order=[0.02, 0.01, 0.005])
    bar_compare(r3, "case", "table3_test3a_spacing.png",
                "Test 3a companion: field-max vs lineout, boundary-point spacing",
                group_order=["naca0012_LTEsparse (ds=4dx)", "naca0012_baseline (ds=dx)",
                             "naca0012_LTEdense (ds=dx/4)"])
    bar_compare(r4, "case", "table4_test3b_shape.png",
                "Test 3b companion: field-max vs lineout, shape/thickness family",
                group_order=["naca0006", "naca0012 (baseline)", "naca0018",
                             "naca0012_roundTE", "cylinder"])

    print("\npy/cpp agreement summary (field-max metric):")
    print_pycpp_summary(r1, "Table 1 (angles)")
    print_pycpp_summary(r2, "Table 2 (grid refinement)")
    print_pycpp_summary(r3, "Table 3 (spacing)")
    print_pycpp_summary(r4, "Table 4 (shape family)")

    ngrid_runs_done = all((HERE / "runs" / "ngrid_sweep" / f"ngrid{n}_{impl}" / f"flow{NGRID_NSTEPS:05d}.bin").exists()
                           for n in (1, 2, 3, 4) for impl in ("py", "cpp"))
    if ngrid_runs_done:
        r5 = table5_ngrid_sweep()
        bar_compare(r5, "ngrid", "table5_ngrid_sweep.png",
                    "Table 5: field-max vs lineout, ngrid=1,2,3,4 (dx=0.02, ds=dx, alpha=0 fixed)",
                    group_order=[1, 2, 3, 4])
        print_pycpp_summary(r5, "Table 5 (ngrid sweep)")
    else:
        print("\nTable 5 (ngrid sweep) skipped -- run run_ngrid_sweep.py first.")


if __name__ == "__main__":
    main()
