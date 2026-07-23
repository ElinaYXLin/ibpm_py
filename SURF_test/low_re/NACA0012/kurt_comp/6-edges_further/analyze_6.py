"""
analyze_6.py

Analysis for recon2 and Groups B-E (Group A: testA_metric_and_geom_audit.py;
Group F: testF_conditioning.py; recon1: recon1_grid_refinement.py -- all
already run). Every LE/TE metric here is computed TWO ways -- the y=0
lineout (matching the original Test 3b numbers) and the full 2-D window
max (per Group A's finding that the lineout metric is the fragile one) --
so every table below can be read either way.

Usage: python3 analyze_6.py [recon2|B|C|D|E|all]
Output: data/*.csv (per group), printed summaries
"""
import csv
import sys

import numpy as np

import common as c

DX = 0.02
LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
TE_XLIM, TE_YLIM = (0.65, 1.35), (-0.25, 0.25)
NSTEPS_002 = 3000


def outdir(group, name, impl):
    return c.RUNS / group / (name if impl == "py" else f"{name}_cpp")


def metrics_for(run_dir, geom_path, dx=DX, nsteps=NSTEPS_002, region="LE"):
    if not (run_dir / f"flow{nsteps:05d}.bin").exists():
        return None
    om = c.load_omega(run_dir, nsteps)
    g = c.load_geom_points(geom_path)
    if region == "LE":
        x0, y0 = g["x"][g["i_le"]], g["y"][g["i_le"]]
    else:
        x0, y0 = g["x"][g["i_te"]], g["y"][g["i_te"]]
    xlim = (x0 - 0.35, x0 + 0.35)
    ylim = (y0 - 0.25, y0 + 0.25)
    X, Y = c.grid_xy(dx)
    Xw, Yw, sub, _, _ = c.window(X, Y, om, xlim, ylim)
    field_max = float(np.abs(sub).max())
    ys = Y[0, :]
    iy0 = c.nearest_index(ys, 0.0)
    xs = X[:, 0]
    m = (xs >= xlim[0]) & (xs <= xlim[1])
    lineout_max = float(np.abs(om[:, iy0][m]).max())
    return dict(field_max=field_max, lineout_max=lineout_max)


def write_csv(name, rows):
    if not rows:
        print(f"{name}: no rows (data not ready)")
        return
    with open(c.DATA / f"{name}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote {name}.csv ({len(rows)} rows)")


# ==================================================================== recon2
def analyze_recon2():
    cases = [
        ("LEonly_dense_Re1000", c.GEOMDIR / "naca0012_dx0.0200_LEonly_dense.geom", "recon2", 1000, "LE-only dense"),
        ("LTEdense_Re500", c.KURT5 / "geom" / "naca0012_dx0.0200_LTEdense.geom", "recon2", 500, "LE+TE dense"),
    ]
    rows = []
    for name, geom, group, re, label in cases:
        for impl in ("py", "cpp"):
            rd = outdir(group, name, impl)
            m = metrics_for(rd, geom)
            if m is None:
                continue
            rows.append(dict(case=label, re=re, impl=impl, **m))
    write_csv("recon2_boundary_density", rows)
    if rows:
        print("\nRecon 2 (LE-only-dense vs LE+TE-dense, crossed with Re) -- LE region:")
        for r in rows:
            print(f"  {r['case']} Re={r['re']} {r['impl']}: field_max={r['field_max']:.2f}, "
                  f"lineout_max={r['lineout_max']:.2f}")
        # context: existing data points
        print("  [existing] LE-only-dense Re=500 (older investigation): field_max=89.9 (per its own report)")
        print("  [existing] LE+TE-dense Re=1000 (5-leading_edge Group3a): "
              "field_max=?, lineout_max=5.802 (see 5-leading_edge/data/test3_3a_spacing.csv)")


# ==================================================================== B1
def analyze_B1():
    shifts = [(0, 0), (0.25, 0), (0.5, 0), (0.75, 0), (0, 0.25), (0, 0.5), (0, 0.75), (0.5, 0.5)]
    rows = []
    for fx, fy in shifts:
        name = "baseline" if (fx, fy) == (0, 0) else f"phase_x{fx:.2f}_y{fy:.2f}"
        for impl in ("py", "cpp"):
            if (fx, fy) == (0, 0):
                rd = c.KURT5 / "runs" / "shape_spacing" / ("naca0012_baseline" if impl == "py" else "naca0012_baseline_cpp")
            else:
                rd = outdir("B1_phase", name, impl)
            m = metrics_for(rd, c.BASE_GEOM_DX002)
            if m is None:
                continue
            rows.append(dict(shift_x=fx, shift_y=fy, impl=impl, **m))
    write_csv("testB1_phase_sweep", rows)
    if rows:
        print("\nTest B1 (phase sweep, NACA0012, fixed shape/grid):")
        for r in rows:
            print(f"  shift=({r['shift_x']},{r['shift_y']}) {r['impl']}: "
                  f"field_max={r['field_max']:.2f}, lineout_max={r['lineout_max']:.2f}")


# ==================================================================== B2
def analyze_B2():
    shapes = {
        "naca0006": c.KURT5 / "geom" / "naca0006_dx0.0200.geom",
        "naca0012": c.BASE_GEOM_DX002,
        "naca0018": c.KURT5 / "geom" / "naca0018_dx0.0200.geom",
    }
    rows = []
    for shape, geom in shapes.items():
        for impl in ("py", "cpp"):
            rd = outdir("B2_phase_equalized", f"{shape}_equalized", impl)
            m = metrics_for(rd, geom)
            if m is None:
                continue
            rows.append(dict(shape=shape, impl=impl, **m))
    write_csv("testB2_phase_equalized", rows)
    if rows:
        print("\nTest B2 (0006/0012/0018, phase-equalized):")
        for r in rows:
            print(f"  {r['shape']} {r['impl']}: field_max={r['field_max']:.2f}, lineout_max={r['lineout_max']:.2f}")


# ==================================================================== C1
def analyze_C1():
    rows = []
    for shape in ("naca0006", "naca0018"):
        for dx, nsteps in [(0.01, 6000), (0.005, 12000)]:
            geom = c.GEOMDIR / f"{shape}_dx{dx:.4f}.geom"
            for impl in ("py", "cpp"):
                rd = outdir("C1_shape_refine", f"{shape}_dx{dx:.4f}", impl)
                m = metrics_for(rd, geom, dx=dx, nsteps=nsteps)
                if m is None:
                    continue
                rows.append(dict(shape=shape, dx=dx, impl=impl, **m))
    write_csv("testC1_shape_refinement", rows)
    if rows:
        print("\nTest C1 (per-shape grid refinement, 0006/0018; 0012's is 5-leading_edge Group 2):")
        for r in rows:
            print(f"  {r['shape']} dx={r['dx']} {r['impl']}: field_max={r['field_max']:.2f}, "
                  f"lineout_max={r['lineout_max']:.2f}")


# ==================================================================== C2
def analyze_C2():
    family = {
        "naca0004": 4 / 12, "naca0006": 0.5, "naca0008": 8 / 12, "naca0010": 10 / 12,
        "naca0012": 1.0, "naca0014": 14 / 12, "naca0016": 16 / 12, "naca0018": 1.5, "naca0020": 20 / 12,
    }
    rows = []
    for name, ratio in family.items():
        r_le = c.R_LE_0012 * ratio ** 2
        for impl in ("py", "cpp"):
            if name == "naca0012":
                rd = c.KURT5 / "runs" / "shape_spacing" / ("naca0012_baseline" if impl == "py" else "naca0012_baseline_cpp")
                geom = c.BASE_GEOM_DX002
            elif name in ("naca0006", "naca0018"):
                rd = c.KURT5 / "runs" / "shape_spacing" / (name if impl == "py" else f"{name}_cpp")
                geom = c.KURT5 / "geom" / f"{name}_dx0.0200.geom"
            else:
                rd = outdir("C2_thickness_family", name, impl)
                geom = c.GEOMDIR / f"{name}_dx0.0200.geom"
            m = metrics_for(rd, geom)
            if m is None:
                continue
            rows.append(dict(shape=name, thickness_ratio=ratio, r_le=r_le, r_le_over_dx=r_le / DX,
                              impl=impl, **m))
    write_csv("testC2_thickness_family", rows)
    if rows:
        print("\nTest C2 (thickness family, r_LE/dx sweep):")
        for r in rows:
            print(f"  {r['shape']} (r_LE/dx={r['r_le_over_dx']:.3f}) {r['impl']}: "
                  f"field_max={r['field_max']:.2f}, lineout_max={r['lineout_max']:.2f}")


# ==================================================================== D1
def analyze_D1():
    factors = [0.5, 1, 2, 4, 8, 16]  # 1=baseline (no density change), 4=recon2's LEonly_dense
    rows = []
    for factor in factors:
        if factor == 1:
            rd_py = c.KURT5 / "runs" / "shape_spacing" / "naca0012_baseline"
            rd_cpp = c.KURT5 / "runs" / "shape_spacing" / "naca0012_baseline_cpp"
            geom = c.BASE_GEOM_DX002
        elif factor == 4:
            rd_py = outdir("recon2", "LEonly_dense_Re1000", "py")
            rd_cpp = outdir("recon2", "LEonly_dense_Re1000", "cpp")
            geom = c.GEOMDIR / "naca0012_dx0.0200_LEonly_dense.geom"
        else:
            rd_py = outdir("D1_point_density", f"LEdensity{factor:g}x", "py")
            rd_cpp = outdir("D1_point_density", f"LEdensity{factor:g}x", "cpp")
            geom = c.GEOMDIR / f"naca0012_dx0.0200_LEdensity{factor:g}x.geom"
        for impl, rd in (("py", rd_py), ("cpp", rd_cpp)):
            m = metrics_for(rd, geom)
            if m is None:
                continue
            rows.append(dict(density_factor=factor, impl=impl, **m))
    write_csv("testD1_point_density", rows)
    if rows:
        print("\nTest D1 (LE point-density levels, NACA0012, LE only):")
        for r in rows:
            print(f"  factor={r['density_factor']} {r['impl']}: field_max={r['field_max']:.2f}, "
                  f"lineout_max={r['lineout_max']:.2f}")


# ==================================================================== E1
def analyze_E1():
    names = ["naca0012_dx0.0200_frontsharp_TEnative", "naca0012_dx0.0200_frontblunt_TEnative",
             "naca0012_dx0.0200_TEsharp_frontnative"]
    rows = []
    for name in names:
        geom = c.GEOMDIR / f"{name}.geom"
        for impl in ("py", "cpp"):
            rd = outdir("E1_decouple", name, impl)
            m_le = metrics_for(rd, geom, region="LE")
            m_te = metrics_for(rd, geom, region="TE")
            if m_le is None:
                continue
            rows.append(dict(case=name, impl=impl, le_field_max=m_le["field_max"], le_lineout_max=m_le["lineout_max"],
                              te_field_max=m_te["field_max"], te_lineout_max=m_te["lineout_max"]))
    write_csv("testE1_decouple", rows)
    if rows:
        print("\nTest E1 (decoupled front-only / TE-only variants):")
        for r in rows:
            print(f"  {r['case']} {r['impl']}: LE field_max={r['le_field_max']:.2f}, "
                  f"TE field_max={r['te_field_max']:.2f}")


# ==================================================================== E2
def analyze_E2():
    shapes = {"naca0006_dx0.0200_roundTE": "naca0006", "naca0018_dx0.0200_roundTE": "naca0018"}
    rows = []
    for name, base in shapes.items():
        geom = c.GEOMDIR / f"{name}.geom"
        for impl in ("py", "cpp"):
            rd = outdir("E2_common_TE", name, impl)
            m_le = metrics_for(rd, geom, region="LE")
            m_te = metrics_for(rd, geom, region="TE")
            if m_le is None:
                continue
            rows.append(dict(base_shape=base, impl=impl, le_field_max=m_le["field_max"],
                              le_lineout_max=m_le["lineout_max"], te_field_max=m_te["field_max"],
                              te_lineout_max=m_te["lineout_max"]))
    write_csv("testE2_common_TE", rows)
    if rows:
        print("\nTest E2 (common-TE resweep, naca0006/0018 with the same rounded TE):")
        for r in rows:
            print(f"  {r['base_shape']}+roundTE {r['impl']}: LE field_max={r['le_field_max']:.2f}, "
                  f"TE field_max={r['te_field_max']:.2f}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    fns = dict(recon2=analyze_recon2, B=lambda: (analyze_B1(), analyze_B2()),
               C=lambda: (analyze_C1(), analyze_C2()), D=analyze_D1,
               E=lambda: (analyze_E1(), analyze_E2()))
    if which == "all":
        for fn in fns.values():
            fn()
    else:
        fns[which]()
