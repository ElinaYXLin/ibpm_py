"""
test2_grid_refinement.py

Group 2 of the LE/TE striping investigation (see README.md): "is the grid
too coarse?" NACA0012, alpha=0deg, steady, Re=1000, at dx=0.01 and
dx=0.005. Runs BOTH py_static and cpp_static at every new dx (per mentor
request, so every plot below shows both implementations, not just one) --
even though Test 1a already showed the two agree to ~1e-13 relative at
dx=0.02, this repeats that check at each new resolution rather than
assuming it holds. dx=0.02 reuses the existing ../1-paper_based alpha=0
runs (both implementations, already computed there).

Domain, dt/dx=0.5 rule, and run duration (t=30) all match the conventions
already established in ../1-paper_based (dx=0.02) and
../../2-leading_edge_investigation/run_grid_refinement.py (dx=0.01/0.005,
same rule, at Re=500).

Usage:
  python3 test2_grid_refinement.py run     # launch missing runs (blocking)
  python3 test2_grid_refinement.py analyze # analyze whatever runs exist
Output: runs/grid_refine/dx<dx>/, figures/test2a_*.png, data/test2a_*.csv
"""
import csv
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

ALPHA = 0.0
T_FINAL = 30.0
LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
TE_XLIM, TE_YLIM = (0.65, 1.35), (-0.25, 0.25)

# dx=0.02 reuses ../1-paper_based's existing alpha=0 py_static run (no rerun)
GRIDS = {
    0.01: dict(dt=0.005, nsteps=6000),
    0.005: dict(dt=0.0025, nsteps=12000),
}


def ensure_geom(dx):
    """dx=0.01 already has a shared geom in SURF_test/geom/ (used by other
    kurt_comp/vortall work); dx=0.005 does not exist yet anywhere in the
    repo, so it's generated here (and kept local to this folder, since a
    dx=0.005 NACA0012 geom isn't part of any other established sweep)."""
    shared = c.REPO / "SURF_test" / "geom" / f"naca0012_dx{dx:.4f}.geom"
    if shared.exists():
        return shared
    local_geom = c.GEOMDIR / f"naca0012_dx{dx:.4f}.geom"
    local_txt = c.GEOMDIR / f"naca0012_dx{dx:.4f}.txt"
    if not local_geom.exists():
        sys.path.insert(0, str(c.REPO / "SURF_test"))
        from make_airfoil_raw import make_raw_for_dx  # noqa: E402
        n, perim = make_raw_for_dx(str(c.NACA0012_DAT), dx, str(local_txt))
        c.write_geom(local_geom, local_txt)
        print(f"  generated {local_geom.name}: {n} points, perimeter={perim:.4f}")
    return local_geom


def outdir_for(dx, impl):
    # py kept its original (pre-mentor-request) unsuffixed dir name so an
    # already-running/completed py job isn't disturbed; cpp is new, so it
    # gets an explicit suffix
    return c.RUNS / "grid_refine" / (f"dx{dx:.4f}" if impl == "py" else f"dx{dx:.4f}_cpp")


def do_run(impls=("py", "cpp")):
    for dx, params in GRIDS.items():
        geom = ensure_geom(dx)
        nx = int(round(6.0 / dx))
        ny = int(round(3.0 / dx))
        for impl in impls:
            outdir = outdir_for(dx, impl)
            if c.is_done(outdir, params["nsteps"]):
                print(f"dx={dx} {impl}: already done ({outdir})")
                continue
            print(f"dx={dx} {impl}: nx={nx} ny={ny} dt={params['dt']} nsteps={params['nsteps']} -> {outdir}", flush=True)
            ok, elapsed = c.run_case(impl, geom, outdir, nx, ny, params["dt"], params["nsteps"],
                                      alpha=ALPHA, restart=params["nsteps"] // 2)
            print(f"  {'OK' if ok else 'FAILED'} in {elapsed:.0f}s", flush=True)


IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}


def peak_amp(X, Y, om, xlim, ylim):
    ys = Y[0, :]
    iy0 = c.nearest_index(ys, 0.0)
    xs = X[:, 0]
    m = (xs >= xlim[0]) & (xs <= xlim[1])
    lineout = om[:, iy0][m]
    return float(np.abs(lineout).max()), float(lineout.max() - lineout.min())


def analyze():
    # dx -> {impl: (X, Y, om)}
    fields = {0.02: {}}
    fields[0.02]["py"] = (*c.grid_xy(0.02)[:2], c.load_omega(c.KURT1 / "runs" / "dx0.020" / "steady_py_a00", 3000))
    fields[0.02]["cpp"] = (*c.grid_xy(0.02)[:2], c.load_omega(c.KURT1 / "runs" / "dx0.020" / "steady_cpp_a00", 3000))
    for dx, params in GRIDS.items():
        fields[dx] = {}
        for impl in ("py", "cpp"):
            outdir = outdir_for(dx, impl)
            if not c.is_done(outdir, params["nsteps"]):
                print(f"dx={dx} {impl}: run not finished yet, skipping in analysis")
                continue
            Xd, Yd, _ = c.grid_xy(dx)
            fields[dx][impl] = (Xd, Yd, c.load_omega(outdir, params["nsteps"]))

    dxs_present = [dx for dx in fields if fields[dx]]
    rows = []
    for dx in dxs_present:
        for impl, (X, Y, om) in fields[dx].items():
            for xlim, ylim, name in [(LE_XLIM, LE_YLIM, "LE"), (TE_XLIM, TE_YLIM, "TE")]:
                peak, amp = peak_amp(X, Y, om, xlim, ylim)
                rows.append(dict(dx=dx, impl=impl, region=name, peak_abs_omega=peak, amplitude=amp))
    with open(c.DATA / "test2a_grid_refinement.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dx", "impl", "region", "peak_abs_omega", "amplitude"])
        w.writeheader(); w.writerows(rows)

    # fields figure: both implementations, every dx, LE+TE
    n_dx = len(dxs_present)
    fig, axes = plt.subplots(2, 2 * n_dx, figsize=(4.0 * n_dx, 8), squeeze=False)
    for col_dx, dx in enumerate(sorted(dxs_present, reverse=True)):
        for col_impl, impl in enumerate(("py", "cpp")):
            if impl not in fields[dx]:
                continue
            X, Y, om = fields[dx][impl]
            col = col_dx * 2 + col_impl
            for row, (xlim, ylim, name) in enumerate([(LE_XLIM, LE_YLIM, "LE"), (TE_XLIM, TE_YLIM, "TE")]):
                peak, _ = peak_amp(X, Y, om, xlim, ylim)
                ax = axes[row, col]
                Xw, Yw, fw, _, _ = c.window(X, Y, om, xlim, ylim)
                V = 8.0
                ax.pcolormesh(Xw, Yw, np.clip(fw[:-1, :-1], -V, V), shading="flat", cmap="jet", vmin=-V, vmax=V)
                ax.set_aspect("equal")
                ax.set_title(f"{name}, dx={dx}, {IMPL_LABEL[impl]}\npeak|omega|={peak:.2f}", fontsize=8)
    fig.suptitle("Test 2a: grid refinement at LE/TE, NACA0012 alpha=0deg steady, Re=1000, "
                 "py_static vs cpp_static (raw cells, flat shading)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = c.FIGS / "test2a_grid_refinement_fields.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    # combined peak-vs-dx plot: py and cpp both in the same panel, one panel per region
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for name, ax in zip(("LE", "TE"), axes):
        for impl in ("py", "cpp"):
            sub = [r for r in rows if r["region"] == name and r["impl"] == impl]
            sub.sort(key=lambda r: -r["dx"])
            if not sub:
                continue
            dxs = [r["dx"] for r in sub]
            peaks = [r["peak_abs_omega"] for r in sub]
            ax.plot(dxs, peaks, "o-", color=IMPL_COLOR[impl], label=IMPL_LABEL[impl], lw=1.6, ms=6)
        ax.set_xscale("log"); ax.invert_xaxis()
        ax.set_xlabel("dx (log scale, refining -->)"); ax.set_ylabel("peak |omega| near " + name)
        ax.set_title(f"{name}: peak vorticity vs grid resolution")
        ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.suptitle("Test 2a: does the LE/TE peak vanish (grid-too-coarse) or persist (not grid-fixable)\n"
                 "as dx -> 0? py_static vs cpp_static, both shown", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "test2a_peak_vs_dx.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    print("\nTest 2a summary:")
    for r in rows:
        print(f"  dx={r['dx']:.4f} {IMPL_LABEL[r['impl']]:<10} {r['region']}: "
              f"peak|omega|={r['peak_abs_omega']:.3f}, amplitude={r['amplitude']:.3f}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        do_run()
    else:
        analyze()
