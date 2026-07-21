"""
test0_data_vs_render.py

Group 0 of the leading/trailing-edge (LE/TE) striping investigation (see
README.md). Uses ONLY the already-existing steady, dx=0.02, Re=1000
snapshots in ../1-paper_based/runs/dx0.020/steady_{py,cpp}_a{00,09,12}/ --
zero new solver runs.

Test 0a -- is the striping in the raw omega array, or only in contourf's
rendering? For each alpha, plot the same LE/TE window three ways: (1) the
usual smooth contourf render (py_static), (2) pcolormesh with flat
(no-interpolation) shading, one rectangle per grid cell (py_static), and
(3) a raw 1-D lineout of omega vs x through the LE/TE at the grid row
nearest y=0, with BOTH py_static and cpp_static overlaid in the same
panel, plus the paper's own cropped field for the same alpha/motion for
direct visual reference. If the per-cell values in (2)/(3) themselves
oscillate cell-to-cell, the striping is a real feature of the solved
field, not a contour-interpolation artifact.

Test 0b -- quantify the LE ringing wavelength/amplitude and the TE blob
size, from the same 1-D lineouts, as the baseline every later group
compares against.

All analysis is done in the solver's own (unrotated) grid frame -- this
solver imposes alpha by rotating the free-stream, not the body (see
../1-paper_based/README.md's "Wake vorticity fields" section), so the
airfoil sits at the same (x,y) location in grid coordinates for every
alpha, and the LE/TE artifact itself is a rigid-frame-independent feature
of the solved field. Rotating to the paper's plotting frame would only
relabel axes, not change any of the conclusions below.

Usage: python3 test0_data_vs_render.py
Output: figures/test0a_*.png, data/test0b_wavelength_amplitude.csv
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

ALPHAS = [0, 9, 12]
IMPL = "py"  # contourf/pcolormesh render panels use py_static; lineout panel overlays both
STEP = 3000  # t=30, developed/final snapshot
DX = 0.02
IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}
PAPER_FIGS = c.KURT1 / "paper_figs"

LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
TE_XLIM, TE_YLIM = (0.65, 1.35), (-0.25, 0.25)


def find_extrema(y):
    """Indices of local maxima/minima of a 1-D array (simple 3-point test)."""
    ext = []
    for i in range(1, len(y) - 1):
        if (y[i] - y[i - 1]) * (y[i + 1] - y[i]) < 0:
            ext.append(i)
    return np.array(ext, dtype=int)


def analyze_lineout(x, om, region_xlim, label, results, alpha):
    m = (x >= region_xlim[0]) & (x <= region_xlim[1])
    xx, oo = x[m], om[m]
    ext = find_extrema(oo)
    if len(ext) >= 2:
        wavelengths_cells = np.diff(ext)
        wavelength_cells = float(np.median(wavelengths_cells))
        wavelength_phys = wavelength_cells * DX
    else:
        wavelength_cells = wavelength_phys = float("nan")
    amp = float(oo[ext].max() - oo[ext].min()) if len(ext) >= 2 else float("nan")
    peak_i = int(np.argmax(np.abs(oo)))
    peak_x = float(xx[peak_i])
    peak_val = float(oo[peak_i])
    n_sign_changes = int(np.sum(np.diff(np.sign(oo)) != 0))
    results.append(dict(alpha=alpha, region=label, n_extrema=len(ext),
                         wavelength_cells=wavelength_cells, wavelength_phys=wavelength_phys,
                         amplitude=amp, peak_x=peak_x, peak_val=peak_val,
                         n_sign_changes_in_window=n_sign_changes))
    return xx, oo, ext


def main():
    X, Y, dx = c.grid_xy(DX)
    xs = X[:, 0]
    ys = Y[0, :]
    iy0 = c.nearest_index(ys, 0.0)

    results = []
    for alpha in ALPHAS:
        run_dir = c.KURT1 / "runs" / "dx0.020" / f"steady_{IMPL}_a{alpha:02d}"
        om = c.load_omega(run_dir, STEP)
        om_cpp = c.load_omega(c.KURT1 / "runs" / "dx0.020" / f"steady_cpp_a{alpha:02d}", STEP)
        lineout = om[:, iy0]  # omega(x, y=0) row
        lineout_cpp = om_cpp[:, iy0]

        xx_le, oo_le, ext_le = analyze_lineout(xs, lineout, LE_XLIM, "LE", results, alpha)
        xx_te, oo_te, ext_te = analyze_lineout(xs, lineout, TE_XLIM, "TE", results, alpha)

        paper_png = PAPER_FIGS / f"steady_a{alpha:02d}.png"
        paper_img = plt.imread(paper_png) if paper_png.exists() else None

        fig, axes = plt.subplots(2, 4, figsize=(19, 7))
        for row, (xlim, ylim, xx, oo, ext, name) in enumerate([
                (LE_XLIM, LE_YLIM, xx_le, oo_le, ext_le, "LE"),
                (TE_XLIM, TE_YLIM, xx_te, oo_te, ext_te, "TE")]):
            Xw, Yw, fw, ix, iyw = c.window(X, Y, om, xlim, ylim)
            V = 8.0

            ax = axes[row, 0]
            ax.contourf(Xw, Yw, np.clip(fw, -V, V), levels=41, cmap="jet", extend="both")
            ax.set_title(f"{name}: contourf (smooth render, {IMPL_LABEL[IMPL]})", fontsize=9)
            ax.set_aspect("equal")

            ax = axes[row, 1]
            # flat shading = one flat-colored rectangle per grid cell, no
            # interpolation between cells -- this is what the DATA looks like
            ax.pcolormesh(Xw, Yw, np.clip(fw[:-1, :-1], -V, V), shading="flat",
                           cmap="jet", vmin=-V, vmax=V)
            ax.set_title(f"{name}: pcolormesh (raw cells, {IMPL_LABEL[IMPL]})", fontsize=9)
            ax.set_aspect("equal")

            ax = axes[row, 2]
            m = (xs >= xlim[0]) & (xs <= xlim[1])
            ax.plot(xx, oo, "o-", ms=3, lw=1, color=IMPL_COLOR["py"], label=IMPL_LABEL["py"])
            ax.plot(xs[m], lineout_cpp[m], "s--", ms=3, lw=1, color=IMPL_COLOR["cpp"],
                    label=IMPL_LABEL["cpp"], alpha=0.85)
            ax.plot(xx[ext], oo[ext], "x", ms=7, color="k", label="py local extrema")
            ax.axhline(0, color="0.6", lw=0.6)
            ax.set_title(f"{name}: raw lineout, omega(x, y=0), py vs cpp", fontsize=9)
            ax.set_xlabel("x"); ax.legend(fontsize=7)

            ax = axes[row, 3]
            if paper_img is not None:
                ax.imshow(paper_img)
                ax.set_title(f"{name}: Kurtulus (2019) crop, alpha={alpha}deg\n(whole-field reference, not zoomed)",
                             fontsize=8.5)
            else:
                ax.text(0.5, 0.5, "no paper crop", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)

        fig.suptitle(f"Test 0a: data vs. render, NACA0012 alpha={alpha}deg steady, Re=1000, dx=0.02, t=30\n"
                     f"py_static vs cpp_static vs Kurtulus (2019)", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        out = c.FIGS / f"test0a_data_vs_render_a{alpha:02d}.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"wrote {out.name}")

    # ---------------- 0b: save the quantified wavelength/amplitude table ----------------
    import csv
    out_csv = c.DATA / "test0b_wavelength_amplitude.csv"
    keys = ["alpha", "region", "n_extrema", "wavelength_cells", "wavelength_phys",
            "amplitude", "peak_x", "peak_val", "n_sign_changes_in_window"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f"wrote {out_csv.name}")

    print("\nTest 0b summary (dx=0.02, cell size = {:.3f}):".format(DX))
    for r in results:
        print(f"  alpha={r['alpha']:>2} {r['region']}: "
              f"wavelength={r['wavelength_cells']:.1f} cells ({r['wavelength_phys']:.4f}c), "
              f"amplitude={r['amplitude']:.3f}, peak|omega|={r['peak_val']:.3f} at x={r['peak_x']:.4f}, "
              f"sign changes in window={r['n_sign_changes_in_window']}")


if __name__ == "__main__":
    main()
