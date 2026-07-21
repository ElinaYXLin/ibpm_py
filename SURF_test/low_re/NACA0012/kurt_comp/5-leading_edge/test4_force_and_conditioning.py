"""
test4_force_and_conditioning.py

Group 4 of the LE/TE striping investigation (see README.md): localizing
the mechanism inside the algorithm, near-free (no new solver runs -- reuses
../1-paper_based's existing dx=0.02 steady snapshots, and the baseline
NACA0012 .geom file directly).

Test 4a -- boundary constraint force f (the Lagrange multiplier the
solver's projection step outputs, stored in every State/.bin restart file)
vs. arc length s, zoomed at the LE and TE, with BOTH py_static and
cpp_static overlaid in the same panel. A spurious sawtooth / high
point-to-point oscillation in f concentrated there would localize the
origin to the constraint (projection) solve itself, rather than the flow
field.

Test 4b -- purely geometric, no solver output involved: boundary
point-to-point spacing ds(s) and local radius of curvature 1/kappa(s)
around the perimeter, from the .geom file alone. Where 1/kappa(s)
approaches or drops below ds(s), the boundary curves faster than the
Lagrangian points can track it -- exactly the "clustering relative to
curvature" conditioning risk this codebase's cholesky_solver.py flags
(a non-SPD/near-singular projection matrix from an "over-resolved
boundary" -- here read as "under-resolved relative to curvature").

Usage: python3 test4_force_and_conditioning.py
Output: figures/test4a_*.png, figures/test4b_*.png, data/test4*.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

ALPHAS = [0, 9, 12]
STEP = 3000
S_WINDOW = 0.15  # arc-length half-window around LE/TE for the zoomed plots


def wrapped_dist(s, s0, perimeter):
    d = np.abs(s - s0)
    return np.minimum(d, perimeter - d)


IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}


def test4a():
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    s, perim, s_le, s_te = g["s"], g["perimeter"], g["s_le"], g["s_te"]
    n = g["n"]

    rows = []
    fig, axes = plt.subplots(len(ALPHAS), 2, figsize=(11, 3.2 * len(ALPHAS)), squeeze=False)
    for row, alpha in enumerate(ALPHAS):
        fmag_by_impl = {}
        for impl in ("py", "cpp"):
            run_dir = c.KURT1 / "runs" / "dx0.020" / f"steady_{impl}_a{alpha:02d}"
            st = c.load_state(run_dir, STEP)
            fx = st.f._data[0 * n:1 * n]
            fy = st.f._data[1 * n:2 * n]
            fmag_by_impl[impl] = np.hypot(fx, fy)

        for col, (s0, label) in enumerate([(s_le, "LE"), (s_te, "TE")]):
            d = wrapped_dist(s, s0, perim)
            m = d < S_WINDOW
            idx = np.array(sorted(np.where(m)[0], key=lambda i: (s[i] - s0 + perim / 2) % perim))
            s_rel = np.array([(s[i] - s0 + perim / 2) % perim - perim / 2 for i in idx])

            ax = axes[row, col]
            for impl in ("py", "cpp"):
                f_here = fmag_by_impl[impl][idx]
                d2 = np.diff(f_here, 2) if len(f_here) > 2 else np.array([np.nan])
                sawtooth = float(np.abs(d2).max()) if len(d2) and not np.all(np.isnan(d2)) else float("nan")
                rows.append(dict(alpha=alpha, impl=impl, region=label, max_f=float(f_here.max()),
                                  sawtooth_2nd_diff=sawtooth))
                ax.plot(s_rel, f_here, "o-" if impl == "py" else "s--", ms=4, lw=1.3,
                        color=IMPL_COLOR[impl], label=IMPL_LABEL[impl], alpha=0.85)
            ax.axvline(0, color="0.6", lw=0.6, ls=":")
            ax.set_title(f"alpha={alpha}deg, {label}: |f| vs arc length, py vs cpp", fontsize=9)
            ax.set_xlabel("s - s_" + label.lower()); ax.legend(fontsize=7)

    fig.suptitle("Test 4a: boundary constraint force |f| vs arc length near LE/TE\n"
                 "(sawtoothing here would localize the artifact to the projection solve)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = c.FIGS / "test4a_force_vs_arclength.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    with open(c.DATA / "test4a_force_vs_arclength.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["alpha", "impl", "region", "max_f", "sawtooth_2nd_diff"])
        w.writeheader(); w.writerows(rows)
    print("\nTest 4a summary:")
    for r in rows:
        print(f"  alpha={r['alpha']:>2} {IMPL_LABEL[r['impl']]:<10} {r['region']}: max|f|={r['max_f']:.3f}, "
              f"sawtooth (max |2nd diff|)={r['sawtooth_2nd_diff']:.3g}")


def curvature(x, y):
    """Discrete curvature at each closed-polyline vertex via the turning
    angle between adjacent segments, divided by the local segment length
    -- kappa ~ dtheta/ds."""
    xp = np.roll(x, -1); yp = np.roll(y, -1)
    xm = np.roll(x, 1); ym = np.roll(y, 1)
    v1 = np.stack([x - xm, y - ym], axis=1)
    v2 = np.stack([xp - x, yp - y], axis=1)
    l1 = np.hypot(v1[:, 0], v1[:, 1])
    l2 = np.hypot(v2[:, 0], v2[:, 1])
    cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
    dot = v1[:, 0] * v2[:, 0] + v1[:, 1] * v2[:, 1]
    dtheta = np.arctan2(cross, dot)
    ds_local = 0.5 * (l1 + l2)
    return np.abs(dtheta) / np.maximum(ds_local, 1e-12), ds_local


def test4b():
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x, y, s, perim = g["x"], g["y"], g["s"], g["perimeter"]
    kappa, ds_local = curvature(x, y)
    r_curv = 1.0 / np.maximum(kappa, 1e-9)

    order = np.argsort(s)
    s_sorted = s[order]; ds_sorted = ds_local[order]; r_sorted = r_curv[order]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    ax.plot(s_sorted, ds_sorted, "-o", ms=3, color="#2980b9", label="point spacing ds(s)")
    ax.axhline(0.02, color="0.5", lw=0.8, ls="--", label="background grid dx=0.02")
    ax.axvline(g["s_le"], color="#c0392b", lw=0.8, ls=":", label="LE")
    ax.axvline(g["s_te"] if g["s_te"] < perim else g["s_te"] - perim, color="#16a085", lw=0.8, ls=":", label="TE")
    ax.set_xlabel("arc length s"); ax.set_ylabel("ds"); ax.set_title("boundary point spacing around the body")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.semilogy(s_sorted, r_sorted, "-o", ms=3, color="#8e44ad", label="local radius of curvature 1/kappa")
    ax.axhline(0.02, color="0.5", lw=0.8, ls="--", label="ds (background dx=0.02)")
    ax.axvline(g["s_le"], color="#c0392b", lw=0.8, ls=":", label="LE")
    ax.axvline(g["s_te"] if g["s_te"] < perim else g["s_te"] - perim, color="#16a085", lw=0.8, ls=":", label="TE")
    ax.set_xlabel("arc length s"); ax.set_ylabel("radius of curvature (log scale)")
    ax.set_title("where does 1/kappa approach or drop below ds?")
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")

    fig.suptitle("Test 4b: boundary-point spacing vs local curvature scale (geometry only, no solve)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = c.FIGS / "test4b_spacing_vs_curvature.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    rows = [dict(s=float(s_sorted[i]), ds=float(ds_sorted[i]), radius_of_curvature=float(r_sorted[i]))
            for i in range(len(s_sorted))]
    with open(c.DATA / "test4b_spacing_vs_curvature.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["s", "ds", "radius_of_curvature"])
        w.writeheader(); w.writerows(rows)

    i_le, i_te = g["i_le"], g["i_te"]
    print("\nTest 4b summary:")
    print(f"  at LE: ds={ds_local[i_le]:.4f}, radius_of_curvature={r_curv[i_le]:.4f} "
          f"(r_curv/ds = {r_curv[i_le] / ds_local[i_le]:.2f})")
    print(f"  at TE: ds={ds_local[i_te]:.4f}, radius_of_curvature={r_curv[i_te]:.4f} "
          f"(r_curv/ds = {r_curv[i_te] / ds_local[i_te]:.2f})")
    print(f"  min radius_of_curvature anywhere on body: {r_curv.min():.4f} "
          f"(at s={s[np.argmin(r_curv)]:.4f})")


if __name__ == "__main__":
    test4a()
    test4b()
