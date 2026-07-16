"""
gen_le_investigation_figs.py

Figure generator for all four tests in the leading-edge (LE) vorticity-
speck investigation (see ../README.md "Leading-edge vorticity
investigation" section). Reads the outputs of:
  - run_grid_refinement.py   (test 1: dx=0.01 vs dx=0.005)
  - run_alpha_sweep.py       (test 2: alpha=0,2,8,10 spatial asymmetry)
  - run_le_densified.py      (test 3: LE-densified boundary points)
  - compute_le_residual.py   (test 4: spatial no-slip residual)

Usage: python3 SURF_test/low_re/NACA0012/leading_edge_investigation/gen_le_investigation_figs.py
Output: SURF_test/low_re/NACA0012/leading_edge_investigation/figures/*.png
"""
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "leading_edge_investigation"
DATA = HERE / "_run_data"
FIGS = HERE / "figures"
FIGS.mkdir(parents=True, exist_ok=True)
DAT_PATH = REPO / "SURF_test" / "low_re" / "NACA0012" / "naca0012.dat.txt"
R_LE = 1.1019 * 0.12 ** 2  # =~ 0.01587


def load_dat_pts(path):
    lines = pathlib.Path(path).read_text().splitlines()
    pts = []
    for l in lines[1:]:
        l = l.strip()
        if l:
            x, y = l.split()
            pts.append((float(x), float(y)))
    return np.array(pts)


AIRFOIL_PTS = load_dat_pts(DAT_PATH)


def field_grid(dx, length=6.0, xoffset=-2.0, yoffset=-1.5, height=3.0):
    nx = int(round(length / dx))
    ny = int(round(height / dx))
    xs = xoffset + np.arange(1, nx) * dx
    ys = yoffset + np.arange(1, ny) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return X, Y


def load_omega(path):
    return State(filename=str(path)).omega._data[0].copy()


def le_window_mask(X, Y, half_width=6 * R_LE):
    """Mask selecting the region within `half_width` of the LE (x=0,y=0)."""
    return (np.abs(X) < half_width) & (np.abs(Y) < half_width)


def peak_near_le(field, X, Y, half_width=6 * R_LE):
    mask = le_window_mask(X, Y, half_width)
    if not mask.any():
        return np.nan, (np.nan, np.nan)
    sub = np.where(mask, np.abs(field), -np.inf)
    idx = np.unravel_index(np.argmax(sub), sub.shape)
    return field[idx], (X[idx], Y[idx])


def draw_field_with_zoom(fig, gs_row, field, X, Y, title, vmax=8.0, zoom_half=6 * R_LE):
    """Draw a full-domain vorticity panel plus a zoomed LE inset, side by side."""
    ax_full = fig.add_subplot(gs_row[0])
    ax_zoom = fig.add_subplot(gs_row[1])

    cf = ax_full.contourf(X, Y, np.clip(field, -vmax, vmax), levels=41, cmap="RdBu_r", extend="both")
    ax_full.fill(AIRFOIL_PTS[:, 0], AIRFOIL_PTS[:, 1], color="0.15", zorder=5)
    ax_full.set_xlim(-2, 4); ax_full.set_ylim(-1.5, 1.5); ax_full.set_aspect("equal")
    ax_full.set_title(title, fontsize=9)
    ax_full.add_patch(plt.Rectangle((-zoom_half, -zoom_half), 2 * zoom_half, 2 * zoom_half,
                                     fill=False, edgecolor="lime", lw=1.2, zorder=6))

    peak_val, (px, py) = peak_near_le(field, X, Y, zoom_half)
    ax_zoom.contourf(X, Y, np.clip(field, -vmax, vmax), levels=41, cmap="RdBu_r", extend="both")
    ax_zoom.fill(AIRFOIL_PTS[:, 0], AIRFOIL_PTS[:, 1], color="0.15", zorder=5)
    ax_zoom.plot([px], [py], "*", color="lime", ms=10, mec="k", mew=0.5, zorder=7)
    ax_zoom.set_xlim(-zoom_half, zoom_half); ax_zoom.set_ylim(-zoom_half, zoom_half)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_title(f"LE zoom (peak $\\omega$={peak_val:+.2f})", fontsize=8)
    return cf, peak_val, (px, py)


# ---------------------------------------------------------------------
# Test 1: grid refinement (dx=0.01 vs dx=0.005), LE peak vs dx
# ---------------------------------------------------------------------
def fig_grid_refinement():
    STEPS_T = [0, 5, 10, 15, 20, 25, 30]
    dxs = [0.01, 0.005]
    # dx=0.02 point reuses test 3's uniform-boundary baseline run (identical
    # Re=500, alpha=0, dx=0.02 setup as ../../run_gridconv.py's dx=0.02
    # point, but WITH snapshots), so the grid-refinement trend can include
    # it without a fourth simulation.
    extra_dirs = {0.02: DATA / "le_uniform_baseline_snap"}
    peaks = {}
    fields_at_peak_time = {}
    for dx in dxs:
        rundir = DATA / f"gridconv_dx{dx}_snap"
        dt = dx / 2
        restart = int(round(5.0 / dt))
        X, Y = field_grid(dx)
        avail = sorted(rundir.glob("flow*.bin"))
        if not avail:
            print(f"grid refinement: no snapshots for dx={dx}, skipping")
            return
        # find, across all snapshots, which one has the largest near-LE peak
        best = None
        for f in avail:
            step = int(f.stem.replace("flow", ""))
            w = load_omega(f)
            pv, loc = peak_near_le(w, X, Y)
            if best is None or abs(pv) > abs(best[0]):
                best = (pv, loc, step, w)
        peaks[dx] = best[:3]
        fields_at_peak_time[dx] = (best[3], best[2])
        print(f"dx={dx}: peak near-LE |omega|={best[0]:+.3f} at step {best[2]} (t={best[2]*dt:.2f})")

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 4, width_ratios=[2, 1.3, 2, 1.3], hspace=0.45, wspace=0.55,
                           top=0.86, bottom=0.06, left=0.05, right=0.98)
    for row, dx in enumerate(dxs):
        w, step = fields_at_peak_time[dx]
        X, Y = field_grid(dx)
        dt = dx / 2
        draw_field_with_zoom(fig, [gs[row, 0], gs[row, 1]], w, X, Y,
                              f"dx={dx}: peak-LE snapshot\n(t={step*dt:.1f})")
        # also show the t=0 (attached, fully developed by t~5-10 anyway) reference panel
        w0 = load_omega(DATA / f"gridconv_dx{dx}_snap" / f"flow{int(round(30.0/dt)):05d}.bin")
        draw_field_with_zoom(fig, [gs[row, 2], gs[row, 3]], w0, X, Y, f"dx={dx}: final state\n(t=30)")
    fig.suptitle("Test 1: leading-edge vorticity speck vs grid refinement (Re=500, $\\alpha$=0)\n"
                 "green box = zoom window (6x LE radius of curvature $r_{LE}$$\\approx$0.016); star = peak location",
                 fontsize=11)
    fig.savefig(FIGS / "fig1_grid_refinement.png", dpi=150)
    plt.close(fig)
    print(f"wrote {FIGS / 'fig1_grid_refinement.png'}")

    # add the dx=0.02 baseline point (reused from test 3) to the trend
    for dx02, rundir02 in extra_dirs.items():
        X02, Y02 = field_grid(dx02)
        avail02 = sorted(rundir02.glob("flow*.bin"))
        best02 = None
        for f in avail02:
            step = int(f.stem.replace("flow", ""))
            w = load_omega(f)
            pv, loc = peak_near_le(w, X02, Y02)
            if best02 is None or abs(pv) > abs(best02[0]):
                best02 = (pv, loc, step)
        peaks[dx02] = best02
        print(f"dx={dx02} (reused baseline): peak near-LE |omega|={best02[0]:+.3f} at step {best02[2]}")

    # summary: peak magnitude and distance-from-LE vs dx (include dx=0.02 baseline from test 3 if available)
    summary_dxs = sorted(peaks.keys(), reverse=True)
    fig2, ax = plt.subplots(1, 2, figsize=(9, 4))
    mags = [abs(peaks[dx][0]) for dx in summary_dxs]
    dists = [np.hypot(*peaks[dx][1]) for dx in summary_dxs]
    ax[0].plot(summary_dxs, mags, "o-", color="C0")
    ax[0].axvline(R_LE, color="gray", ls="--", lw=1, label=f"$r_{{LE}}$={R_LE:.4f}")
    ax[0].set_xlabel("dx"); ax[0].set_ylabel("peak $|\\omega|$ near LE"); ax[0].invert_xaxis()
    ax[0].set_title("LE peak vorticity magnitude vs dx"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(summary_dxs, dists, "o-", color="C3")
    ax[1].axvline(R_LE, color="gray", ls="--", lw=1, label=f"$r_{{LE}}$={R_LE:.4f}")
    ax[1].set_xlabel("dx"); ax[1].set_ylabel("distance of peak from LE (chord units)"); ax[1].invert_xaxis()
    ax[1].set_title("LE peak localization vs dx"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig2.suptitle("Test 1 summary: does the LE speck shrink & localize as dx < $r_{LE}$?")
    fig2.tight_layout(rect=(0, 0, 1, 0.92))
    fig2.savefig(FIGS / "fig1b_grid_refinement_summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"wrote {FIGS / 'fig1b_grid_refinement_summary.png'}")

    with open(HERE / "data" / "grid_refinement_peaks.txt", "w") as f:
        f.write("dx, peak|omega|_near_LE, peak_x, peak_y, dist_from_LE, step, t\n")
        for dx in summary_dxs:
            pv, (px, py), step = peaks[dx]
            dt = dx / 2
            f.write(f"{dx}, {pv:+.4f}, {px:+.4f}, {py:+.4f}, {np.hypot(px,py):.4f}, {step}, {step*dt:.2f}\n")


# ---------------------------------------------------------------------
# Test 2: alpha sweep, spatial top/bottom LE asymmetry
# ---------------------------------------------------------------------
def fig_alpha_sweep():
    alphas = [0, 2, 8, 10]
    dx = 0.02
    X, Y = field_grid(dx)
    step_final = 3000

    fig = plt.figure(figsize=(11, 9))
    gs = fig.add_gridspec(len(alphas), 2, width_ratios=[2, 1], hspace=0.4, wspace=0.2)
    asym = []
    for i, a in enumerate(alphas):
        rundir = DATA / f"alpha_a{a:+03.0f}_snap"
        fp = rundir / f"flow{step_final:05d}.bin"
        if not fp.exists():
            print(f"alpha sweep: missing {fp}, skipping"); return
        w = load_omega(fp)
        cf, peak_val, (px, py) = draw_field_with_zoom(fig, [gs[i, 0], gs[i, 1]], w, X, Y,
                                                        f"$\\alpha$={a}°, t=30")
        # top/bottom asymmetry: max |omega| in upper-LE window vs lower-LE window
        zoom_half = 6 * R_LE
        mask_top = le_window_mask(X, Y, zoom_half) & (Y > 0)
        mask_bot = le_window_mask(X, Y, zoom_half) & (Y < 0)
        top_peak = np.abs(w[mask_top]).max() if mask_top.any() else 0.0
        bot_peak = np.abs(w[mask_bot]).max() if mask_bot.any() else 0.0
        asym.append((a, top_peak, bot_peak))
        print(f"alpha={a}: top LE peak={top_peak:.3f}, bottom LE peak={bot_peak:.3f}, "
              f"asymmetry={(top_peak-bot_peak)/(top_peak+bot_peak):+.3f}")
    fig.suptitle(f"Test 2: LE vorticity field vs angle of attack (Re=500, dx={dx})\n"
                 "green box = LE zoom window; star = peak location", fontsize=11)
    fig.savefig(FIGS / "fig2_alpha_sweep.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIGS / 'fig2_alpha_sweep.png'}")

    a_vals = [r[0] for r in asym]
    top_vals = [r[1] for r in asym]
    bot_vals = [r[2] for r in asym]
    rel_asym = [(t - b) / (t + b) for t, b in zip(top_vals, bot_vals)]
    fig2, ax = plt.subplots(1, 2, figsize=(9, 4))
    ax[0].plot(a_vals, top_vals, "o-", label="top LE peak $|\\omega|$")
    ax[0].plot(a_vals, bot_vals, "s-", label="bottom LE peak $|\\omega|$")
    ax[0].set_xlabel(r"$\alpha$ (deg)"); ax[0].set_ylabel("peak $|\\omega|$ near LE")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3); ax[0].set_title("Top/bottom LE peaks vs $\\alpha$")
    ax[1].plot(a_vals, rel_asym, "o-", color="C2")
    ax[1].axhline(0, color="gray", lw=0.8)
    ax[1].set_xlabel(r"$\alpha$ (deg)"); ax[1].set_ylabel("relative asymmetry (top-bot)/(top+bot)")
    ax[1].grid(alpha=0.3); ax[1].set_title("LE peak asymmetry vs $\\alpha$")
    fig2.suptitle("Test 2 summary: does top/bottom LE asymmetry grow with angle of attack?")
    fig2.tight_layout(rect=(0, 0, 1, 0.92))
    fig2.savefig(FIGS / "fig2b_alpha_asymmetry_summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"wrote {FIGS / 'fig2b_alpha_asymmetry_summary.png'}")

    with open(HERE / "data" / "alpha_asymmetry.txt", "w") as f:
        f.write("alpha, top_LE_peak, bottom_LE_peak, relative_asymmetry\n")
        for a, t, b in asym:
            f.write(f"{a}, {t:+.4f}, {b:+.4f}, {(t-b)/(t+b):+.4f}\n")


# ---------------------------------------------------------------------
# Test 3: LE-densified boundary points vs uniform baseline (dx=0.02 fixed)
# ---------------------------------------------------------------------
def fig_le_densified():
    dx = 0.02
    X, Y = field_grid(dx)
    step_final = 3000
    cases = {"uniform (baseline)": DATA / "le_uniform_baseline_snap",
             "LE-densified boundary": DATA / "le_densified_snap"}
    fig = plt.figure(figsize=(9, 7))
    gs = fig.add_gridspec(2, 2, wspace=0.25, hspace=0.35)
    results = {}
    for i, (label, rundir) in enumerate(cases.items()):
        fp = rundir / f"flow{step_final:05d}.bin"
        if not fp.exists():
            print(f"LE-densified: missing {fp}, skipping"); return
        w = load_omega(fp)
        cf, peak_val, loc = draw_field_with_zoom(fig, [gs[i, 0], gs[i, 1]], w, X, Y, label)
        results[label] = peak_val
    fig.suptitle("Test 3: LE peak vorticity, uniform vs LE-densified boundary points\n"
                 "(background grid dx=0.02 held fixed in both -- isolates boundary-point density)", fontsize=10)
    fig.savefig(FIGS / "fig3_le_densified.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIGS / 'fig3_le_densified.png'}")
    print(f"peaks: {results}")
    with open(HERE / "data" / "le_densified_peaks.txt", "w") as f:
        for label, pv in results.items():
            f.write(f"{label}: peak|omega| near LE = {pv:+.4f}\n")


# ---------------------------------------------------------------------
# Test 4: spatial no-slip residual map
# ---------------------------------------------------------------------
def fig_residual_spatial():
    csv_path = HERE / "data" / "le_residual_spatial.csv"
    if not csv_path.exists():
        print("residual spatial: no data, skipping"); return
    d = np.genfromtxt(csv_path, delimiter=",", names=True)
    steps = sorted(set(d["step"]))

    fig, axes = plt.subplots(1, len(steps), figsize=(3.2 * len(steps), 3.4), sharey=True)
    if len(steps) == 1:
        axes = [axes]
    for ax, step in zip(axes, steps):
        row = d[d["step"] == step]
        order = np.argsort(row["s"])
        ax.semilogy(row["d_to_le"][order], np.maximum(row["res_mag"][order], 1e-18), "o", ms=3)
        ax.set_xlabel("arc-length distance from LE")
        ax.set_title(f"t={row['time'][0]:.1f}", fontsize=9)
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("no-slip residual $|C(\\omega)-b|$")
    fig.suptitle("Test 4: no-slip constraint residual vs distance from LE, Re=500 $\\alpha$=0 dx=0.02\n"
                 "(t=0 excluded: pre-solve initial-condition residual, not a solver accuracy measure)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(FIGS / "fig4_residual_vs_distance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIGS / 'fig4_residual_vs_distance.png'}")

    # spatial map (scatter over the body outline) at the last recorded step
    last = d[d["step"] == steps[-1]]
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    sc = ax2.scatter(last["x"], last["y"], c=np.maximum(last["res_mag"], 1e-18), cmap="viridis",
                      norm=matplotlib.colors.LogNorm(), s=40, edgecolor="k", linewidth=0.3)
    ax2.plot(0, 0, "r*", ms=12, label="LE (x=0,y=0)")
    ax2.set_aspect("equal"); ax2.set_xlabel("x"); ax2.set_ylabel("y")
    ax2.set_title(f"Test 4: spatial no-slip residual map, t={last['time'][0]:.1f}")
    ax2.legend(fontsize=8)
    fig2.colorbar(sc, ax=ax2, label="$|C(\\omega)-b|$ (log scale)")
    fig2.savefig(FIGS / "fig4b_residual_spatial_map.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"wrote {FIGS / 'fig4b_residual_spatial_map.png'}")


def main():
    fig_le_densified()
    fig_alpha_sweep()
    fig_residual_spatial()
    fig_grid_refinement()


if __name__ == "__main__":
    main()
