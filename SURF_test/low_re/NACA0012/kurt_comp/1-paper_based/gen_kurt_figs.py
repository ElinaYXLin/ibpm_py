"""
gen_kurt_figs.py

Comparison figures for the Kurtulus (2019) reproduction. Reads the reduced
tables from analyze_kurt.py (kurt_comp/data/) plus the raw runs, and the
digitized paper data, and writes figures into kurt_comp/figures/.

All of Kurtulus's data is graphical (no tables), so per the task each figure
shows the paper's data redrawn in its own grayscale style side-by-side with the
IBPM (py_static / cpp_static) result on matched axes. Vorticity fields use a
jet colormap (blue=negative, green~0, red=positive) matching the paper's Figs.

Usage: python3 SURF_test/low_re/NACA0012/kurt_comp/gen_kurt_figs.py
"""
import io
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "1-paper_based"
RUNS = KURT / "runs"
DATA = KURT / "data"
FIGS = KURT / "figures"
FIGS.mkdir(exist_ok=True)


def read_csv_with_comments(path):
    """np.genfromtxt(..., names=True, comments='#') mis-parses files that
    have '#'-prefixed description lines ABOVE a real header row (it tries to
    use an early comment line's shape for column count) -- strip full-line
    comments ourselves first, then hand genfromtxt a clean names=True table."""
    lines = [l for l in pathlib.Path(path).read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    return np.genfromtxt(io.StringIO("\n".join(lines)), delimiter=",",
                         names=True, dtype=None, encoding="utf-8")

sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py_static")
pkg.__path__ = [str(REPO / "py_static")]
sys.modules["py_static"] = pkg
from py_static.state import State  # noqa: E402

FREQ = {"f1hz": 0.684931506849315, "f4hz": 2.73972602739726}
IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}
MOTION_LS = {"steady": "-", "f1hz": "--", "f4hz": ":"}
GRID = "dx0.020"


def load_mc(grid=GRID):
    p = DATA / f"mean_coefficients_{grid}.csv"
    if not p.exists():
        return None
    return np.genfromtxt(p, delimiter=",", names=True,
                         dtype=None, encoding="utf-8")


def sub(mc, motion, impl):
    m = (mc["motion"] == motion) & (mc["impl"] == impl) & (mc["nan"] == 0)
    a = mc["alpha_deg"][m]
    o = np.argsort(a)
    return a[o], mc["cl_mean"][m][o], mc["cd_mean"][m][o]


# ---------------------------------------------------------------- Figure 1
def fig1_mean():
    mc = load_mc()
    if mc is None:
        print("fig1: no data yet"); return
    k = read_csv_with_comments(DATA / "kurtulus_fig1_digitized.csv")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    # --- left column: Kurtulus paper style ---
    for row, (qty, klc, kdc, ylab) in enumerate([
            ("cl", ["cl_steady", "cl_f1hz", "cl_f4hz"], None, r"$\overline{C_l}$"),
            ("cd", ["cd_steady", "cd_f1hz", "cd_f4hz"], None, r"$\overline{C_d}$")]):
        axK = axes[row, 0]
        styles = [("steady", "k-", "^", "white"),
                  ("f1hz", "k-", "o", "white"),
                  ("f4hz", "k-", "o", "0.6")]
        for (mot, ls, mk, mfc), col in zip(styles, klc):
            axK.plot(k["alpha_deg"], k[col], ls, marker=mk, mfc=mfc, mec="k",
                     ms=5, lw=0.8, label=f"NACA0012 {mot}")
        axK.set_title(f"Kurtulus (2019) Fig 1 (digitized) -- {ylab}", fontsize=10)
        axK.set_ylabel(ylab); axK.set_xlabel(r"$\alpha_0$ [deg]")
        axK.legend(fontsize=7); axK.grid(alpha=0.3)
        # --- right column: IBPM ---
        axI = axes[row, 1]
        for mot in ("steady", "f1hz", "f4hz"):
            for impl in ("py", "cpp"):
                a, cl, cd = sub(mc, mot, impl)
                y = cl if qty == "cl" else cd
                axI.plot(a, y, MOTION_LS[mot], color=IMPL_COLOR[impl], lw=1.3,
                         marker="." if impl == "py" else "x", ms=4,
                         label=f"{mot} {IMPL_LABEL[impl]}")
        axI.set_title(f"IBPM (this solver) -- {ylab}", fontsize=10)
        axI.set_xlabel(r"$\alpha_0$ [deg]"); axI.grid(alpha=0.3)
        axI.legend(fontsize=6, ncol=2)
        lo = min(axK.get_ylim()[0], axI.get_ylim()[0])
        hi = max(axK.get_ylim()[1], axI.get_ylim()[1])
        axK.set_ylim(lo, hi); axI.set_ylim(lo, hi)
    fig.suptitle("Mean aerodynamic coefficients, NACA0012 Re=1000: "
                 "Kurtulus CFD vs this immersed-boundary solver", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIGS / "fig1_mean_coefficients.png", dpi=140)
    plt.close(fig)
    print("wrote fig1_mean_coefficients.png")


# ---------------------------------------------------------------- Figure 19
def fig19_shedding():
    p = DATA / f"shedding_strouhal_{GRID}.csv"
    if not p.exists():
        print("fig19: no data yet"); return
    st = np.genfromtxt(p, delimiter=",", names=True, dtype=None, encoding="utf-8")
    k = read_csv_with_comments(DATA / "kurtulus_fig19_digitized.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    axes[0].plot(k["alpha_deg"], k["strouhal"], "k-o", mfc="0.6", ms=5, lw=0.9,
                 label="Kurtulus Fig 19 (digitized)")
    axes[0].set_title("Kurtulus (2019) Fig 19 (digitized), as Strouhal", fontsize=10)
    for impl in ("py", "cpp"):
        m = st["impl"] == impl
        a = st["alpha_deg"][m]; o = np.argsort(a)
        axes[1].plot(a[o], st["strouhal_nondim"][m][o], "-o", color=IMPL_COLOR[impl],
                     ms=4, lw=1.2, label=IMPL_LABEL[impl])
    axes[1].set_title("IBPM (this solver): non-dim shedding freq = Strouhal", fontsize=10)
    for ax in axes:
        ax.set_xlabel(r"$\alpha_0$ [deg]"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        ax.axvspan(0, 7.5, color="0.9", zorder=0)
    axes[0].set_ylabel("Strouhal number $St = f c / U$")
    fig.suptitle("Vortex-shedding Strouhal number vs angle of attack, steady "
                 "NACA0012 Re=1000 (shaded: paper reports no shedding below ~8 deg)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIGS / "fig19_shedding_strouhal.png", dpi=140)
    plt.close(fig)
    print("wrote fig19_shedding_strouhal.png")


# ------------------------------------------------ Figure 11: instantaneous
def fig11_instantaneous():
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for col, mot in enumerate(("f1hz", "f4hz")):
        f = FREQ[mot]
        period = 1.0 / f
        for impl in ("py", "cpp"):
            run = RUNS / GRID / f"{mot}_{impl}_a00"
            fp = run / "flow.force"
            if not fp.exists():
                continue
            d = np.loadtxt(fp)
            t, cd, cl = d[:, 1], d[:, 2], d[:, 3]
            # last 2 periods, re-zero time
            tmax = t[-1]
            m = t > (tmax - 2 * period)
            tt = (t[m] - t[m][0]) / period
            axes[0, col].plot(tt, cl[m], color=IMPL_COLOR[impl], lw=1.3,
                              label=IMPL_LABEL[impl])
            axes[1, col].plot(tt, cd[m], color=IMPL_COLOR[impl], lw=1.3,
                              label=IMPL_LABEL[impl])
        axes[0, col].set_title(f"pitching {mot}, $\\alpha_0=0$", fontsize=10)
        axes[0, col].set_ylabel("$C_l$"); axes[1, col].set_ylabel("$C_d$")
        axes[1, col].set_xlabel("t / pitch period")
        axes[1, col].axhline(0, color="0.5", lw=0.7, ls=":")
        for r in (0, 1):
            axes[r, col].grid(alpha=0.3); axes[r, col].legend(fontsize=8)
    fig.suptitle("Instantaneous lift/drag over 2 pitch periods (developed state), "
                 "NACA0012 Re=1000, $\\alpha_0=0$ (cf. Kurtulus Fig 11)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGS / "fig11_instantaneous_pitch.png", dpi=140)
    plt.close(fig)
    print("wrote fig11_instantaneous_pitch.png")


# --------------------------------- Figure 13/14: hysteresis loop vs table
def _hysteresis_error_metrics(kt, t, cd_m, cl_m, alpha, f):
    """Peak-to-peak amplitude ratio (ibpm/paper) and RMS absolute error
    (ibpm interpolated onto the paper's own alpha, matched by branch so the
    double-valued loop isn't averaged across both halves), for both Cl and
    Cd. Written out and printed so a "which one matches better" claim in
    the README is checked against numbers instead of a visual impression --
    see kurt_comp/1-paper_based/README.md's "Instantaneous forces and
    hysteresis" section for why this matters (a prior version of that
    section had Cl/Cd's relative agreement backwards).
    """
    dadt = np.cos(2 * np.pi * f * t)
    rows = []
    for name, paper_col, arr in [("Cl", "cl", cl_m), ("Cd", "cd", cd_m)]:
        errs = []
        for i in range(len(kt)):
            a0 = kt["alpha_deg"][i]
            br = kt["branch"][i]
            mask = (dadt > 0) if br == "up" else (dadt < 0)
            if mask.sum() < 2:
                continue
            idx = np.argsort(alpha[mask])
            val = np.interp(a0, alpha[mask][idx], arr[mask][idx])
            errs.append(val - kt[paper_col][i])
        errs = np.array(errs)
        rms = float(np.sqrt(np.mean(errs ** 2)))
        paper_range = float(kt[paper_col].max() - kt[paper_col].min())
        ibpm_range = float(arr.max() - arr.min())
        rows.append(dict(coef=name, ibpm_pk2pk=ibpm_range, paper_pk2pk=paper_range,
                          pk2pk_ratio=ibpm_range / paper_range,
                          rms_abs_error=rms, rms_error_pct_of_paper_range=100 * rms / paper_range))
    return rows


def fig1314_hysteresis():
    tabp = DATA / "fig13_14_comparison.csv"
    kt = read_csv_with_comments(DATA / "kurtulus_fig13_14_table.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    f = FREQ["f4hz"]; period = 1.0 / f
    metrics = None
    for impl in ("py", "cpp"):
        run = RUNS / GRID / f"f4hz_{impl}_a00"
        fp = run / "flow.force"
        if not fp.exists():
            continue
        d = np.loadtxt(fp)
        t, cd, cl = d[:, 1], d[:, 2], d[:, 3]
        m = t > (t[-1] - 2 * period)
        alpha = np.sin(2 * np.pi * f * t[m])  # deg, amplitude 1
        axes[0].plot(alpha, cl[m], color=IMPL_COLOR[impl], lw=1.2, label=IMPL_LABEL[impl])
        axes[1].plot(alpha, cd[m], color=IMPL_COLOR[impl], lw=1.2, label=IMPL_LABEL[impl])
        if impl == "py":  # py/cpp are identical here (see 5-leading_edge's Test 1a); one is enough
            metrics = _hysteresis_error_metrics(kt, t[m], cd[m], cl[m], alpha, f)
    # paper points
    for br, mk in [("down", "v"), ("up", "^")]:
        mb = kt["branch"] == br
        axes[0].plot(kt["alpha_deg"][mb], kt["cl"][mb], mk, color="k", mfc="0.7", ms=7,
                     ls="none", label=f"Kurtulus ({br})")
        axes[1].plot(kt["alpha_deg"][mb], kt["cd"][mb], mk, color="k", mfc="0.7", ms=7,
                     ls="none")
    if metrics:
        with open(DATA / "hysteresis_error_metrics.csv", "w", newline="") as fh:
            import csv as _csv
            w = _csv.DictWriter(fh, fieldnames=list(metrics[0].keys()))
            w.writeheader(); w.writerows(metrics)
        by_coef = {r["coef"]: r for r in metrics}
        axes[0].set_ylabel("$C_l$")
        axes[1].set_ylabel("$C_d$")
        axes[0].set_title(f"pk-pk: ibpm/paper = {by_coef['Cl']['pk2pk_ratio']:.2f}x "
                          f"(+{100*(by_coef['Cl']['pk2pk_ratio']-1):.0f}%), "
                          f"RMS err = {by_coef['Cl']['rms_error_pct_of_paper_range']:.0f}% of paper range", fontsize=9)
        axes[1].set_title(f"pk-pk: ibpm/paper = {by_coef['Cd']['pk2pk_ratio']:.2f}x "
                          f"(+{100*(by_coef['Cd']['pk2pk_ratio']-1):.0f}%), "
                          f"RMS err = {by_coef['Cd']['rms_error_pct_of_paper_range']:.0f}% of paper range", fontsize=9)
        print(f"  Cl: pk2pk ratio={by_coef['Cl']['pk2pk_ratio']:.3f}, "
              f"RMS error={by_coef['Cl']['rms_error_pct_of_paper_range']:.1f}% of paper range")
        print(f"  Cd: pk2pk ratio={by_coef['Cd']['pk2pk_ratio']:.3f}, "
              f"RMS error={by_coef['Cd']['rms_error_pct_of_paper_range']:.1f}% of paper range")
    else:
        axes[0].set_ylabel("$C_l$"); axes[1].set_ylabel("$C_d$")
    for ax in axes:
        ax.set_xlabel(r"instantaneous $\alpha$ [deg]"); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.suptitle("Dynamic hysteresis loop, f=4Hz pitching, $\\alpha_0=0$, NACA0012 "
                 "Re=1000: IBPM (lines) vs Kurtulus Fig 13/14 points", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIGS / "fig13_14_hysteresis.png", dpi=140)
    plt.close(fig)
    print("wrote fig13_14_hysteresis.png")


# ---------------------------------------------- thrust (negative Cd) check
def thrust_fig():
    p = DATA / f"thrust_check_{GRID}.csv"
    if not p.exists():
        print("thrust: no data yet"); return
    th = np.genfromtxt(p, delimiter=",", names=True, dtype=None, encoding="utf-8")
    fig, ax = plt.subplots(figsize=(9, 5))
    for impl in ("py", "cpp"):
        m = th["impl"] == impl
        a = th["alpha_deg"][m]; o = np.argsort(a)
        ax.plot(a[o], th["cd_min"][m][o], "-o", color=IMPL_COLOR[impl], ms=4,
                label=IMPL_LABEL[impl])
    ax.axhline(0, color="k", lw=0.8)
    ax.axvspan(3, 37, color="orange", alpha=0.15,
               label="Kurtulus: thrust (min $C_d<0$) here")
    ax.set_xlabel(r"$\alpha_0$ [deg]"); ax.set_ylabel("min instantaneous $C_d$")
    ax.set_title("Minimum instantaneous drag over a cycle, f=4Hz pitching, "
                 "NACA0012 Re=1000\n(paper: $C_d<0$ = thrust for $3\\leq\\alpha_0\\leq37$)",
                 fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "thrust_check.png", dpi=140)
    plt.close(fig)
    print("wrote thrust_check.png")


# ------------------------------------------------------ wake vorticity fields
def wake_contours():
    NX, NY = 300, 150
    dx = 6.0 / NX
    xs = -2 + np.arange(1, NX) * dx
    ys = -1.5 + np.arange(1, NY) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.genfromtxt(REPO / "SURF_test" / "low_re" / "NACA0012" / "1-basics" /
                        "naca0012.dat.txt", skip_header=1)
    angles = [0, 9, 12]
    # third column = the paper's own figure, cropped by extract_paper_figs.py
    # (Figure 2 for steady, Figure 6 for f4hz -- see paper_figs/README.md)
    paper_dir = pathlib.Path(__file__).resolve().parent / "paper_figs"
    motions = [("steady", "steady"), ("f4hz", "f4hz (t=dev)")]
    for mot, mlabel in motions:
        fig, axes = plt.subplots(len(angles), 3, figsize=(15.5, 2.4 * len(angles)))
        axes = np.atleast_2d(axes)
        for r, ang in enumerate(angles):
            for c, impl in enumerate(("py", "cpp")):
                ax = axes[r, c]
                run = RUNS / GRID / f"{mot}_{impl}_a{ang:02d}"
                snaps = sorted(run.glob("flow[0-9]*.bin"))
                if not snaps:
                    ax.text(0.5, 0.5, "no snapshot", ha="center", va="center")
                    ax.set_xticks([]); ax.set_yticks([]); continue
                w = State(filename=str(snaps[-1])).omega._data[0]
                V = 8.0
                ax.contourf(X, Y, np.clip(w, -V, V), levels=41, cmap="jet",
                            extend="both")
                ax.fill(pts[:, 0], pts[:, 1], color="0.1", zorder=5)
                ax.set_xlim(-1, 4); ax.set_ylim(-1.2, 1.2); ax.set_aspect("equal")
                ax.set_title(f"{IMPL_LABEL[impl]}, $\\alpha_0$={ang}deg", fontsize=9)
            # third column: what we're actually comparing against
            ax = axes[r, 2]
            paper_png = paper_dir / f"{mot}_a{ang:02d}.png"
            if paper_png.exists():
                ax.imshow(plt.imread(paper_png))
                ax.set_title(f"Kurtulus (2019), $\\alpha_0$={ang}deg", fontsize=9)
            else:
                ax.text(0.5, 0.5, "paper crop missing", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
        fig.suptitle(f"Wake vorticity, {mlabel}, NACA0012 Re=1000 -- ibpm (py_static/cpp_static) "
                     f"vs. the paper's own figure\n(jet: blue=-, green=0, red=+, matching Kurtulus Figs; "
                     f"third column cropped directly from Kurtulus 2019 Fig. {'2' if mot == 'steady' else '6'})",
                     fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        out = FIGS / f"wake_{mot}.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"wrote {out.name}")


# ---------------------------------------- wake vorticity, paper's own frame
def wake_contours_paperframe():
    """Same wake data as wake_contours(), but rotated into the frame the
    paper actually plots in.

    The solver imposes alpha0 by rotating the *free-stream* to angle alpha0
    relative to the fixed grid/body (BaseFlow builds Flux.UniformFlow(grid,
    mag, alpha), and ibpm.py's drag/lift split -- drag = Fx*cos(a)+Fy*sin(a),
    lift = -Fx*sin(a)+Fy*cos(a) -- confirms (cos a, sin a) is the free-stream
    direction in grid coordinates); the .geom file itself is never rotated
    per alpha0 (run_kurt_suite.py's ensure_geom reuses the same raw geometry
    for every angle). Kurtulus's own figures do the opposite: free-stream
    horizontal, body pitched to alpha0. Both are the same flow field up to a
    rigid rotation of the whole picture by alpha0, so we undo it here by
    applying R(-alpha0) to the grid coordinates and the airfoil outline
    (vorticity is a scalar and doesn't transform) before plotting -- this
    should reproduce the paper's layout if (and only if) the discrepancy the
    mentor flagged is a plotting convention, not a solver difference.
    """
    NX, NY = 300, 150
    dx = 6.0 / NX
    xs = -2 + np.arange(1, NX) * dx
    ys = -1.5 + np.arange(1, NY) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.genfromtxt(REPO / "SURF_test" / "low_re" / "NACA0012" / "1-basics" /
                        "naca0012.dat.txt", skip_header=1)
    angles = [0, 9, 12]
    paper_dir = pathlib.Path(__file__).resolve().parent / "paper_figs"
    motions = [("steady", "steady"), ("f4hz", "f4hz (t=dev)")]
    for mot, mlabel in motions:
        fig, axes = plt.subplots(len(angles), 3, figsize=(15.5, 3.2 * len(angles)))
        axes = np.atleast_2d(axes)
        for r, ang in enumerate(angles):
            a = ang * np.pi / 180.0
            ca, sa = np.cos(a), np.sin(a)
            # R(-alpha0): undoes the solver's flow-rotation convention so the
            # free-stream direction (cos a, sin a) in grid coords becomes (1, 0)
            Xr = X * ca + Y * sa
            Yr = -X * sa + Y * ca
            ptsr = np.column_stack([pts[:, 0] * ca + pts[:, 1] * sa,
                                     -pts[:, 0] * sa + pts[:, 1] * ca])
            for c, impl in enumerate(("py", "cpp")):
                ax = axes[r, c]
                run = RUNS / GRID / f"{mot}_{impl}_a{ang:02d}"
                snaps = sorted(run.glob("flow[0-9]*.bin"))
                if not snaps:
                    ax.text(0.5, 0.5, "no snapshot", ha="center", va="center")
                    ax.set_xticks([]); ax.set_yticks([]); continue
                w = State(filename=str(snaps[-1])).omega._data[0]
                V = 8.0
                ax.contourf(Xr, Yr, np.clip(w, -V, V), levels=41, cmap="jet",
                            extend="both")
                ax.fill(ptsr[:, 0], ptsr[:, 1], color="0.1", zorder=5)
                ax.set_xlim(-1, 4); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
                ax.set_title(f"{IMPL_LABEL[impl]}, $\\alpha_0$={ang}deg "
                             f"(rotated -{ang}deg to paper frame)", fontsize=9)
            ax = axes[r, 2]
            paper_png = paper_dir / f"{mot}_a{ang:02d}.png"
            if paper_png.exists():
                ax.imshow(plt.imread(paper_png))
                ax.set_title(f"Kurtulus (2019), $\\alpha_0$={ang}deg", fontsize=9)
            else:
                ax.text(0.5, 0.5, "paper crop missing", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
        fig.suptitle(f"Wake vorticity, {mlabel}, NACA0012 Re=1000 -- ibpm rotated into the "
                     f"paper's plotting frame (free-stream horizontal, body pitched)\n"
                     f"(jet: blue=-, green=0, red=+; third column cropped directly from "
                     f"Kurtulus 2019 Fig. {'2' if mot == 'steady' else '6'})",
                     fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        out = FIGS / f"wake_{mot}_paperframe.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"wrote {out.name}")


if __name__ == "__main__":
    fig1_mean()
    fig19_shedding()
    fig11_instantaneous()
    fig1314_hysteresis()
    thrust_fig()
    wake_contours()
    wake_contours_paperframe()
