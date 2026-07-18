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
def fig1314_hysteresis():
    tabp = DATA / "fig13_14_comparison.csv"
    kt = read_csv_with_comments(DATA / "kurtulus_fig13_14_table.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    f = FREQ["f4hz"]; period = 1.0 / f
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
    # paper points
    for br, mk in [("down", "v"), ("up", "^")]:
        mb = kt["branch"] == br
        axes[0].plot(kt["alpha_deg"][mb], kt["cl"][mb], mk, color="k", mfc="0.7", ms=7,
                     ls="none", label=f"Kurtulus ({br})")
        axes[1].plot(kt["alpha_deg"][mb], kt["cd"][mb], mk, color="k", mfc="0.7", ms=7,
                     ls="none")
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
    motions = [("steady", "steady"), ("f4hz", "f4hz (t=dev)")]
    for mot, mlabel in motions:
        fig, axes = plt.subplots(len(angles), 2, figsize=(11, 2.4 * len(angles)))
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
        fig.suptitle(f"Wake vorticity, {mlabel}, NACA0012 Re=1000 "
                     f"(jet: blue=-, green=0, red=+, matching Kurtulus Figs)", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out = FIGS / f"wake_{mot}.png"
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
