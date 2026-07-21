"""
gen_followup_figs.py

Figures for the 7 follow-up tests (A-G) in 2-follow_up/README.md. One figure
per test question (A/B/C combined into per-test line plots against alpha;
D/E combined into one "blockage" bar chart; F/G combined into one "alpha=0
offset" bar chart, since they're the same quantity under different knobs).

Usage: python3 SURF_test/low_re/NACA0012/kurt_comp/2-follow_up/gen_followup_figs.py
Output: 2-follow_up/figures/*.png
"""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures"
FIGS.mkdir(exist_ok=True)

# fixed categorical colors, used consistently across every figure --
# never reassigned/cycled per-figure
C_RAW = "#c0392b"       # red: the ORIGINAL 1-paper_based measurement
C_ALT = "#2980b9"       # blue: the follow-up's alternative measurement
C_ALT2 = "#16a085"      # teal: a third series where needed
C_REF = "#7f8c8d"       # gray: reference/target value (paper, pi, zero)


def load_csv(name):
    return np.genfromtxt(DATA / name, delimiter=",", names=True, dtype=None, encoding="utf-8")


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------- Test A
def fig_a_thrust_window():
    d = load_csv("thrust_window_reanalysis.csv")
    m = d["impl"] == "py"
    a, cd_min, cd_p05, cd_pa = d["alpha_deg"][m], d["cd_min"][m], d["cd_p05"][m], d["cd_phase_avg_min"][m]
    order = np.argsort(a)
    a, cd_min, cd_p05, cd_pa = a[order], cd_min[order], cd_p05[order], cd_pa[order]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    style_ax(ax)
    ax.axhline(0, color="#333333", linewidth=1)
    ax.axvspan(3, 37, color=C_REF, alpha=0.12, zorder=0, label="paper's thrust window (3–37°)")
    ax.plot(a, cd_min, color=C_RAW, linewidth=2, marker="o", markersize=3.5,
            label="raw minimum (1-paper_based's measurement)")
    ax.plot(a, cd_pa, color=C_ALT, linewidth=2, marker="o", markersize=3.5,
            label="phase-averaged minimum (this follow-up)")
    ax.set_xlabel("mean angle of attack α₀ (deg)")
    ax.set_ylabel("minimum instantaneous $C_d$ over the cycle")
    ax.set_title("Test A: is \"minimum instantaneous drag\" a fair thrust statistic?\n"
                  "f=4Hz pitching, py_static", fontsize=11)
    ax.legend(loc="lower right", fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "test_A_thrust_window.png", dpi=150)
    plt.close(fig)
    print("wrote test_A_thrust_window.png")


# ---------------------------------------------------------------- Test B
def fig_b_period_locked():
    d = load_csv("period_locked_mean_reanalysis.csv")
    m = d["impl"] == "py"
    a = d["alpha_deg"][m]
    cl_fixed = d["cl_mean_last50pct"][m]
    cl_locked = d["cl_mean_periodlocked"][m]
    order = np.argsort(a)
    a, cl_fixed, cl_locked = a[order], cl_fixed[order], cl_locked[order]
    band = (a >= 15) & (a <= 40)

    fig, ax = plt.subplots(figsize=(8.5, 5))
    style_ax(ax)
    ax.axvspan(15, 40, color=C_REF, alpha=0.12, zorder=0, label="post-stall region (15–40°)")
    ax.plot(a, cl_fixed, color=C_RAW, linewidth=2, marker="o", markersize=3.5,
            label="fixed last-50%-of-run window (1-paper_based)")
    ax.plot(a, cl_locked, color=C_ALT, linewidth=1.5, marker="x", markersize=5,
            linestyle="--", label="whole-shedding-period window (this follow-up)")
    ax.set_xlabel("mean angle of attack α₀ (deg)")
    ax.set_ylabel(r"mean lift coefficient $\overline{C_l}$ (steady)")
    ax.set_title("Test B: does the averaging window explain the post-stall jaggedness?\n"
                  "steady, py_static — the two curves nearly overlap: it doesn't", fontsize=11)
    ax.legend(loc="upper left", fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "test_B_period_locked_mean.png", dpi=150)
    plt.close(fig)
    print("wrote test_B_period_locked_mean.png")


# ---------------------------------------------------------------- Test C
def fig_c_strouhal_resolution():
    d = load_csv("strouhal_fine_reanalysis.csv")
    m = d["impl"] == "py"
    a, st_raw, st_fine = d["alpha_deg"][m], d["st_raw_bin"][m], d["st_zeropad_interp"][m]
    order = np.argsort(a)
    a, st_raw, st_fine = a[order], st_raw[order], st_fine[order]

    kurt_path = HERE.parents[0] / "1-paper_based" / "data" / "kurtulus_fig19_digitized.csv"
    # np.genfromtxt(..., names=True, comments='#') mis-parses this file's
    # '#'-prefixed description lines ABOVE the real header row (same issue
    # ../1-paper_based/gen_kurt_figs.py's read_csv_with_comments() works
    # around) -- strip full-line comments first, then parse the clean table
    import io as _io
    _lines = [l for l in kurt_path.read_text().splitlines()
              if l.strip() and not l.lstrip().startswith("#")]
    kurt = np.genfromtxt(_io.StringIO("\n".join(_lines)), delimiter=",",
                          names=True, dtype=None, encoding="utf-8")

    fig, ax = plt.subplots(figsize=(8.5, 5))
    style_ax(ax)
    # both ibpm series use the SAME drawing method (direct point-to-point
    # line, no drawstyle="steps-mid") -- an earlier version of this figure
    # rendered the raw-bin series with steps-mid and the fine series with a
    # direct line, which manufactured visual "stair-steps" for the raw
    # series that were partly a rendering choice, not purely the data (see
    # README's correction note in the Test C section)
    ax.plot(a, st_raw, color=C_RAW, linewidth=1.5, marker="o", markersize=3,
            label="raw FFT bin (1-paper_based)")
    ax.plot(a, st_fine, color=C_ALT, linewidth=2, marker="o", markersize=3,
            label="zero-padded + interpolated (this follow-up)")
    ax.plot(kurt["alpha_deg"], kurt["strouhal"], color=C_REF, linewidth=2, marker="s",
            markersize=3, linestyle="-", label="Kurtulus (2019) Fig 19 (digitized)")
    ax.set_xlabel("mean angle of attack α₀ (deg)")
    ax.set_ylabel("vortex-shedding Strouhal number (steady)")
    ax.set_title("Test C: finer frequency resolution vs. the paper's own curve\n"
                  "steady, py_static", fontsize=11)
    ax.legend(loc="upper right", fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "test_C_strouhal_resolution.png", dpi=150)
    plt.close(fig)
    print("wrote test_C_strouhal_resolution.png")


# ------------------------------------------------------------- Tests D+E
def fig_de_blockage():
    dn = load_csv("ngrid_sweep_reanalysis.csv")
    dd = load_csv("domain_sweep_reanalysis.csv")

    def slope(alpha, cl):
        order = np.argsort(alpha)
        return np.polyfit(np.radians(alpha[order]), cl[order], 1)[0]

    labels, slopes, colors = [], [], []
    for ngrid, label in [(1, "ngrid=1\n(baseline)"), (2, "ngrid=2"), (3, "ngrid=3")]:
        m = dn["ngrid"] == ngrid
        labels.append(label)
        slopes.append(slope(dn["alpha_deg"][m], dn["cl_mean"][m]))
        colors.append(C_RAW if ngrid == 1 else C_ALT)
    for dom, label in [("baseline_L6", None), ("large_L10", "larger\ndomain (L=10)")]:
        if label is None:
            continue
        m = dd["domain"] == dom
        labels.append(label)
        slopes.append(slope(dd["alpha_deg"][m], dd["cl_mean"][m]))
        colors.append(C_ALT2)

    fig, ax = plt.subplots(figsize=(8, 5.2))
    style_ax(ax)
    x = np.arange(len(labels))
    ax.axhline(np.pi, color="#333333", linewidth=1.3, linestyle="--", zorder=1)
    ax.text(len(labels) - 0.55, np.pi + 0.03, "π (textbook thin-airfoil value)",
            fontsize=8.5, color="#333333", ha="right", va="bottom")
    bars = ax.bar(x, slopes, color=colors, width=0.6, zorder=2)
    for xi, s in zip(x, slopes):
        ax.text(xi, s + 0.03, f"{s:.2f}", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(r"lift-curve slope $dC_l/d\alpha$ (rad$^{-1}$)")
    ax.set_ylim(0, max(slopes) * 1.18)
    ax.set_title("Tests D & E: does more far-field domain reduce the lift-slope excess?\n"
                  "steady, α=0–5°, py_static, dx=0.02", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "test_DE_blockage_liftslope.png", dpi=150)
    plt.close(fig)
    print("wrote test_DE_blockage_liftslope.png")


# ------------------------------------------------------------- Tests F+G
def fig_fg_alpha0_offset():
    labels = ["dx=0.02\n(baseline)", "dx=0.02,\nhalf-cell shifted", "dx=0.01\n(refined)"]
    values = [-0.00648, -0.02325, -0.00089]
    colors = [C_RAW, C_ALT, C_ALT2]

    fig, ax = plt.subplots(figsize=(7, 5.2))
    style_ax(ax)
    ax.axhline(0, color="#333333", linewidth=1.3, zorder=1)
    x = np.arange(len(labels))
    ax.bar(x, values, color=colors, width=0.55, zorder=2)
    for xi, v in zip(x, values):
        ax.text(xi, v - 0.0016 if v < 0 else v + 0.0009, f"{v:+.5f}",
                ha="center", va="top" if v < 0 else "bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(r"lift coefficient $C_l$ at $\alpha$=0° (should be exactly 0)")
    ax.set_title("Tests F & G: is the nonzero $C_l(0)$ a grid artifact?\n"
                  "steady, py_static", fontsize=11)
    ax.set_ylim(min(values) * 1.35, max(0, max(values)) + 0.003)
    fig.tight_layout()
    fig.savefig(FIGS / "test_FG_alpha0_offset.png", dpi=150)
    plt.close(fig)
    print("wrote test_FG_alpha0_offset.png")


if __name__ == "__main__":
    fig_a_thrust_window()
    fig_b_period_locked()
    fig_c_strouhal_resolution()
    fig_de_blockage()
    fig_fg_alpha0_offset()
