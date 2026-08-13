"""
gen_further_figs.py

Figures for the 7 further tests (1a/1b/2a/3a zero-new-run; 2b/2c/3b new-run)
in 3-further/README.md.

Usage: python3 gen_further_figs.py
Output: 3-further/figures/*.png
"""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures"
FIGS.mkdir(exist_ok=True)

C_A = "#c0392b"
C_B = "#2980b9"
C_C = "#16a085"
C_D = "#8e44ad"
C_REF = "#7f8c8d"


def load_csv(name):
    return np.genfromtxt(DATA / name, delimiter=",", names=True, dtype=None, encoding="utf-8")


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def fig_1a():
    d = load_csv("test1a_full_phase_average.csv")
    a = d["alpha_deg"]
    order = np.argsort(a)
    a = a[order]
    m4 = d["cd_phase_avg_min_4period"][order]
    mf = d["cd_phase_avg_min_fullrecord"][order]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    style_ax(ax)
    ax.axhline(0, color="#333333", linewidth=1)
    ax.plot(a, m4, color=C_A, linewidth=2, marker="o", markersize=4,
            label="phase-avg over last 4 periods (2-follow_up's Test A)")
    ax.plot(a, mf, color=C_B, linewidth=2, marker="s", markersize=4,
            label="phase-avg over the full developed record (this test)")
    ax.set_xlabel("mean angle of attack α₀ (deg)")
    ax.set_ylabel("minimum phase-averaged $C_d$")
    ax.set_title("Test 1a: does more averaging remove the residual dip at 34-39°?\n"
                 "f=4Hz pitching, py_static", fontsize=11)
    ax.legend(loc="upper left", fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "test_1a_full_phase_average.png", dpi=150)
    plt.close(fig)
    print("wrote test_1a_full_phase_average.png")


def fig_1b():
    d = load_csv("test1b_dip_vs_transitions.csv")
    a = d["alpha_deg"]
    order = np.argsort(a)
    a = a[order]
    dip = d["thrust_dip_cd"][order]
    st = d["steady_strouhal_fine"][order]
    cl = d["steady_cl_mean"][order]

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 8), sharex=True)
    for ax, y, color, label, ylab in [
        (axes[0], dip, C_A, "thrust dip (f4Hz, phase-avg min $C_d$)", "$C_d$"),
        (axes[1], st, C_B, "steady vortex-shedding Strouhal", "St"),
        (axes[2], cl, C_C, "steady mean lift", "$\\overline{C_l}$"),
    ]:
        style_ax(ax)
        ax.plot(a, y, color=color, linewidth=2, marker="o", markersize=4)
        ax.set_ylabel(ylab, fontsize=9.5)
        ax.axvline(34.5, color=C_REF, linewidth=1.2, linestyle="--", zorder=1)
        ax.set_title(label, fontsize=9.5, loc="left")
    axes[0].axhline(0, color="#333333", linewidth=1)
    axes[-1].set_xlabel("mean angle of attack α₀ (deg)")
    fig.suptitle("Test 1b: three signals, same angle axis -- do they line up at 34-35°?",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGS / "test_1b_dip_vs_transitions.png", dpi=150)
    plt.close(fig)
    print("wrote test_1b_dip_vs_transitions.png")


def fig_2a():
    d = np.load(DATA / "test2a_spectrogram.npz")
    angles, freq_grid, amp_map = d["angles"], d["freq_grid"], d["amp_map"]
    keep = angles <= 45
    angles, amp_map = angles[keep], amp_map[keep]

    cmap = LinearSegmentedColormap.from_list("amp", ["#ffffff", "#2980b9", "#16a085", "#c0392b"])
    fig, ax = plt.subplots(figsize=(9, 6))
    im = ax.pcolormesh(angles, freq_grid, np.nan_to_num(amp_map).T, cmap=cmap,
                        shading="nearest", vmin=0, vmax=1)
    ax.set_xlabel("mean angle of attack α₀ (deg)")
    ax.set_ylabel("Strouhal number St")
    ax.set_ylim(0, 1.0)
    ax.set_title("Test 2a: full spectral content vs. angle (steady, py_static)\n"
                 "not just the single dominant peak -- normalized amplitude per angle",
                 fontsize=11)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("FFT amplitude (normalized to that angle's own peak)")
    fig.tight_layout()
    fig.savefig(FIGS / "test_2a_spectrogram.png", dpi=150)
    plt.close(fig)
    print("wrote test_2a_spectrogram.png")


def fig_3a():
    d = load_csv("test3a_running_mean_convergence.csv")
    a = d["alpha_deg"]
    order = np.argsort(a)
    a = a[order]
    pct = d["pct_change_last_third_vs_first_two_thirds"][order]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    style_ax(ax)
    ax.axhline(0, color="#333333", linewidth=1)
    ax.axhspan(-5, 5, color=C_REF, alpha=0.15, zorder=0, label="within ±5% (roughly converged)")
    colors = [C_A if abs(p) > 10 else (C_D if abs(p) > 5 else C_C) for p in pct]
    ax.bar(a, pct, color=colors, width=0.8, zorder=2)
    ax.set_xlabel("mean angle of attack α₀ (deg)")
    ax.set_ylabel("% change: last-third-of-run mean vs. first-two-thirds mean")
    ax.set_title("Test 3a: has the mean $C_l$ actually converged by t=30?\n"
                 "steady, py_static -- large bars mean the run is still drifting",
                 fontsize=11)
    ax.legend(loc="upper right", fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "test_3a_running_mean_convergence.png", dpi=150)
    plt.close(fig)
    print("wrote test_3a_running_mean_convergence.png")


def _paper_strouhal_at(angles):
    """Interpolate Kurtulus (2019) Fig 19's digitized Strouhal-vs-alpha
    curve at the given angles, for a "ground truth" reference line on the
    2b/2c bar charts. Same file used by 1-paper_based's fig19_shedding()
    and 2-follow_up's fig_c_strouhal_resolution()."""
    kurt_path = HERE.parents[0] / "1-paper_based" / "data" / "kurtulus_fig19_digitized.csv"
    lines = [l for l in kurt_path.read_text().splitlines()
              if l.strip() and not l.lstrip().startswith("#")]
    import io as _io
    kurt = np.genfromtxt(_io.StringIO("\n".join(lines)), delimiter=",",
                          names=True, dtype=None, encoding="utf-8")
    order = np.argsort(kurt["alpha_deg"])
    return np.interp(angles, kurt["alpha_deg"][order], kurt["strouhal"][order])


def fig_2b_2c():
    d2b = load_csv("test2b_dx_refine_strouhal.csv")
    d2c = load_csv("test2c_ngrid_strouhal.csv")
    angles = [15, 20, 30, 40]
    paper_st = _paper_strouhal_at(angles)

    fig, axes = plt.subplots(1, 4, figsize=(13, 4.3), sharey=False)
    for i, a in enumerate(angles):
        ax = axes[i]
        style_ax(ax)
        row2b = d2b[d2b["alpha_deg"] == a]
        col005 = row2b["st_dx0005"]
        st005 = (float(col005[0]) if len(row2b) and np.issubdtype(col005.dtype, np.floating)
                  and not np.isnan(col005[0]) else None)
        labels = ["dx=0.02\n(baseline)", "dx=0.01"] + (["dx=0.005"] if st005 is not None else [])
        vals = [float(row2b["st_dx002"][0]), float(row2b["st_dx001"][0])] + ([st005] if st005 is not None else [])
        colors = [C_A, C_B, C_C][:len(vals)]
        x = np.arange(len(vals))
        ax.bar(x, vals, color=colors, width=0.6, zorder=2)
        for xi, v in zip(x, vals):
            ax.text(xi, v + 0.01, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
        ax.axhline(paper_st[i], color=C_REF, linewidth=1.6, linestyle="--", zorder=3,
                   label="Kurtulus (2019)" if i == 0 else None)
        ax.text(len(vals) - 0.5, paper_st[i], f" paper: {paper_st[i]:.3f}", color=C_REF,
                fontsize=7, va="bottom", ha="right")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_title(f"α={a}°", fontsize=10)
        if i == 0:
            ax.set_ylabel("Strouhal number St")
    fig.suptitle("Test 2b: does the shedding-frequency reading change under grid refinement,\n"
                 "and does it move toward the paper's own value (dashed gray)?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(FIGS / "test_2b_dx_refine_strouhal.png", dpi=150)
    plt.close(fig)
    print("wrote test_2b_dx_refine_strouhal.png")

    fig, axes = plt.subplots(1, 4, figsize=(13, 4.3), sharey=False)
    for i, a in enumerate(angles):
        ax = axes[i]
        style_ax(ax)
        row = d2c[d2c["alpha_deg"] == a]
        vals = [float(row["st_ngrid1"][0]), float(row["st_ngrid2"][0]), float(row["st_ngrid3"][0])]
        labels = ["ngrid=1\n(baseline)", "ngrid=2", "ngrid=3"]
        colors = [C_A, C_B, C_C]
        x = np.arange(3)
        ax.bar(x, vals, color=colors, width=0.6, zorder=2)
        for xi, v in zip(x, vals):
            ax.text(xi, v + 0.01, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
        ax.axhline(paper_st[i], color=C_REF, linewidth=1.6, linestyle="--", zorder=3,
                   label="Kurtulus (2019)" if i == 0 else None)
        ax.text(2.5, paper_st[i], f" paper: {paper_st[i]:.3f}", color=C_REF,
                fontsize=7, va="bottom", ha="right")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_title(f"α={a}°", fontsize=10)
        if i == 0:
            ax.set_ylabel("Strouhal number St")
    fig.suptitle("Test 2c: does the shedding-frequency reading change with more far-field domain,\n"
                 "and does it move toward the paper's own value (dashed gray)?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(FIGS / "test_2c_ngrid_strouhal.png", dpi=150)
    plt.close(fig)
    print("wrote test_2c_ngrid_strouhal.png")


def fig_2d():
    d2d = load_csv("test2d_dxngrid_strouhal.csv")
    angles = [15, 20, 30, 40]
    paper_st = _paper_strouhal_at(angles)

    fig, axes = plt.subplots(1, 4, figsize=(13, 4.3), sharey=False)
    for i, a in enumerate(angles):
        ax = axes[i]
        style_ax(ax)
        row = d2d[d2d["alpha_deg"] == a]
        vals = [float(row["st_dx001_ngrid1"][0]), float(row["st_dx001_ngrid2"][0]),
                float(row["st_dx001_ngrid3"][0])]
        labels = ["dx=0.01\nngrid=1\n(2b baseline)", "dx=0.01\nngrid=2", "dx=0.01\nngrid=3"]
        colors = [C_A, C_B, C_C]
        x = np.arange(3)
        ax.bar(x, vals, color=colors, width=0.6, zorder=2)
        for xi, v in zip(x, vals):
            ax.text(xi, v + 0.01, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
        ax.axhline(paper_st[i], color=C_REF, linewidth=1.6, linestyle="--", zorder=3,
                   label="Kurtulus (2019)" if i == 0 else None)
        ax.text(2.5, paper_st[i], f" paper: {paper_st[i]:.3f}", color=C_REF,
                fontsize=7, va="bottom", ha="right")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_title(f"α={a}°", fontsize=10)
        if i == 0:
            ax.set_ylabel("Strouhal number St")
    fig.suptitle("Test 2d: does combining finer dx (0.01) with more far-field domain (ngrid=2,3)\n"
                 "synergize, beyond what dx=0.01 alone (2b) or ngrid alone (2c) achieved?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(FIGS / "test_2d_dxngrid_strouhal.png", dpi=150)
    plt.close(fig)
    print("wrote test_2d_dxngrid_strouhal.png")


def fig_3b():
    d = load_csv("test3b_ic_ensemble.csv")
    angles = [25, 28, 30]
    neighbors = {25: (24, 26), 28: (27, 29), 30: (29, 31)}
    ic_order = ["impulsive", "from_below", "from_above",
                "perturb_relm1em02", "perturb_relm1em03", "perturb_relp1em02"]
    colors = [C_A, C_D, C_D, C_B, C_B, C_B]

    fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=False)
    for i, a in enumerate(angles):
        ax = axes[i]
        style_ax(ax)
        below, above = neighbors[a]
        ic_labels = ["impulsive\n(baseline)", f"from below\n({below}°)", f"from above\n({above}°)",
                     "Re-1%", "Re-0.1%", "Re+1%"]
        rows = d[d["alpha_deg"] == a]
        vals = []
        for ic in ic_order:
            r = rows[rows["ic_type"] == ic]
            vals.append(float(r["cl_mean"][0]) if len(r) else np.nan)
        x = np.arange(len(ic_order))
        ax.bar(x, vals, color=colors, width=0.65, zorder=2)
        base = vals[0]
        ax.axhline(base, color="#333333", linewidth=1, linestyle=":", zorder=1)
        for xi, v in zip(x, vals):
            if not np.isnan(v):
                ax.text(xi, v + 0.01, f"{v:.2f}", ha="center", va="bottom", fontsize=7.5)
        ax.set_xticks(x)
        ax.set_xticklabels(ic_labels, fontsize=7, rotation=30, ha="right")
        ax.set_title(f"α={a}°", fontsize=10)
        if i == 0:
            ax.set_ylabel(r"mean lift coefficient $\overline{C_l}$")
    fig.suptitle("Test 3b: does the developed-state mean depend on how the run was started?\n"
                 "dashed line = impulsive-start baseline", fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(FIGS / "test_3b_ic_ensemble.png", dpi=150)
    plt.close(fig)
    print("wrote test_3b_ic_ensemble.png")


if __name__ == "__main__":
    fig_1a()
    fig_1b()
    fig_2a()
    fig_3a()
    fig_2b_2c()
    fig_2d()
    fig_3b()
