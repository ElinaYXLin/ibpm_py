"""
gen_followup_figs_hm.py

Figures for Tests H-M (see README.md's "H-M" section and
analyze_followup_hm.py). Run analyze_followup_hm.py first.

Usage: python3 gen_followup_figs_hm.py
Output: figures/test_H_*.png, figures/test_IJKL_*.png, figures/test_M_*.png
"""
import csv
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures"
FIGS.mkdir(exist_ok=True)

C_A = "#c0392b"
C_B = "#2980b9"
C_REF = "#7f8c8d"


def load_csv(name):
    with open(DATA / name) as f:
        return list(csv.DictReader(f))


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def fig_h():
    rows = load_csv("test_H_steady_vs_dynamic.csv")
    a = [int(r["alpha_deg"]) for r in rows]
    ibpm_steady = [float(r["ibpm_cd_steady"]) for r in rows]
    ibpm_f4hz = [float(r["ibpm_cd_f4hz"]) for r in rows]
    paper = [float(r["paper_cd_steady"]) for r in rows]  # steady==f4hz in the digitized data at these alpha

    fig, ax = plt.subplots(figsize=(8, 5))
    style_ax(ax)
    ax.fill_between(a, [p - 0.05 for p in paper], [p + 0.05 for p in paper],
                     color=C_REF, alpha=0.18, label="paper +/-0.05 digitization band")
    ax.plot(a, paper, "--", color=C_REF, lw=1.6, marker="s", ms=5, label="Kurtulus (2019) Fig 1")
    ax.plot(a, ibpm_steady, "-", color=C_A, lw=1.8, marker="o", ms=5, label="ibpm steady mean $C_d$")
    ax.plot(a, ibpm_f4hz, "-", color=C_B, lw=1.8, marker="o", ms=5, label="ibpm f4Hz mean $C_d$")
    ax.set_xlabel(r"mean angle of attack $\alpha_0$ (deg)")
    ax.set_ylabel(r"mean $C_d$")
    ax.set_title("Test H: is the Cd excess already in the steady baseline,\n"
                  "or specific to the pitching oscillation?", fontsize=11)
    ax.legend(fontsize=8.5, loc="upper left")
    fig.tight_layout()
    out = FIGS / "test_H_steady_vs_dynamic.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_ijkl():
    p = DATA / "test_IJKL_knob_sensitivity.csv"
    if not p.exists():
        print("test_IJKL: no data yet"); return
    rows = [r for r in load_csv(p.name) if r["impl"] == "py"]  # cpp is identical, see data/ CSV
    labels = [r["case"] for r in rows]
    cl = [float(r["cl_pk2pk_ratio"]) for r in rows]
    cd = [float(r["cd_pk2pk_ratio"]) for r in rows]

    fig, ax = plt.subplots(figsize=(1.6 * len(labels) + 2, 5.2))
    style_ax(ax)
    x = np.arange(len(labels))
    ax.axhline(1.0, color="#333333", lw=1.3, ls="--", zorder=1, label="perfect match (ibpm=paper)")
    ax.bar(x - 0.18, cl, width=0.36, color=C_A, label="$C_l$ peak-to-peak ratio", zorder=2)
    ax.bar(x + 0.18, cd, width=0.36, color=C_B, label="$C_d$ peak-to-peak ratio", zorder=2)
    for xi, v in zip(x, cl):
        ax.text(xi - 0.18, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=7.5)
    for xi, v in zip(x, cd):
        ax.text(xi + 0.18, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("ibpm / paper peak-to-peak ratio")
    ax.set_title("Tests I-L: does any single knob shrink Cd's excess toward the paper "
                  "(ratio -> 1.0),\nwithout also changing Cl's? (py_static shown; cpp_static "
                  "identical in every case, see data/)", fontsize=10.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = FIGS / "test_IJKL_knob_sensitivity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out.name}")


def fig_m():
    rows = load_csv("test_M_paper_selfcheck.csv")
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    style_ax(ax)
    x = np.arange(len(rows))
    table_vals = [float(r["fig13_14_table_mean"]) for r in rows]
    fig1_vals = [float(r["fig1_digitized_value"]) for r in rows]
    labels = [r["quantity"] for r in rows]
    ax.bar(x - 0.18, table_vals, width=0.36, color=C_A, label="Fig 13/14 table\n(time-weighted mean)")
    ax.bar(x + 0.18, fig1_vals, width=0.36, color=C_REF, label="Fig 1 (digitized,\n+/-0.05)")
    for xi, v in zip(x, table_vals):
        ax.text(xi - 0.18, v + 0.02, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    for xi, v in zip(x, fig1_vals):
        ax.text(xi + 0.18, v + 0.02, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("value at alpha0=0, f=4Hz")
    ax.set_title("Test M: does the paper's own Fig 13/14 table agree with\nits own Fig 1, at the same condition?",
                  fontsize=10.5)
    ax.legend(fontsize=8)
    ax.text(1, 0.35, "Cl comparison isn't apples-to-apples:\nthe table only covers ~72% of one\n"
            "cycle (down+up branches), so its\ntime-average needn't match Fig 1's\ntrue full-cycle mean of 0.\n"
            "Not evidence of an inconsistency.", fontsize=7.5, color="#555555", ha="center", va="center",
            style="italic")
    fig.tight_layout()
    out = FIGS / "test_M_paper_selfcheck.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    fig_h()
    fig_ijkl()
    fig_m()
