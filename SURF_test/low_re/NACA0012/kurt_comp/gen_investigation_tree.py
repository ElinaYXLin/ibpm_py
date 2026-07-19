"""
gen_investigation_tree.py

Map of the whole kurt_comp/ investigation as a left-to-right tree: root =
Kurtulus (2019)'s published results; layer 1 = every comparison chart in
1-paper_based/ that measures ibpm against them; layer 2 = the 2-follow_up/
test that investigated an anomaly found in a layer-1 chart (only charts
that actually flagged an anomaly get a layer-2 child -- charts that simply
confirmed agreement are leaves); layer 3 = the 3-further/ test that dug
into an open question left by a layer-2 test. Two layer-1 charts
(fig1_mean_coefficients.png and fig13_14_hysteresis.png) both feed the same
layer-2 node (test_DE_blockage_liftslope.png) -- the README explicitly ties
the hysteresis loop's oversized amplitude to the same ~14% lift-slope
excess, so that's a real shared cause, not a tree-drawing convenience.

Node label = "(#) file_name", # = which test-batch folder (1/2/3) the file
is under.

Usage: python3 gen_investigation_tree.py
Output: investigation_tree.png
"""
import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "investigation_tree.png")

COLOR_ROOT = "#333333"
COLOR_1 = "#c0392b"   # 1-paper_based
COLOR_2 = "#2980b9"   # 2-follow_up
COLOR_3 = "#16a085"   # 3-further
COLOR_LEAF = "#95a5a6"  # confirmed, no anomaly -- no children

COL_X = [1.5, 10.0, 20.5, 31.0]
BODY_FS = 8.6
TAG_FS = 7.6

DATA_PER_INCH = 1.5
_PT_PER_INCH = 72.0
_CHAR_W_PT_FRAC = 0.56
_LINE_H_PT_FRAC = 1.32

ALL_BOXES = []
PENDING_ARROWS = []


def _dpp():
    return DATA_PER_INCH / _PT_PER_INCH


def _text_size(lines, fontsize):
    dpp = _dpp()
    w = max((len(s) for s in lines), default=0) * fontsize * _CHAR_W_PT_FRAC * dpp
    h = len(lines) * fontsize * _LINE_H_PT_FRAC * dpp
    return w, h


def box(ax, col, y, lines, edgecolor, facecolor="white", linewidth=1.6, x=None):
    cx = COL_X[col] if x is None else x
    tw, th = _text_size(lines, BODY_FS)
    w = max(tw + 0.6, 2.2)
    h = max(th + 0.34, 0.66)
    b = FancyBboxPatch((cx - w / 2, y - h / 2), w, h,
                         boxstyle="round,pad=0.02,rounding_size=0.08",
                         linewidth=linewidth, edgecolor=edgecolor,
                         facecolor=facecolor, zorder=2)
    ax.add_patch(b)
    ax.text(cx, y, "\n".join(lines), ha="center", va="center", fontsize=BODY_FS,
             color="#222222", zorder=3, linespacing=_LINE_H_PT_FRAC)
    rec = (cx, y, w, h)
    ALL_BOXES.append(rec)
    return rec


def right(b):
    x, y, w, h = b
    return (x + w / 2, y)


def left(b):
    x, y, w, h = b
    return (x - w / 2, y)


def arrow(ax, src, dst, color, lw=1.4):
    a = FancyArrowPatch(right(src), left(dst), arrowstyle="-|>", mutation_scale=13,
                         color=color, linewidth=lw, zorder=1,
                         connectionstyle="arc3,rad=0.05")
    ax.add_patch(a)


def main():
    fig, ax = plt.subplots(figsize=(10, 10))

    root = box(ax, 0, 20, ["Kurtulus (2019)'s", "Published Results"], COLOR_ROOT,
               facecolor="#f2f2f2")

    # ---------------- layer 1: every 1-paper_based comparison chart ----------------
    l1_y = {}
    specs1 = [
        ("fig1_mean_coefficients.png", 34),
        ("fig11_instantaneous_pitch.png", 30),
        ("fig13_14_hysteresis.png", 26.5),
        ("fig19_shedding_strouhal.png", 22),
        ("thrust_check.png", 17),
        ("wake_steady.png", 12),
        ("wake_f4hz.png", 8),
    ]
    for name, y in specs1:
        l1_y[name] = box(ax, 1, y, [f"(1) {name}"], COLOR_1)
        arrow(ax, root, l1_y[name], COLOR_ROOT)

    # ---------------- layer 2: 2-follow_up tests, only for flagged anomalies ----------------
    l2_y = {}
    specs2 = [
        ("test_DE_blockage_liftslope.png", 32, ["fig1_mean_coefficients.png",
                                                  "fig13_14_hysteresis.png"]),
        ("test_FG_alpha0_offset.png", 27, ["fig1_mean_coefficients.png"]),
        ("test_B_period_locked_mean.png", 22.5, ["fig1_mean_coefficients.png"]),
        ("test_C_strouhal_resolution.png", 17.5, ["fig19_shedding_strouhal.png"]),
        ("test_A_thrust_window.png", 12, ["thrust_check.png"]),
    ]
    for name, y, parents in specs2:
        l2_y[name] = box(ax, 2, y, [f"(2) {name}"], COLOR_2)
        for p in parents:
            arrow(ax, l1_y[p], l2_y[name], COLOR_1)

    # leaves at layer 1 with no anomaly to chase (drawn as slightly muted)
    for name in ("fig11_instantaneous_pitch.png", "wake_steady.png", "wake_f4hz.png"):
        b = l1_y[name]
        ax.text(right(b)[0] + 0.3, b[1], "confirmed --\nno anomaly", ha="left", va="center",
                fontsize=7, color=COLOR_LEAF, style="italic")

    # ---------------- layer 3: 3-further tests, digging into open questions ----------------
    specs3 = [
        ("test_3a_running_mean_convergence.png", 25, "test_B_period_locked_mean.png"),
        ("test_3b_ic_ensemble.png", 21, "test_B_period_locked_mean.png"),
        ("test_2a_spectrogram.png", 19.5, "test_C_strouhal_resolution.png"),
        ("test_2b_dx_refine_strouhal.png", 16.5, "test_C_strouhal_resolution.png"),
        ("test_2c_ngrid_strouhal.png", 13.5, "test_C_strouhal_resolution.png"),
        ("test_1a_full_phase_average.png", 10.5, "test_A_thrust_window.png"),
        ("test_1b_dip_vs_transitions.png", 7.5, "test_A_thrust_window.png"),
    ]
    for name, y, parent in specs3:
        b3 = box(ax, 3, y, [f"(3) {name}"], COLOR_3)
        arrow(ax, l2_y[parent], b3, COLOR_2)

    # no further branch for test_DE / test_FG -- both fully explained their anomaly
    for name in ("test_DE_blockage_liftslope.png", "test_FG_alpha0_offset.png"):
        b = l2_y[name]
        ax.text(right(b)[0] + 0.3, b[1], "resolved --\nno open question", ha="left", va="center",
                fontsize=7, color=COLOR_LEAF, style="italic")

    # ---------------- legend ----------------
    leg_x, leg_y = COL_X[3] + 6.5, 34
    ax.text(leg_x, leg_y, "Legend", fontsize=10, fontweight="bold")
    items = [
        (COLOR_ROOT, "Kurtulus (2019) -- the published paper"),
        (COLOR_1, "(1) 1-paper_based/figures/ -- ibpm vs. the paper"),
        (COLOR_2, "(2) 2-follow_up/figures/ -- is an anomaly real or an artifact?"),
        (COLOR_3, "(3) 3-further/figures/ -- digging into what's left open"),
        (COLOR_LEAF, "italic note = investigation ended here"),
    ]
    yy = leg_y - 1.1
    for color, label in items:
        ax.add_patch(FancyBboxPatch((leg_x, yy - 0.22), 0.5, 0.34, boxstyle="round,pad=0.02",
                                      linewidth=1.4, edgecolor=color, facecolor="white", zorder=2))
        ax.text(leg_x + 0.7, yy - 0.05, label, fontsize=7.8, va="center")
        yy -= 0.9

    ax.set_title(
        "kurt_comp/ investigation map: from Kurtulus (2019)'s published results,\n"
        "through every ibpm comparison, to the follow-up tests each flagged anomaly spawned",
        fontsize=13, pad=10)

    xs = [b[0] - b[2] / 2 for b in ALL_BOXES] + [b[0] + b[2] / 2 for b in ALL_BOXES] + [leg_x + 6]
    ys = [b[1] - b[3] / 2 for b in ALL_BOXES] + [b[1] + b[3] / 2 for b in ALL_BOXES]
    xmin, xmax = min(xs) - 0.5, max(xs) + 0.5
    ymin, ymax = min(ys) - 1.0, max(ys) + 1.5
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")

    title_in = 0.6
    figw = (xmax - xmin) / DATA_PER_INCH
    figh = (ymax - ymin) / DATA_PER_INCH + title_in
    fig.set_size_inches(figw, figh)
    fig.subplots_adjust(left=0.005, right=0.995, top=1 - title_in / figh, bottom=0.005)

    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT} ({figw:.1f}in x {figh:.1f}in, {len(ALL_BOXES)} boxes)")


if __name__ == "__main__":
    main()
