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

Layout note: column x-positions and row y-positions are computed from the
actual rendered box sizes (which scale with BODY_FS) rather than hardcoded,
so cranking the font up packs everything tighter instead of leaving the
boxes marooned in a hardcoded, small-font-era grid.

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

# 1 data unit == 1 inch (kept explicit so every spacing constant below can be
# read directly as "this many inches"), which is what makes it possible to
# size the whole layout off of BODY_FS instead of hand-picked coordinates.
DATA_PER_INCH = 1.0
_PT_PER_INCH = 72.0
_CHAR_W_PT_FRAC = 0.56
_LINE_H_PT_FRAC = 1.32

BODY_FS = 17          # node label font size (was 8.6)
LEAF_FS = 13.5        # "confirmed/resolved" side-note font size (was 7/7.8)
TITLE_FS = 21         # figure title (was 13)
LEGEND_HEAD_FS = 16   # "Legend" header (was 10)
LEGEND_FS = 13.5      # legend item labels (was 7.8)

ROW_GAP = 0.32        # vertical gap between stacked boxes in the same column
COL_GAP_PLAIN = 1.0   # column gap where boxes have no right-side leaf note
COL_GAP_ANNOT = 2.6   # column gap where boxes may carry a right-side leaf note
LEFT_MARGIN = 0.3

ALL_BOXES = []


def _dpp():
    return DATA_PER_INCH / _PT_PER_INCH


def _text_size(lines, fontsize):
    dpp = _dpp()
    w = max((len(s) for s in lines), default=0) * fontsize * _CHAR_W_PT_FRAC * dpp
    h = len(lines) * fontsize * _LINE_H_PT_FRAC * dpp
    return w, h


def box_dims(lines, fontsize=BODY_FS):
    tw, th = _text_size(lines, fontsize)
    w = max(tw + 0.6, 2.2)
    h = max(th + 0.34, 0.66)
    return w, h


def box(ax, x, y, lines, edgecolor, facecolor="white", linewidth=1.6):
    w, h = box_dims(lines)
    b = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                         boxstyle="round,pad=0.02,rounding_size=0.08",
                         linewidth=linewidth, edgecolor=edgecolor,
                         facecolor=facecolor, zorder=2)
    ax.add_patch(b)
    ax.text(x, y, "\n".join(lines), ha="center", va="center", fontsize=BODY_FS,
             color="#222222", zorder=3, linespacing=_LINE_H_PT_FRAC)
    rec = (x, y, w, h)
    ALL_BOXES.append(rec)
    return rec


def right(b):
    x, y, w, h = b
    return (x + w / 2, y)


def left(b):
    x, y, w, h = b
    return (x - w / 2, y)


def arrow(ax, src, dst, color, lw=1.4):
    a = FancyArrowPatch(right(src), left(dst), arrowstyle="-|>", mutation_scale=15,
                         color=color, linewidth=lw, zorder=1,
                         connectionstyle="arc3,rad=0.05")
    ax.add_patch(a)


def stacked_ys(n, h, gap, center=0.0):
    """n boxes of height h, top-to-bottom, centered on `center`."""
    total = n * h + (n - 1) * gap
    y = center + total / 2 - h / 2
    ys = []
    for _ in range(n):
        ys.append(y)
        y -= h + gap
    return ys


def main():
    fig, ax = plt.subplots()

    root_lines = ["Kurtulus (2019)'s", "Published Results"]

    l1_names = [
        "fig1_mean_coefficients.png",
        "fig11_instantaneous_pitch.png",
        "fig13_14_hysteresis.png",
        "fig19_shedding_strouhal.png",
        "thrust_check.png",
        "wake_steady.png",
        "wake_f4hz.png",
    ]
    leaf1 = {
        "fig11_instantaneous_pitch.png": "confirmed --\nno anomaly",
        "wake_steady.png": "confirmed --\nno anomaly",
        "wake_f4hz.png": "confirmed --\nno anomaly",
    }

    l2_defs = [
        ("test_DE_blockage_liftslope.png", ["fig1_mean_coefficients.png",
                                             "fig13_14_hysteresis.png"]),
        ("test_FG_alpha0_offset.png", ["fig1_mean_coefficients.png"]),
        ("test_B_period_locked_mean.png", ["fig1_mean_coefficients.png"]),
        ("test_C_strouhal_resolution.png", ["fig19_shedding_strouhal.png"]),
        ("test_A_thrust_window.png", ["thrust_check.png"]),
    ]
    leaf2 = {
        "test_DE_blockage_liftslope.png": "resolved --\nno open question",
        "test_FG_alpha0_offset.png": "resolved --\nno open question",
    }

    l3_defs = [
        ("test_3a_running_mean_convergence.png", "test_B_period_locked_mean.png"),
        ("test_3b_ic_ensemble.png", "test_B_period_locked_mean.png"),
        ("test_2a_spectrogram.png", "test_C_strouhal_resolution.png"),
        ("test_2b_dx_refine_strouhal.png", "test_C_strouhal_resolution.png"),
        ("test_2c_ngrid_strouhal.png", "test_C_strouhal_resolution.png"),
        ("test_1a_full_phase_average.png", "test_A_thrust_window.png"),
        ("test_1b_dip_vs_transitions.png", "test_A_thrust_window.png"),
    ]

    # ---------------- sizes (drive column x-positions off actual box sizes) ----------------
    root_w, root_h = box_dims(root_lines)
    l1_sizes = {n: box_dims([f"(1) {n}"]) for n in l1_names}
    l2_sizes = {n: box_dims([f"(2) {n}"]) for n, _ in l2_defs}
    l3_sizes = {n: box_dims([f"(3) {n}"]) for n, _ in l3_defs}
    l1_w_max = max(w for w, h in l1_sizes.values())
    l2_w_max = max(w for w, h in l2_sizes.values())
    l3_w_max = max(w for w, h in l3_sizes.values())
    l1_h = next(iter(l1_sizes.values()))[1]
    l2_h = next(iter(l2_sizes.values()))[1]
    l3_h = next(iter(l3_sizes.values()))[1]

    col_x0 = LEFT_MARGIN + root_w / 2
    col_x1 = col_x0 + root_w / 2 + COL_GAP_PLAIN + l1_w_max / 2
    col_x2 = col_x1 + l1_w_max / 2 + COL_GAP_ANNOT + l2_w_max / 2
    col_x3 = col_x2 + l2_w_max / 2 + COL_GAP_ANNOT + l3_w_max / 2

    l1_ys = stacked_ys(len(l1_names), l1_h, ROW_GAP)
    l2_ys = stacked_ys(len(l2_defs), l2_h, ROW_GAP)
    l3_ys = stacked_ys(len(l3_defs), l3_h, ROW_GAP)

    # ---------------- draw ----------------
    root = box(ax, col_x0, 0.0, root_lines, COLOR_ROOT, facecolor="#f2f2f2")

    l1_y = {}
    for name, y in zip(l1_names, l1_ys):
        l1_y[name] = box(ax, col_x1, y, [f"(1) {name}"], COLOR_1)
        arrow(ax, root, l1_y[name], COLOR_ROOT)
    for name, note in leaf1.items():
        b = l1_y[name]
        ax.text(right(b)[0] + 0.3, b[1], note, ha="left", va="center",
                fontsize=LEAF_FS, color=COLOR_LEAF, style="italic")

    l2_y = {}
    for (name, parents), y in zip(l2_defs, l2_ys):
        l2_y[name] = box(ax, col_x2, y, [f"(2) {name}"], COLOR_2)
        for p in parents:
            arrow(ax, l1_y[p], l2_y[name], COLOR_1)
    for name, note in leaf2.items():
        b = l2_y[name]
        ax.text(right(b)[0] + 0.3, b[1], note, ha="left", va="center",
                fontsize=LEAF_FS, color=COLOR_LEAF, style="italic")

    for (name, parent), y in zip(l3_defs, l3_ys):
        b3 = box(ax, col_x3, y, [f"(3) {name}"], COLOR_3)
        arrow(ax, l2_y[parent], b3, COLOR_2)

    # ---------------- legend ----------------
    col_top = max(root_h / 2,
                  (l1_h * len(l1_names) + ROW_GAP * (len(l1_names) - 1)) / 2,
                  (l2_h * len(l2_defs) + ROW_GAP * (len(l2_defs) - 1)) / 2,
                  (l3_h * len(l3_defs) + ROW_GAP * (len(l3_defs) - 1)) / 2)
    leg_x = col_x3 + l3_w_max / 2 + 1.3
    leg_y = col_top
    ax.text(leg_x, leg_y, "Legend", fontsize=LEGEND_HEAD_FS, fontweight="bold")
    items = [
        (COLOR_ROOT, "Kurtulus (2019) -- the published paper"),
        (COLOR_1, "(1) 1-paper_based/figures/ -- ibpm vs. the paper"),
        (COLOR_2, "(2) 2-follow_up/figures/ -- is an anomaly real or an artifact?"),
        (COLOR_3, "(3) 3-further/figures/ -- digging into what's left open"),
        (COLOR_LEAF, "italic note = investigation ended here"),
    ]
    yy = leg_y - 0.85
    leg_label_w = max(len(label) for _, label in items) * LEGEND_FS * _CHAR_W_PT_FRAC / _PT_PER_INCH
    for color, label in items:
        ax.add_patch(FancyBboxPatch((leg_x, yy - 0.2), 0.55, 0.36, boxstyle="round,pad=0.02",
                                      linewidth=1.4, edgecolor=color, facecolor="white", zorder=2))
        ax.text(leg_x + 0.75, yy - 0.02, label, fontsize=LEGEND_FS, va="center")
        yy -= 0.75

    ax.set_title(
        "kurt_comp/ investigation map: from Kurtulus (2019)'s published results,\n"
        "through every ibpm comparison, to the follow-up tests each flagged anomaly spawned",
        fontsize=TITLE_FS, pad=10)

    xs = [b[0] - b[2] / 2 for b in ALL_BOXES] + [b[0] + b[2] / 2 for b in ALL_BOXES]
    xs += [leg_x + 0.75 + leg_label_w]
    ys = [b[1] - b[3] / 2 for b in ALL_BOXES] + [b[1] + b[3] / 2 for b in ALL_BOXES]
    ys += [leg_y + 0.5, yy + 0.5]
    xmin, xmax = min(xs) - 0.3, max(xs) + 0.3
    ymin, ymax = min(ys) - 0.4, max(ys) + 0.4
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")

    title_in = 1.15
    figw = (xmax - xmin) / DATA_PER_INCH
    figh = (ymax - ymin) / DATA_PER_INCH + title_in
    fig.set_size_inches(figw, figh)
    fig.subplots_adjust(left=0.005, right=0.995, top=1 - title_in / figh, bottom=0.005)

    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT} ({figw:.1f}in x {figh:.1f}in, {len(ALL_BOXES)} boxes)")


if __name__ == "__main__":
    main()
