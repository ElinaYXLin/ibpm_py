"""
generate_execution_flowchart_v2.py

"More user-friendly" simplified version of main_execution_flowchart.png,
built from a draft outline (flowchart_draft.docx) rather than an exhaustive
call-graph trace. Same underlying run (py/ibpm.py's main()), but pruned to
~29 boxes instead of V1's ~60+, and only 3 layers of call depth (vs. V1's up
to 8). Depth is conveyed purely by column position and rightward arrows (no
number/letter/roman-numeral labels on the boxes themselves): column 0 is a
top-level step (source order in ibpm.py's main()), column 1 is one call
deeper, column 2 is one call deeper still.

Box text is the draft's own wording verbatim, EXCEPT boxes that were marked
"(unsure)" in the draft -- those were corrected against the actual source
(the correction, and why, is in flowchart/output_V2/README.md). Every box
still carries a file:line tag, positioned close above the box it cites,
matching the citation style of V1's diagram (see
generate_execution_flowchart.py) -- with enough vertical clearance from the
box above it in the same column that the two never overlap.

Usage:
    python3 flowchart/generate_execution_flowchart_v2.py

Output:
    flowchart/output_V2/main_execution_flowchart_v2.png
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output_V2")

COLOR_DRIVER = "#c0392b"      # ibpm.py:main()
COLOR_SOLVER = "#8e44ad"      # ib_solver.py
COLOR_MODEL = "#16a085"       # navier_stokes_model.py
COLOR_IO = "#546e7a"          # logger.py

LEGEND = [
    (COLOR_DRIVER, "ibpm.py (driver / main())"),
    (COLOR_MODEL, "navier_stokes_model.py (NavierStokesModel)"),
    (COLOR_SOLVER, "ib_solver.py (IBSolver + subclasses)"),
    (COLOR_IO, "logger.py"),
]

FILE_IBPM = "ibpm.py"
FILE_IBS = "ib_solver.py"
FILE_NSM = "navier_stokes_model.py"
FILE_LOGGER = "logger.py"

COL_X = [3.0, 13.5, 24.0]

BODY_FS = 8.4
TAG_FS = 6.8
TAG_GAP = 0.24         # vertical gap between a box's top edge and its file:line tag
                        # (large enough that an incoming arrowhead landing at
                        # top-center never pokes into the tag text on narrow boxes)

DATA_PER_INCH = 1.55
_PT_PER_INCH = 72.0
_CHAR_W_PT_FRAC = 0.56
_LINE_H_PT_FRAC = 1.32

ALL_BOXES = []
PENDING_ARROWS = []


def _data_per_pt():
    return DATA_PER_INCH / _PT_PER_INCH


def _text_size_data(lines, fontsize):
    dpp = _data_per_pt()
    max_chars = max((len(s) for s in lines), default=0)
    w = max_chars * fontsize * _CHAR_W_PT_FRAC * dpp
    h = len(lines) * fontsize * _LINE_H_PT_FRAC * dpp
    return w, h


def box_height(lines):
    _, th = _text_size_data(lines, BODY_FS)
    return max(th + 0.34, 0.66)


def tag_height():
    return TAG_FS * _LINE_H_PT_FRAC * _data_per_pt()


def box(ax, col, y, lines, edgecolor, file, tag_lines, facecolor="white",
        linewidth=1.5, zorder=2, x=None):
    """One flowchart box. `file:tag_lines` is drawn just above the box's
    top-left corner (small gap, TAG_GAP), citing exactly where the code this
    box shows is written (same convention as V1) -- no number/letter/roman
    label; depth is conveyed by column position and the rightward arrows
    alone. Callers are responsible for leaving enough vertical room in the
    same column that this tag doesn't collide with the box placed above it
    (see the dy spacing used in main(), sized for this)."""
    cx = COL_X[col] if x is None else x
    tw, th = _text_size_data(lines, BODY_FS)
    w = max(tw + 0.55, 2.0)
    h = max(th + 0.34, 0.66)
    b = FancyBboxPatch((cx - w / 2, y - h / 2), w, h,
                         boxstyle="round,pad=0.02,rounding_size=0.08",
                         linewidth=linewidth, edgecolor=edgecolor,
                         facecolor=facecolor, zorder=zorder)
    ax.add_patch(b)
    ax.text(cx, y, "\n".join(lines), ha="center", va="center", fontsize=BODY_FS,
             color="#222222", zorder=zorder + 1, linespacing=_LINE_H_PT_FRAC)
    tag = f"{file}:{tag_lines}"
    ax.text(cx - w / 2 + 0.08, y + h / 2 + TAG_GAP, tag, ha="left", va="bottom",
             fontsize=TAG_FS, style="italic", color=edgecolor, zorder=zorder + 1)
    rec = (cx, y, w, h)
    ALL_BOXES.append(rec)
    return rec


def section_header(ax, y, text):
    ax.text(COL_X[0] - 1.6, y, text, ha="left", va="center", fontsize=11.5,
             fontweight="bold", color="#111111")
    ax.plot([COL_X[0] - 1.6, COL_X[2] + 3.2], [y - 0.42, y - 0.42],
             color="#bbbbbb", linewidth=1.0, zorder=0)
    return y


def top(b):
    x, y, w, h = b
    return (x, y + h / 2)


def bottom(b):
    x, y, w, h = b
    return (x, y - h / 2)


def left(b):
    x, y, w, h = b
    return (x - w / 2, y)


def right(b):
    x, y, w, h = b
    return (x + w / 2, y)


def call_arrow(ax, src, dst, color="#333333", lw=1.3, mutation_scale=12):
    assert dst[0] > src[0]
    PENDING_ARROWS.append(dict(kind="call", src=src, dst=dst, color=color, lw=lw,
                                mutation_scale=mutation_scale))


def seq_arrow(ax, src, dst, color, lw=1.4, ls="solid"):
    assert dst[0] == src[0]
    PENDING_ARROWS.append(dict(kind="seq", src=src, dst=dst, color=color, lw=lw, ls=ls))


def loop_arrow(ax, src, dst, color, x_offset, label, fontsize=7.5):
    lx = src[0] + x_offset
    a1 = FancyArrowPatch(right(src) if x_offset > 0 else left(src),
                          (lx, src[1]), arrowstyle="-", mutation_scale=1,
                          color=color, linewidth=1.2, zorder=1)
    a2 = FancyArrowPatch((lx, src[1]), (lx, dst[1]), arrowstyle="-|>",
                          mutation_scale=12, color=color, linewidth=1.2,
                          linestyle="dashed", zorder=1)
    a3 = FancyArrowPatch((lx, dst[1]), right(dst) if x_offset > 0 else left(dst),
                          arrowstyle="-", mutation_scale=1, color=color,
                          linewidth=1.2, zorder=1)
    for a in (a1, a2, a3):
        ax.add_patch(a)
    ax.text(lx + (0.4 if x_offset > 0 else -0.4), (src[1] + dst[1]) / 2, label,
             ha="center", va="center", fontsize=fontsize, color=color, rotation=90)


def flush_arrows(ax):
    for spec in PENDING_ARROWS:
        src, dst = spec["src"], spec["dst"]
        if spec["kind"] == "call":
            p1, p2 = right(src), left(dst)
        else:
            p1, p2 = bottom(src), top(dst)
        a = FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=spec.get("mutation_scale", 12),
                             color=spec["color"], linewidth=spec["lw"],
                             linestyle=spec.get("ls", "solid"), zorder=1)
        ax.add_patch(a)


def draw_legend(ax, x, y_top):
    row_h = 0.75
    swatch_w, swatch_h = 0.5, 0.32
    ax.text(x, y_top, "Legend: box color = source file", ha="left", va="bottom",
             fontsize=9.5, fontweight="bold", color="#222222")
    y = y_top - 0.6
    for color, label in LEGEND:
        ax.add_patch(Rectangle((x, y - swatch_h / 2), swatch_w, swatch_h,
                                 facecolor="white", edgecolor=color, linewidth=1.8, zorder=2))
        ax.text(x + swatch_w + 0.2, y, label, ha="left", va="center", fontsize=8.2,
                 color="#222222", zorder=2)
        y -= row_h
    ax.text(x, y - 0.25,
            "Column 0 = a top-level step (source order in ibpm.py's main()).\n"
            "Column 1 = one layer of function-call depth to the right.\n"
            "Column 2 = one more layer of call depth to the right.\n"
            "Every box cites file:line for the code IT shows, just above it.\n"
            "Dashed arrow = loop-back (repeat).\n\n"
            "Box text = the draft's own wording, verbatim, EXCEPT boxes\n"
            "originally marked \"(unsure)\" -- those were corrected against\n"
            "the source code; see output_V2/README.md for what changed and why,\n"
            "and for other draft inaccuracies that were fixed here too.",
            ha="left", va="top", fontsize=7.6, color="#444444", style="italic")


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 10))

    y = 46.0
    dy = 2.0

    # ============================== SECTION 1: preparation ==============================
    section_header(ax, y, "SECTION 1: preparation")
    y -= 1.4

    b1 = box(ax, 0, y, ["Read arguments"], COLOR_DRIVER, FILE_IBPM, "166")
    y -= dy
    b2 = box(ax, 0, y, ["Build grid"], COLOR_DRIVER, FILE_IBPM, "280")
    y -= dy
    b3 = box(ax, 0, y, ["Load .geom file"], COLOR_DRIVER, FILE_IBPM, "285")
    y -= dy
    b4 = box(ax, 0, y, ["Create the free-stream", "(\"base\") flow, at the", "given angle of attack"],
             COLOR_DRIVER, FILE_IBPM, "298")
    y -= dy
    b5 = box(ax, 0, y, ["Build model and solver"], COLOR_DRIVER, FILE_IBPM, "318-347")
    b5a = box(ax, 1, y + 0.85, ["Build a Navier-Stokes model"], COLOR_MODEL, FILE_NSM, "37-64")
    call_arrow(ax, b5, b5a, color=COLOR_MODEL)
    b5b = box(ax, 1, y - 0.85,
              ["Build a solver, like Nonlinear,", "Linearized, Adjoint, Periodic, or SFD"],
              COLOR_SOLVER, FILE_IBS, "84-104,230-408")
    call_arrow(ax, b5, b5b, color=COLOR_SOLVER)

    # ============================== SECTION 2: initialization ==============================
    y -= dy + 1.5
    section_header(ax, y, "SECTION 2: initialization")
    y -= 1.4

    b6 = box(ax, 0, y, ["Initialize state and load", "initial conditions"], COLOR_DRIVER,
             FILE_IBPM, "361-369")
    y -= dy
    b7 = box(ax, 0, y, ["Move bodies to their", "position at the current time"], COLOR_DRIVER,
             FILE_IBPM, "394")
    y -= dy
    b8 = box(ax, 0, y, ["Initialize the model"], COLOR_DRIVER, FILE_IBPM, "397")
    y -= dy
    b9 = box(ax, 0, y, ["Load solver if available,", "else initialize and save it"],
             COLOR_DRIVER, FILE_IBPM, "400-405")
    y -= dy
    b10 = box(ax, 0, y, ["Update operators: re-sync the", "base flow / body positions /",
              "regularizer to the current time"], COLOR_DRIVER, FILE_IBPM, "411")
    y -= dy
    b11 = box(ax, 0, y, ["Refresh state in the", "Navier-Stokes model"], COLOR_DRIVER,
              FILE_IBPM, "412")
    b11a = box(ax, 1, y, ["Compute the flux"], COLOR_MODEL, FILE_NSM, "158-166")
    call_arrow(ax, b11, b11a, color=COLOR_MODEL)

    # ============================== SECTION 3: first round of outputs ==============================
    y -= dy + 1.5
    section_header(ax, y, "SECTION 3: first round of outputs")
    y -= 1.4

    b12 = box(ax, 0, y, ["Set up output objects", "(Tecplot/Restart/Force/Energy)", "and add to logger"],
              COLOR_DRIVER, FILE_IBPM, "417-443")
    y -= dy
    b13 = box(ax, 0, y, ["Initialize the logger"], COLOR_DRIVER, FILE_IBPM, "445")
    y -= dy
    b14 = box(ax, 0, y, ["Perform initial output"], COLOR_DRIVER, FILE_IBPM, "446")
    b14a = box(ax, 1, y, ["Call each entry that needs", "to be outputted, and do the output"],
               COLOR_IO, FILE_LOGGER, "66-70")
    call_arrow(ax, b14, b14a, color=COLOR_IO)

    # ============================== SECTION 4: solver loop ==============================
    y -= dy + 1.5
    section_header(ax, y, "SECTION 4: solver loop, for every step up to numSteps")
    y -= 1.4

    b_loophead = box(ax, 0, y, ["For each step,", "1 to numSteps:"], COLOR_DRIVER, FILE_IBPM,
                      "449", facecolor="#fdf2ef")
    y -= dy

    b15 = box(ax, 0, y, ["Advance the solver"], COLOR_DRIVER, FILE_IBPM, "452",
              facecolor="#f5eefb")
    b15a = box(ax, 1, y + 0.85, ["For every substep, set the", "nonlinear term based on the solver"],
               COLOR_SOLVER, FILE_IBS, "175,177")
    call_arrow(ax, b15, b15a, color=COLOR_SOLVER)
    b15b = box(ax, 1, y - 1.25, ["For each substep,", "do the following:"], COLOR_SOLVER,
               FILE_IBS, "185-216")
    call_arrow(ax, b15, b15b, color=COLOR_SOLVER)

    ry = y + 1.6
    romans = [
        (["Update operators, if the model", "is time dependent"], "187-188"),
        (["Calculate the \"a\" and \"b\"", "right-hand-side terms for the", "constrained (projection) solve"],
         "191-203,206"),
        (["Solve using the projection", "solver, given a and b"], "209"),
        (["Refresh the state"], "212"),
        (["If SFD solver: also integrate", "the filtered state ωhat"], "354-386"),
    ]
    roman_boxes = []
    for lines_, ln in romans:
        rb = box(ax, 2, ry, lines_, COLOR_SOLVER, FILE_IBS, ln)
        call_arrow(ax, b15b, rb, color=COLOR_SOLVER, mutation_scale=10)
        roman_boxes.append(rb)
        ry -= 1.55

    branch_bottom = roman_boxes[-1][1] - roman_boxes[-1][3] / 2

    y = min(branch_bottom, b15b[1] - b15b[3] / 2) - dy
    b16 = box(ax, 0, y, ["Compute net force"], COLOR_DRIVER, FILE_IBPM, "453")
    y -= dy
    # Split from a single "output + cleanup" box (flagged as an error in the
    # draft): logger.doOutput (462) runs EVERY step, inside the loop -- it's
    # what the loop-back arrow wraps around, along with b_loophead/b15/b16.
    # logger.cleanup (480) runs ONCE, after the loop is entirely done, so
    # it's drawn OUTSIDE the loop-back and not part of the repeated cycle.
    b17 = box(ax, 0, y, ["Output results for this step"], COLOR_DRIVER,
              FILE_IBPM, "462")

    seq_arrow(ax, b_loophead, b15, COLOR_DRIVER)
    seq_arrow(ax, b15, b16, COLOR_DRIVER)
    seq_arrow(ax, b16, b17, COLOR_DRIVER)
    loop_arrow(ax, b17, b_loophead, COLOR_DRIVER, x_offset=-2.2,
               label="repeat until step = numSteps")

    y -= dy + 1.3
    b18 = box(ax, 0, y, ["Clean up the logger", "(once, after the loop)"], COLOR_DRIVER,
              FILE_IBPM, "480")
    seq_arrow(ax, b17, b18, COLOR_DRIVER)

    # ============================== column-0 sequencing (top to bottom, section by section) ==============================
    for p, q in [(b1, b2), (b2, b3), (b3, b4), (b4, b5), (b5, b6), (b6, b7), (b7, b8),
                 (b8, b9), (b9, b10), (b10, b11), (b11, b12), (b12, b13), (b13, b14),
                 (b14, b_loophead)]:
        seq_arrow(ax, p, q, COLOR_DRIVER)

    draw_legend(ax, COL_X[2] + 1.6, b1[1] + 0.3)

    flush_arrows(ax)

    ax.set_title(
        "ibpm_py -- execution flow (simplified / \"user-friendly\" version, V2)\n"
        "Based on a draft outline, cross-checked against the source. Column 0 = a top-level "
        "step; column 1 = one call-layer right; column 2 = one more call-layer right. See "
        "legend for what was corrected vs. the original draft.",
        fontsize=12.5, pad=8,
    )

    xs = [b[0] - b[2] / 2 for b in ALL_BOXES] + [b[0] + b[2] / 2 for b in ALL_BOXES]
    ys = [b[1] - b[3] / 2 for b in ALL_BOXES] + [b[1] + b[3] / 2 for b in ALL_BOXES]
    xmin, xmax = min(xs) - 1.0, max(xs + [COL_X[2] + 8.5]) + 0.5
    ymin, ymax = min(ys) - 1.0, max(ys) + 1.0
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")

    title_in = 0.55
    figw = (xmax - xmin) / DATA_PER_INCH
    figh = (ymax - ymin) / DATA_PER_INCH + title_in
    fig.set_size_inches(figw, figh)
    fig.subplots_adjust(left=0.003, right=0.997, top=1 - title_in / figh, bottom=0.003)

    out_png = os.path.join(OUT_DIR, "main_execution_flowchart_v2.png")
    fig.savefig(out_png, dpi=150)
    print(f"wrote {out_png}  ({figw:.1f}in x {figh:.1f}in, {len(ALL_BOXES)} boxes)")


if __name__ == "__main__":
    main()
