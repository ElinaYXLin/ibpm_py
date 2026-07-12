"""
generate_execution_flowchart.py

Hand-traced "what calls what" flowchart for actually running a simulation
with py/ibpm.py's main() -- i.e. the numerical algorithm, not just the
import graph (see generate_module_graph.py for the mechanical version).

Layout: a strict left-to-right layered call graph. Column 0 is every
statement of interest in ibpm.py's main(), TOP TO BOTTOM IN ACTUAL SOURCE
ORDER. Call depth increases left to right. Every box carries a `file:line`
tag (top-left, outside the box) citing exactly where the code IT SHOWS is
written -- not where it's called from. Every rightward ARROW is additionally
labeled with the file:line of the CALL SITE it leaves from (a line inside
the source box's own code). Arrows are drawn as straight lines wherever
geometrically possible; where a straight line would pass behind an unrelated
box, that portion is rendered dashed (still visible, clearly marked as
"passing behind", not implying a connection to that box). The only
intentionally backward (right-to-left) arrows are the two dashed loop-backs
for this code's two real loops (outer numSteps, inner RK/AB2 substep).

Box size is content-driven (not a fixed per-column size): every box uses the
SAME font size, and its width/height are computed from its own text so
there's no wasted interior whitespace and no overflow.

Nothing here is inferred: each box/edge was found by opening the cited file
and reading the code. See flowchart/README.md for how to re-verify each box
yourself, and flowchart/execution_flow_refs.csv for the same citations in a
table.

Usage:
    python3 flowchart/generate_execution_flowchart.py

Output:
    flowchart/output/main_execution_flowchart.png
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

COLOR_DRIVER = "#c0392b"      # ibpm.py:main() -- column 0
COLOR_SOLVER = "#8e44ad"      # ib_solver.py
COLOR_PROJ = "#2980b9"        # projection_solver.py
COLOR_MODEL = "#16a085"       # navier_stokes_model.py
COLOR_LOWLEVEL = "#7f8c8d"    # elliptic_solver.py / cholesky_solver.py / conjugate_gradient_solver.py
COLOR_GEOM = "#c2185b"        # geometry.py / rigid_body.py
COLOR_IO = "#546e7a"          # logger.py / output_*.py
COLOR_NOTE = "#d35400"        # polymorphic dispatch / SFD-only / conditional notes

LEGEND = [
    (COLOR_DRIVER, "ibpm.py (driver / main())"),
    (COLOR_SOLVER, "ib_solver.py (IBSolver + subclasses)"),
    (COLOR_PROJ, "projection_solver.py (ProjectionSolver.solve)"),
    (COLOR_MODEL, "navier_stokes_model.py (NavierStokesModel)"),
    (COLOR_LOWLEVEL, "elliptic_solver.py / cholesky_solver.py / conjugate_gradient_solver.py"),
    (COLOR_GEOM, "geometry.py / rigid_body.py"),
    (COLOR_IO, "logger.py / output_*.py"),
    (COLOR_NOTE, "polymorphic dispatch or conditional-only note"),
]

FILE_IBPM = "ibpm.py"
FILE_IBS = "ib_solver.py"
FILE_NSM = "navier_stokes_model.py"
FILE_PROJ = "projection_solver.py"
FILE_CHOL_CG = "cholesky_solver.py | conjugate_gradient_solver.py"
FILE_GEOM = "geometry.py"
FILE_RB = "rigid_body.py"
FILE_LOGGER = "logger.py"
FILE_OUT_TEC = "output_tecplot.py"
FILE_OUT_RESTART = "output_restart.py"
FILE_OUT_FORCE = "output_force.py"
FILE_OUT_ENERGY = "output_energy.py"

# column x-centers (call depth 0 = ibpm.py, increasing = deeper calls);
# spacing is generous because box WIDTH is now content-driven (varies per
# box, not fixed per column) -- this just needs to comfortably fit the
# widest box likely to appear in each column plus arrow-label room.
COL_X = [3.0, 11.5, 20.0, 29.0, 38.0, 47.0, 56.0, 65.0, 74.0]

BODY_FS = 7.6      # ONE font size for every box's body text
TAG_FS = 6.2        # ONE font size for every file:line tag
ARROW_LABEL_FS = 6.4

# Empirical text-metric constants (data units per point), calibrated so
# auto-sized boxes snugly fit their text at BODY_FS without overflow or
# excess padding -- tuned against this diagram's actual figsize/DPI.
DATA_PER_INCH = 1.55
_PT_PER_INCH = 72.0
_CHAR_W_PT_FRAC = 0.56          # avg glyph advance width as a fraction of font size, in points
_LINE_H_PT_FRAC = 1.32          # line height as a fraction of font size, in points


def _data_per_pt():
    return DATA_PER_INCH / _PT_PER_INCH


def _text_size_data(lines, fontsize):
    dpp = _data_per_pt()
    max_chars = max((len(s) for s in lines), default=0)
    w = max_chars * fontsize * _CHAR_W_PT_FRAC * dpp
    h = len(lines) * fontsize * _LINE_H_PT_FRAC * dpp
    return w, h


ALL_BOXES = []       # every (x, y, w, h) placed so far, for arrow/box-crossing detection
PENDING_ARROWS = []   # deferred arrow specs, drawn after all boxes exist


def box(ax, col, y, lines, edgecolor, file, tag_lines, facecolor="white",
        linewidth=1.4, zorder=2, style="round,pad=0.02,rounding_size=0.07",
        x=None):
    """Draw one flowchart box, auto-sized to its text at the single global
    BODY_FS font size. `file`/`tag_lines` cite exactly where the code THIS
    BOX SHOWS is written (its own file:line), rendered as a small tag above
    the box's top-left corner -- outside the box, so it never collides with
    an incoming arrow (which lands at top-center) or the body text."""
    cx = COL_X[col] if x is None else x
    tw, th = _text_size_data(lines, BODY_FS)
    w = max(tw + 0.45, 1.7)
    h = max(th + 0.30, 0.62)
    b = FancyBboxPatch((cx - w / 2, y - h / 2), w, h, boxstyle=style,
                         linewidth=linewidth, edgecolor=edgecolor,
                         facecolor=facecolor, zorder=zorder)
    ax.add_patch(b)
    ax.text(cx, y, "\n".join(lines), ha="center", va="center", fontsize=BODY_FS,
             color="#222222", zorder=zorder + 1, linespacing=_LINE_H_PT_FRAC)
    tag = f"{file}:{tag_lines}" if tag_lines else file
    ax.text(cx - w / 2 + 0.10, y + h / 2 + 0.06, tag, ha="left", va="bottom",
             fontsize=TAG_FS, style="italic", color=edgecolor, zorder=zorder + 1)
    rec = (cx, y, w, h)
    ALL_BOXES.append(rec)
    return rec


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


def call_arrow(ax, src, dst, line_label, color="#333333", lw=1.2,
               mutation_scale=11, label_frac=0.52, fontsize=None):
    """Queue an arrow for an actual function call: src (shallower column) ->
    dst (strictly deeper column), labeled with the file:line of the call
    site (a line inside src's own code) -- NOT the line dst is defined on.
    Drawn later, straight, once every box exists (see flush_arrows)."""
    assert dst[0] > src[0], "call_arrow must always advance to a deeper (righter) column"
    PENDING_ARROWS.append(dict(kind="call", src=src, dst=dst, label=line_label, color=color,
                                lw=lw, mutation_scale=mutation_scale, label_frac=label_frac,
                                fontsize=fontsize or ARROW_LABEL_FS))


def seq_arrow(ax, src, dst, color, lw=1.3, ls="solid"):
    """Queue a top-to-bottom sequencing arrow within the same column (or at
    least never moving left) -- "what happens next" at the same call depth."""
    assert dst[0] >= src[0], "seq_arrow must never move to a shallower (lefter) column"
    PENDING_ARROWS.append(dict(kind="seq", src=src, dst=dst, color=color, lw=lw, ls=ls))


def loop_arrow(ax, src, dst, color, x_offset, label, fontsize=7.0, ls="dashed"):
    """Dashed feedback loop routed through a vertical line to the side --
    the ONLY intentional backward-pointing construct in this diagram,
    reserved for the two real loops (outer numSteps, inner substep). Drawn
    immediately (not deferred): loop arrows are routed outside the box
    columns entirely, so they can't cross through other boxes."""
    lx = src[0] + x_offset
    a1 = FancyArrowPatch(left(src) if x_offset < 0 else right(src), (lx, src[1] - src[3] / 2),
                          arrowstyle="-", mutation_scale=1, color=color, linewidth=1.1, zorder=1)
    a2 = FancyArrowPatch((lx, src[1] - src[3] / 2), (lx, dst[1]), arrowstyle="-|>",
                          mutation_scale=11, color=color, linewidth=1.1, linestyle=ls, zorder=1)
    a3 = FancyArrowPatch((lx, dst[1]), left(dst) if x_offset < 0 else right(dst), arrowstyle="-",
                          mutation_scale=1, color=color, linewidth=1.1, zorder=1)
    for a in (a1, a2, a3):
        ax.add_patch(a)
    ax.text(lx + (0.35 if x_offset > 0 else -0.35), (src[1] - src[3] / 2 + dst[1]) / 2, label,
             ha="center", va="center", fontsize=fontsize, color=color, rotation=90)


def _clip_segment_to_rect(p1, p2, rect):
    """Liang-Barsky clip of segment p1->p2 against axis-aligned rect
    (x0,y0,x1,y1). Returns (tmin,tmax) in [0,1] where the segment is INSIDE
    the rect's interior, or None if it never enters."""
    x0, y0, x1, y1 = rect
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    tmin, tmax = 0.0, 1.0
    for p_, q_ in ((-dx, p1[0] - x0), (dx, x1 - p1[0]), (-dy, p1[1] - y0), (dy, y1 - p1[1])):
        if abs(p_) < 1e-12:
            if q_ < 0:
                return None
            continue
        t = q_ / p_
        if p_ < 0:
            if t > tmax:
                return None
            tmin = max(tmin, t)
        else:
            if t < tmin:
                return None
            tmax = min(tmax, t)
    if tmin >= tmax:
        return None
    return (max(tmin, 0.0), min(tmax, 1.0))


def _dashed_intervals(p1, p2, exclude):
    """Merged, sorted list of (t0,t1) in [0,1] along segment p1->p2 where it
    passes through the interior of any placed box OTHER than those in
    `exclude` (the arrow's own endpoints)."""
    intervals = []
    for b in ALL_BOXES:
        if any(b is e for e in exclude):
            continue
        x, y, w, h = b
        hit = _clip_segment_to_rect(p1, p2, (x - w / 2, y - h / 2, x + w / 2, y + h / 2))
        if hit is not None and hit[1] > 1e-6 and hit[0] < 1 - 1e-6:
            intervals.append(hit)
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for t0, t1 in intervals[1:]:
        if t0 <= merged[-1][1] + 1e-6:
            merged[-1][1] = max(merged[-1][1], t1)
        else:
            merged.append([t0, t1])
    return [(a, b_) for a, b_ in merged]


def _draw_straight_arrow(ax, p1, p2, color, lw, mutation_scale, dashed_intervals):
    """Draw p1->p2 as a straight line, split into solid/dashed segments per
    `dashed_intervals`; only the final segment carries the arrowhead."""
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    segs, cursor = [], 0.0
    for t0, t1 in dashed_intervals:
        if t0 > cursor + 1e-9:
            segs.append((cursor, t0, False))
        segs.append((max(t0, cursor), t1, True))
        cursor = t1
    if cursor < 1.0 - 1e-9:
        segs.append((cursor, 1.0, False))
    if not segs:
        segs = [(0.0, 1.0, False)]
    for i, (t0, t1, dashed) in enumerate(segs):
        q1 = (p1[0] + dx * t0, p1[1] + dy * t0)
        q2 = (p1[0] + dx * t1, p1[1] + dy * t1)
        is_last = i == len(segs) - 1
        a = FancyArrowPatch(q1, q2, arrowstyle="-|>" if is_last else "-",
                             mutation_scale=mutation_scale if is_last else 1,
                             color=color, linewidth=lw,
                             linestyle="dashed" if dashed else "solid", zorder=1)
        ax.add_patch(a)


def flush_arrows(ax):
    """Draw every queued arrow now that ALL_BOXES is complete, so
    box-crossing detection sees the final layout."""
    for spec in PENDING_ARROWS:
        src, dst = spec["src"], spec["dst"]
        if spec["kind"] == "call":
            p1, p2 = right(src), left(dst)
        else:
            p1, p2 = bottom(src), top(dst)
        intervals = _dashed_intervals(p1, p2, exclude=(src, dst))
        _draw_straight_arrow(ax, p1, p2, spec["color"], spec["lw"],
                              spec.get("mutation_scale", 11), intervals)
        label = spec.get("label")
        if label:
            frac = spec.get("label_frac", 0.5)
            lx, ly = p1[0] + (p2[0] - p1[0]) * frac, p1[1] + (p2[1] - p1[1]) * frac
            ax.text(lx, ly, label, ha="center", va="center", fontsize=spec["fontsize"],
                     color=spec["color"], zorder=3,
                     bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                               edgecolor=spec["color"], linewidth=0.6, alpha=0.95))


def draw_movebodies_subtree(ax, source, col0, label):
    """geometry.moveBodies(t) -> loop bodies -> RigidBody.moveBody(t).
    Shared helper: this exact call chain is reached from two places in
    main() (directly, and again inside NavierStokesModel.updateOperators)."""
    g = box(ax, col0, source[1], ["geometry.moveBodies(t)"], COLOR_GEOM, FILE_GEOM, "128-131")
    call_arrow(ax, source, g, label, color=COLOR_GEOM)
    rb = box(ax, col0 + 1, source[1], ["for body in bodies:", "body.moveBody(t)"], COLOR_GEOM,
             FILE_GEOM, "130-131")
    call_arrow(ax, g, rb, "130-131", color=COLOR_GEOM)
    return g, rb


def draw_update_operators_subtree(ax, source, col0, label):
    """model.updateOperators(t) -> {baseFlow.moveFlow(t) if bfTimeDependent,
    geometry.moveBodies(t)+regularizer.update() if geTimeDependent}."""
    uo = box(ax, col0, source[1], ["NavierStokesModel", ".updateOperators(t)"], COLOR_MODEL,
             FILE_NSM, "110-116")
    call_arrow(ax, source, uo, label, color=COLOR_MODEL)

    mf = box(ax, col0 + 1, source[1] + 0.85, ["if bfTimeDependent():", "baseFlow.moveFlow(t)"],
             COLOR_MODEL, FILE_NSM, "112-113")
    call_arrow(ax, uo, mf, "112-113", color=COLOR_MODEL)

    mb_note = box(ax, col0 + 1, source[1] - 0.95,
                  ["if geTimeDependent():", "geometry.moveBodies(t)", "regularizer.update()"],
                  COLOR_MODEL, FILE_NSM, "114-116")
    call_arrow(ax, uo, mb_note, "114-116", color=COLOR_MODEL)

    draw_movebodies_subtree(ax, mb_note, col0 + 2, "115")
    return uo


def draw_refresh_state_subtree(ax, source, col0, label):
    """model.refreshState(x) -> computeFlux(omega, q) ->
    computeFluxWithoutBaseFlow(omega, q) [vorticityToStreamfunction+Curl]
    + q += baseFlow.getFlux()."""
    rs = box(ax, col0, source[1], ["NavierStokesModel", ".refreshState(x)"], COLOR_MODEL,
             FILE_NSM, "168-170")
    call_arrow(ax, source, rs, label, color=COLOR_MODEL)

    cf = box(ax, col0 + 1, source[1], ["computeFlux(omega, q)"], COLOR_MODEL, FILE_NSM, "158-166")
    call_arrow(ax, rs, cf, "170", color=COLOR_MODEL)

    cfwbf = box(ax, col0 + 2, source[1] + 0.8,
                ["computeFluxWithoutBaseFlow(omega,q):", "psi = vorticityToStreamfunction(omega)",
                 "Curl(psi, q)"], COLOR_MODEL, FILE_NSM, "145-148")
    call_arrow(ax, cf, cfwbf, "165", color=COLOR_MODEL)

    getflux = box(ax, col0 + 2, source[1] - 0.8, ["q += baseFlow.getFlux()"], COLOR_MODEL,
                  FILE_NSM, "166")
    call_arrow(ax, cf, getflux, "166", color=COLOR_MODEL)
    return rs


def draw_logger_dooutput_subtree(ax, source, col0, label):
    """logger.doOutput(q_potential, x) -> for each registered entry due
    this step: entry.output.doOutput(...) -> the concrete Output
    subclass's own doOutput (each just writes its own file)."""
    ld = box(ax, col0, source[1], ["Logger.doOutput", "(q_potential, x)"], COLOR_IO, FILE_LOGGER,
             "45-71")
    call_arrow(ax, source, ld, label, color=COLOR_IO)

    loop = box(ax, col0 + 1, source[1], ["for entry in _outputs:", "if entry.shouldBeCalled(x):",
               "entry.output.doOutput(...)"], COLOR_IO, FILE_LOGGER, "66-70")
    call_arrow(ax, ld, loop, "66-70", color=COLOR_IO)

    targets = [
        (["OutputTecplot.doOutput:", "writes .plt Tecplot file"], FILE_OUT_TEC, "59"),
        (["OutputRestart.doOutput:", "writes .bin restart file"], FILE_OUT_RESTART, "29"),
        (["OutputForce.doOutput:", "x.computeNetForce() -> writes .force line"], FILE_OUT_FORCE, "45-100"),
        (["OutputEnergy.doOutput:", "writes .energy line"], FILE_OUT_ENERGY, "46"),
    ]
    yy = source[1] + 1.5
    for txt, f, ln in targets:
        tbox = box(ax, col0 + 2, yy, txt, COLOR_IO, f, ln)
        call_arrow(ax, loop, tbox, "69", color=COLOR_IO, fontsize=6.0)
        yy -= 1.15
    return ld


def draw_legend(ax, x, y_top):
    """Small color-key panel: one swatch + label per source-file color."""
    row_h = 0.72
    swatch_w, swatch_h = 0.5, 0.32
    title_y = y_top
    ax.text(x, title_y, "Legend: box color = source file", ha="left", va="bottom",
             fontsize=8.6, fontweight="bold", color="#222222")
    y = title_y - 0.55
    for color, label in LEGEND:
        ax.add_patch(Rectangle((x, y - swatch_h / 2), swatch_w, swatch_h,
                                 facecolor="white", edgecolor=color, linewidth=1.6, zorder=2))
        ax.text(x + swatch_w + 0.18, y, label, ha="left", va="center", fontsize=7.4,
                 color="#222222", zorder=2)
        y -= row_h
    ax.text(x, y - 0.15,
            "Solid arrow = call (labeled with the CALL SITE's file:line).\n"
            "Dashed segment = arrow passing behind an unrelated box (not connected to it).\n"
            "Dashed loop-back = the outer numSteps loop / inner substep loop (the only\n"
            "backward-pointing arrows in this diagram).",
            ha="left", va="top", fontsize=6.9, color="#444444", style="italic")


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 10))  # placeholder; resized at the end to fit content

    # ===================== COLUMN 0: ibpm.py main(), in ACTUAL SOURCE ORDER =====================
    y = 60.0
    dy = 1.6
    n01 = box(ax, 0, y, ["ParmParser(argc, argv)"], COLOR_DRIVER, FILE_IBPM, "166")
    y -= dy
    n02 = box(ax, 0, y, ["Grid(nx, ny, ngrid, ...)"], COLOR_DRIVER, FILE_IBPM, "280")
    y -= dy
    n03 = box(ax, 0, y, ["geom.load(geomFile)"], COLOR_DRIVER, FILE_IBPM, "285")
    y -= dy
    n04 = box(ax, 0, y, ["BaseFlow(grid, mag, alpha)"], COLOR_DRIVER, FILE_IBPM, "298")
    y -= dy
    n05 = box(ax, 0, y, ["Build model + solver", "(branches on -model flag)"],
              COLOR_DRIVER, FILE_IBPM, "318-347")
    n05_y = y

    y -= dy + 1.1
    n06 = box(ax, 0, y, ["x = State(grid, ...);", "x.load(icFile)"], COLOR_DRIVER,
              FILE_IBPM, "361-369")

    y -= dy + 0.15
    n07 = box(ax, 0, y, ["geom.moveBodies(x.time)"], COLOR_DRIVER, FILE_IBPM, "394")
    draw_movebodies_subtree(ax, n07, 1, "394")

    y -= dy + 0.15
    n08 = box(ax, 0, y, ["model.init()"], COLOR_DRIVER, FILE_IBPM, "397")
    mi = box(ax, 1, y, ["NavierStokesModel.init()"], COLOR_MODEL, FILE_NSM, "66-72")
    call_arrow(ax, n08, mi, "397", color=COLOR_MODEL)
    reg = box(ax, 2, y, ["if not hasBeenInitialized:", "regularizer.update()"], COLOR_MODEL,
              FILE_NSM, "68,71")
    call_arrow(ax, mi, reg, "68-71", color=COLOR_MODEL)

    y -= dy + 0.15
    n09 = box(ax, 0, y, ["solver.load(...) or", "solver.init()+save(...)"],
              COLOR_DRIVER, FILE_IBPM, "400-405")
    n09_y = y
    sl = box(ax, 1, y + 0.9, ["solver.load(basename)"], COLOR_DRIVER, FILE_IBPM, "400")
    call_arrow(ax, n09, sl, "400", color=COLOR_SOLVER)
    si = box(ax, 1, y - 0.9, ["if load failed:", "solver.init(); solver.save(basename)"],
             COLOR_DRIVER, FILE_IBPM, "403-405")
    call_arrow(ax, n09, si, "403-405", color=COLOR_SOLVER)
    ibs_loop = box(ax, 2, y, ["for i in range(nsteps):", "_solver[i].load/init/save(...)"],
                   COLOR_SOLVER, FILE_IBS, "112-139")
    call_arrow(ax, sl, ibs_loop, "122-128", color=COLOR_SOLVER)
    call_arrow(ax, si, ibs_loop, "112-114,130-139", color=COLOR_SOLVER)
    chol = box(ax, 3, y, ["CholeskySolver.init():", "computeMatrixM()", "-> computeFactorization()",
               "ConjugateGradientSolver:", "init() is a no-op (base class)"], COLOR_LOWLEVEL,
               FILE_CHOL_CG, "110-122 | n/a")
    call_arrow(ax, ibs_loop, chol, "110-122", color=COLOR_SOLVER)

    y -= dy + 0.4
    n10 = box(ax, 0, y, ["model.updateOperators(x.time)"], COLOR_DRIVER, FILE_IBPM, "411")
    draw_update_operators_subtree(ax, n10, 1, "411")

    y -= dy + 1.35
    n11 = box(ax, 0, y, ["model.refreshState(x)"], COLOR_DRIVER, FILE_IBPM, "412")
    draw_refresh_state_subtree(ax, n11, 1, "412")

    y -= dy + 0.4
    n12 = box(ax, 0, y, ["Setup Output objects", "(Tecplot/Restart/Force/Energy);",
              "logger.addOutput(...)"], COLOR_DRIVER, FILE_IBPM, "417-443")

    y -= dy + 0.2
    n13 = box(ax, 0, y, ["logger.init()"], COLOR_DRIVER, FILE_IBPM, "445")
    li = box(ax, 1, y, ["Logger.init()"], COLOR_IO, FILE_LOGGER, "73-80")
    call_arrow(ax, n13, li, "445", color=COLOR_IO)
    li_loop = box(ax, 2, y, ["for entry in _outputs:", "entry.output.init()"], COLOR_IO,
                  FILE_LOGGER, "76-79")
    call_arrow(ax, li, li_loop, "76-79", color=COLOR_IO)

    y -= dy + 0.15
    n14 = box(ax, 0, y, ["logger.doOutput(q_potential, x)", "(initial output)"],
              COLOR_DRIVER, FILE_IBPM, "446")
    draw_logger_dooutput_subtree(ax, n14, 1, "446")

    y -= dy + 1.6
    loop_top_y = y
    n15 = box(ax, 0, y, ["for i in 1..numSteps:"], COLOR_DRIVER, FILE_IBPM, "449",
              facecolor="#fdf2ef")

    y -= dy
    n16 = box(ax, 0, y, ["solver.advance(x)"], COLOR_SOLVER, FILE_IBPM, "452",
              facecolor="#f5eefb")

    # ================= BRANCH B: n16 "solver.advance(x)" (the big one) =================
    b1_y = n16[1]
    b1 = box(ax, 1, b1_y, ["IBSolver.advance(x)"], COLOR_SOLVER, FILE_IBS, "171-183")
    call_arrow(ax, n16, b1, "452", color=COLOR_SOLVER, lw=1.8, mutation_scale=15)

    b2 = box(ax, 2, b1_y, ["for i in range(nsteps):", "nonlinear = self.N(x)"], COLOR_SOLVER,
             FILE_IBS, "175,177")
    call_arrow(ax, b1, b2, "175,177", color=COLOR_SOLVER)

    b3 = box(ax, 3, b1_y + 1.7,
             ["N(x) polymorphic override:", "Nonlinear: Curl(Cross(q,ω))  [241]",
              "Linearized: Curl(Cross(q0,ω)+Cross(q,ω0))  [261]",
              "Adjoint: Lap(Curl(Cross(q0,q)))-Curl(Cross(q,ω0))  [282]",
              "Periodic: as Linearized, ω0->ω0periodic[k]  [311]",
              "SFD: as Nonlinear, - chi*(ω-ωhat)  [347]"],
             COLOR_NOTE, FILE_IBS, "241,261,282,311,347")
    call_arrow(ax, b2, b3, "177", color=COLOR_SOLVER)

    b4 = box(ax, 3, b1_y - 1.7, ["advanceSubstep(x, nonlinear, i)"], COLOR_SOLVER, FILE_IBS,
             "185-216")
    call_arrow(ax, b2, b4, "180", color=COLOR_SOLVER)

    b5 = box(ax, 4, b1_y + 2.9, ["if isTimeDependent():", "model.updateOperators(t)"],
             COLOR_SOLVER, FILE_IBS, "187-188")
    call_arrow(ax, b4, b5, "187-188", color=COLOR_SOLVER)
    draw_update_operators_subtree(ax, b5, 5, "188")

    b6 = box(ax, 4, b1_y + 0.65,
             ["a = Laplacian(ω)*coef", "+ coef*nonlinear (+bn*Nprev)"], COLOR_SOLVER,
             FILE_IBS, "191-203")
    call_arrow(ax, b4, b6, "191-203", color=COLOR_SOLVER)

    b7 = box(ax, 4, b1_y - 0.85, ["b = model.getConstraints()"], COLOR_SOLVER, FILE_IBS, "206")
    call_arrow(ax, b4, b7, "206", color=COLOR_SOLVER)
    # placed clear ABOVE the wide ProjectionSolver.solve() breakdown (c_labels,
    # drawn below in column 5 across nearly all of branch B's height) so the
    # two column-5 subtrees don't collide
    gc = box(ax, 5, b1_y + 1.5, ["geometry.getVelocities()", "- toBoundary(baseFlow.getFlux())"],
             COLOR_MODEL, FILE_NSM, "105-107")
    call_arrow(ax, b7, gc, "105-107", color=COLOR_SOLVER)

    b8 = box(ax, 4, b1_y - 2.25, ["self._solver[i].solve(", "a, b, x.omega, x.f)"], COLOR_SOLVER,
             FILE_IBS, "209")
    call_arrow(ax, b4, b8, "209", color=COLOR_SOLVER)

    b9 = box(ax, 4, b1_y - 5.0, ["model.refreshState(x)"], COLOR_SOLVER, FILE_IBS, "212")
    call_arrow(ax, b4, b9, "212", color=COLOR_SOLVER)
    # skip column 5 entirely here (occupied top-to-bottom by the
    # ProjectionSolver.solve() breakdown below) -- start this subtree one
    # column further right instead of colliding with it
    draw_refresh_state_subtree(ax, b9, 6, "212")

    b_sfd = box(ax, 4, b1_y - 6.9, ["[SFDSolver only] also", "integrates filtered state ωhat"],
                COLOR_NOTE, FILE_IBS, "354-386")
    call_arrow(ax, b4, b_sfd, "354-386", color=COLOR_SOLVER)

    # ProjectionSolver.solve() breakdown -- one layer deeper than b8. All six
    # steps below are sequential lines inside ProjectionSolver.solve() itself
    # (projection_solver.py), even though several delegate out to
    # NavierStokesModel/HelmholtzSolver/CholeskySolver -- the "->" in each
    # box's text names that delegate, but the line cited (and the box's own
    # file tag) is always projection_solver.py's.
    c_y = b8[1] + 2.2
    c_labels_lines = [
        (["Ainv(a,ω*) ->", "HelmholtzSolver.solve"], "128,144-146"),
        (["C(ω*,rhs) -> model.C:", "Poisson.solve+Curl+toBoundary"], "132"),
        (["Minv(rhs,f): Cholesky", "back-sub | CG iterate"], "134"),
        (["B(f,c) -> model.B:", "toFlux(f); Curl(q,ω)"], "138"),
        (["Ainv(c,c) ->", "HelmholtzSolver.solve"], "139"),
        (["ω.assign(ω* - c)"], "140"),
    ]
    for lines_, lbl in c_labels_lines:
        cc = box(ax, 5, c_y, lines_, COLOR_PROJ, FILE_PROJ, lbl)
        call_arrow(ax, b8, cc, lbl, color=COLOR_PROJ, fontsize=6.0, label_frac=0.4)
        c_y -= 1.5

    branchB_bottom = min(b_sfd[1] - b_sfd[3] / 2, c_y + 1.5 - 0.95 / 2) - 0.6

    b10 = box(ax, 2, branchB_bottom - 1.2, ["after all substeps:", "x.time += dt; x.timestep += 1"],
              COLOR_SOLVER, FILE_IBS, "182-183")
    # SAME column as b2 (both are statements in IBSolver.advance's own body,
    # b10 sequentially after the substep loop) -- straight DOWN, never left.
    seq_arrow(ax, b2, b10, COLOR_SOLVER, lw=1.0, ls="dashed")
    loop_arrow(ax, b4, b2, COLOR_SOLVER, x_offset=1.8, label="repeat per\nsubstep i", fontsize=6.5)

    # ===================== back in column 0, below branch B entirely =====================
    y = branchB_bottom - dy - 1.0
    n17 = box(ax, 0, y, ["x.computeNetForce()"], COLOR_DRIVER, FILE_IBPM, "453")

    y -= dy + 0.15
    n18 = box(ax, 0, y, ["logger.doOutput(q_potential, x)"], COLOR_DRIVER, FILE_IBPM, "462")
    draw_logger_dooutput_subtree(ax, n18, 1, "462")
    loop_bottom_y = y

    y -= dy + 1.6
    n19 = box(ax, 0, y, ["logger.cleanup()"], COLOR_DRIVER, FILE_IBPM, "480")
    lc = box(ax, 1, y, ["Logger.cleanup()"], COLOR_IO, FILE_LOGGER, "82-89")
    call_arrow(ax, n19, lc, "480", color=COLOR_IO)
    lc_loop = box(ax, 2, y, ["for entry in _outputs:", "entry.output.cleanup()"], COLOR_IO,
                  FILE_LOGGER, "84-88")
    call_arrow(ax, lc, lc_loop, "84-88", color=COLOR_IO)

    # ===================== column-0 sequencing arrows (straight down only) =====================
    for p, q in [(n01, n02), (n02, n03), (n03, n04), (n04, n05), (n05, n06),
                 (n06, n07), (n07, n08), (n08, n09), (n09, n10), (n10, n11),
                 (n11, n12), (n12, n13), (n13, n14), (n14, n15), (n15, n16),
                 (n16, n17), (n17, n18), (n18, n19)]:
        seq_arrow(ax, p, q, COLOR_DRIVER)

    # ===================== BRANCH A: n05 "build model + solver" =====================
    a1 = box(ax, 1, n05_y + 1.0, ["NavierStokesModel", "(grid, geom, Re[, q_potential])"],
             COLOR_MODEL, FILE_NSM, "37-64")
    a2 = box(ax, 1, n05_y - 1.0, ["<Solver>(grid, model, dt,", "scheme, ...) -- Nonlinear/",
             "Linearized/Adjoint/Periodic/SFD"], COLOR_SOLVER, FILE_IBS, "84-104,230-408")
    call_arrow(ax, n05, a1, "318,324,330,340,346")
    call_arrow(ax, n05, a2, "319,325,331,341,347")

    a3 = box(ax, 2, a2[1], ["IBSolver.__init__", "-> self.createAllSolvers()"], COLOR_SOLVER,
             FILE_IBS, "84-104")
    call_arrow(ax, a2, a3, "84-104")

    a4 = box(ax, 3, a2[1], ["createAllSolvers(): loop substeps", "-> createSolver(beta)"],
             COLOR_SOLVER, FILE_IBS, "141-145")
    call_arrow(ax, a3, a4, "104")

    a5 = box(ax, 4, a2[1], ["CholeskySolver (stationary)", "ConjugateGradientSolver (moving)"],
             COLOR_LOWLEVEL, FILE_CHOL_CG, "153-165")
    call_arrow(ax, a4, a5, "141-165")

    # ===================== loop-back arrows (the only intentionally backward arrows) ============
    loop_arrow(ax, n18, n15, COLOR_DRIVER, x_offset=-2.4, label="repeat\nnumSteps\ntimes")

    # ===================== legend, in the empty space above/right of branch A =====================
    draw_legend(ax, COL_X[6], n01[1] + 0.3)

    # ===================== finalize: draw all deferred arrows, then size the figure =====================
    flush_arrows(ax)

    ax.set_title(
        "ibpm_py -- execution flow of a simulation run (py/ibpm.py: main())\n"
        "Column 0 = ibpm.py's own statements, TOP TO BOTTOM IN ACTUAL SOURCE ORDER. Every box "
        "is tagged file:line for the code IT shows; every rightward arrow is additionally "
        "labeled with the CALL SITE's file:line. Straight arrows; dashed where they pass "
        "behind an unrelated box. See legend for colors.", fontsize=12.5, pad=8,
    )

    xs = [b[0] - b[2] / 2 for b in ALL_BOXES] + [b[0] + b[2] / 2 for b in ALL_BOXES]
    ys = [b[1] - b[3] / 2 for b in ALL_BOXES] + [b[1] + b[3] / 2 for b in ALL_BOXES]
    xmin, xmax = min(xs) - 1.0, max(xs) + 1.0
    ymin, ymax = min(ys) - 1.0, max(ys) + 1.0
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")

    title_in = 0.55
    figw = (xmax - xmin) / DATA_PER_INCH
    figh = (ymax - ymin) / DATA_PER_INCH + title_in
    fig.set_size_inches(figw, figh)
    fig.subplots_adjust(left=0.003, right=0.997, top=1 - title_in / figh, bottom=0.003)

    out_png = os.path.join(OUT_DIR, "main_execution_flowchart.png")
    fig.savefig(out_png, dpi=150)
    print(f"wrote {out_png}  ({figw:.1f}in x {figh:.1f}in, {len(ALL_BOXES)} boxes, "
          f"{len(PENDING_ARROWS)} arrows)")


if __name__ == "__main__":
    main()
