"""
generate_execution_flowchart.py

Hand-traced "what calls what" flowchart for actually running a simulation
with py/ibpm.py's main() -- i.e. the numerical algorithm, not just the
import graph (see generate_module_graph.py for the mechanical version).

Layout: a strict left-to-right layered call graph. Column 0 is every
statement of interest in ibpm.py's main(), TOP TO BOTTOM IN ACTUAL SOURCE
ORDER -- so "where does ibpm.py call into X" is always answered by scanning
straight down column 0. Every rightward arrow is labeled with the file:line
of the CALL SITE (i.e. the line in the box the arrow leaves, not the line
the callee is defined on). Call depth increases left to right: column 1 is
what column-0 code calls directly, column 2 is what column-1 code calls,
and so on -- and this now goes as deep as the call chain actually goes
(e.g. model.updateOperators() -> geometry.moveBodies() -> RigidBody.moveBody(),
or logger.doOutput() -> each registered Output.doOutput()), rather than
stopping after one hop. Every arrow either goes strictly left-to-right (a
call) or straight down within one column (sequencing); the only
right-to-left segments are the explicitly dashed loop-back arrows for the
two real loops in this code (the outer numSteps loop, the inner RK/AB2
substep loop) -- nothing else should ever point backward.

Nothing here is inferred: each box/edge was found by opening the cited file
and reading the call. See flowchart/README.md for how to re-verify each box
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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

COLOR_DRIVER = "#c0392b"      # ibpm.py:main() -- column 0
COLOR_SOLVER = "#8e44ad"      # ib_solver.py
COLOR_PROJ = "#2980b9"        # projection_solver.py
COLOR_MODEL = "#16a085"       # navier_stokes_model.py
COLOR_LOWLEVEL = "#7f8c8d"    # elliptic_solver.py / cholesky/CG solvers
COLOR_GEOM = "#c2185b"        # geometry.py / rigid_body.py / base_flow.py
COLOR_IO = "#546e7a"          # logger.py / output*.py
COLOR_NOTE = "#d35400"        # polymorphic dispatch / SFD-only / conditional notes

FILE_IBPM = "ibpm.py"
FILE_IBS = "ib_solver.py"
FILE_NSM = "navier_stokes_model.py"
FILE_PROJ = "projection_solver.py"
FILE_CHOL_CG = "cholesky_solver.py | conjugate_gradient_solver.py"
FILE_GEOM = "geometry.py"
FILE_RB = "rigid_body.py"
FILE_BF = "base_flow.py"
FILE_LOGGER = "logger.py"
FILE_OUTPUT = "output_*.py"

# column x-centers (call depth 0 = ibpm.py, increasing = deeper calls)
COL_X = [3.0, 10.6, 18.4, 26.6, 35.0, 43.4, 51.8, 60.2, 68.6]
COL_W = [6.0, 6.6, 7.2, 7.8, 6.6, 6.2, 6.2, 6.2, 6.2]


def box(ax, col, y, lines, edgecolor, file, h=1.0, fontsize=8.0, facecolor="white",
        linewidth=1.4, zorder=2, style="round,pad=0.02,rounding_size=0.07"):
    """Draw one flowchart box. `file` is the py/ source file this box's code
    belongs to; rendered as a small tag above the box's top-left corner,
    outside the box, so it never collides with an incoming arrow (which
    lands at top-center) or body text."""
    x, w = COL_X[col], COL_W[col]
    b = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=style,
                         linewidth=linewidth, edgecolor=edgecolor,
                         facecolor=facecolor, zorder=zorder)
    ax.add_patch(b)
    ax.text(x, y, "\n".join(lines), ha="center", va="center", fontsize=fontsize,
             color="#222222", zorder=zorder + 1, linespacing=1.25)
    ax.text(x - w / 2 + 0.12, y + h / 2 + 0.07, file, ha="left", va="bottom",
             fontsize=6.3, style="italic", color=edgecolor, zorder=zorder + 1)
    return (x, y, w, h)


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
               rad=0.0, mutation_scale=11, label_frac=0.52, fontsize=6.6):
    """Arrow for an actual function call: src (left column) -> dst (strictly
    deeper column), labeled with the file:line of the call site (a line
    inside src's own code) -- NOT the line dst is defined on. dst MUST be in
    a column to the right of src; this never points backward."""
    assert dst[0] > src[0], "call_arrow must always advance to a deeper (righter) column"
    p1, p2 = right(src), left(dst)
    a = FancyArrowPatch(p1, p2, connectionstyle=f"arc3,rad={rad}",
                         arrowstyle="-|>", mutation_scale=mutation_scale,
                         color=color, linewidth=lw, zorder=1)
    ax.add_patch(a)
    if line_label:
        lx = p1[0] + (p2[0] - p1[0]) * label_frac
        ly = p1[1] + (p2[1] - p1[1]) * label_frac
        ax.text(lx, ly, line_label, ha="center", va="center", fontsize=fontsize,
                 color=color, zorder=3,
                 bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                           edgecolor=color, linewidth=0.6, alpha=0.92))


def seq_arrow(ax, src, dst, color, lw=1.3, ls="solid"):
    """Plain top-to-bottom sequencing arrow within the SAME column (or at
    least never moving left) -- used for "what happens next" at the same
    call depth, e.g. consecutive statements in the same function body."""
    assert dst[0] >= src[0], "seq_arrow must never move to a shallower (lefter) column"
    a = FancyArrowPatch(bottom(src), top(dst), arrowstyle="-|>",
                         mutation_scale=11, color=color, linewidth=lw,
                         linestyle=ls, zorder=1)
    ax.add_patch(a)


def loop_arrow(ax, src, dst, color, x_offset, label, fontsize=7.0, ls="dashed"):
    """Dashed feedback loop routed through a vertical line to the side --
    the ONLY intentional backward-pointing construct in this diagram,
    reserved for the two real loops (outer numSteps, inner substep)."""
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


def draw_movebodies_subtree(ax, source, col0, label):
    """geometry.moveBodies(t) -> loop bodies -> RigidBody.moveBody(t).
    Shared helper: this exact call chain is reached from two places in
    main() (directly, and again inside NavierStokesModel.updateOperators)."""
    g = box(ax, col0, source[1], ["geometry.moveBodies(t)"], COLOR_GEOM, FILE_GEOM,
            h=0.85, fontsize=7.3)
    call_arrow(ax, source, g, label, color=COLOR_GEOM)
    rb = box(ax, col0 + 1, source[1], ["for body in bodies:", "body.moveBody(t)"], COLOR_GEOM,
             FILE_RB, h=0.9, fontsize=7.3)
    call_arrow(ax, g, rb, "128-131", color=COLOR_GEOM)
    return g, rb


def draw_update_operators_subtree(ax, source, col0, label):
    """model.updateOperators(t) -> {baseFlow.moveFlow(t) if bfTimeDependent,
    geometry.moveBodies(t)+regularizer.update() if geTimeDependent}."""
    uo = box(ax, col0, source[1], ["NavierStokesModel", ".updateOperators(t)"], COLOR_MODEL,
             FILE_NSM, h=0.95, fontsize=7.4)
    call_arrow(ax, source, uo, label, color=COLOR_MODEL)

    mf = box(ax, col0 + 1, source[1] + 0.75, ["if bfTimeDependent():", "baseFlow.moveFlow(t)"],
             COLOR_GEOM, FILE_BF, h=0.9, fontsize=7.2)
    call_arrow(ax, uo, mf, "112-113", color=COLOR_MODEL, rad=0.15)

    mb_note = box(ax, col0 + 1, source[1] - 0.85,
                  ["if geTimeDependent():", "geometry.moveBodies(t)", "regularizer.update()"],
                  COLOR_MODEL, FILE_NSM, h=1.2, fontsize=7.2)
    call_arrow(ax, uo, mb_note, "114-116", color=COLOR_MODEL, rad=-0.15)

    g, rb = draw_movebodies_subtree(ax, mb_note, col0 + 2, "115")
    return uo


def draw_refresh_state_subtree(ax, source, col0, label):
    """model.refreshState(x) -> computeFlux(omega, q) ->
    computeFluxWithoutBaseFlow(omega, q) [vorticityToStreamfunction+Curl]
    + q += baseFlow.getFlux()."""
    rs = box(ax, col0, source[1], ["NavierStokesModel", ".refreshState(x)"], COLOR_MODEL,
             FILE_NSM, h=0.95, fontsize=7.4)
    call_arrow(ax, source, rs, label, color=COLOR_MODEL)

    cf = box(ax, col0 + 1, source[1], ["computeFlux(omega, q)"], COLOR_MODEL, FILE_NSM,
             h=0.85, fontsize=7.4)
    call_arrow(ax, rs, cf, "170", color=COLOR_MODEL)

    cfwbf = box(ax, col0 + 2, source[1] + 0.7,
                ["computeFluxWithoutBaseFlow(omega,q):", "psi = vorticityToStreamfunction(omega)",
                 "Curl(psi, q)"], COLOR_MODEL, FILE_NSM, h=1.15, fontsize=7.0)
    call_arrow(ax, cf, cfwbf, "165", color=COLOR_MODEL, rad=0.12)

    getflux = box(ax, col0 + 2, source[1] - 0.7, ["q += baseFlow.getFlux()"], COLOR_GEOM,
                  FILE_BF, h=0.85, fontsize=7.3)
    call_arrow(ax, cf, getflux, "166", color=COLOR_MODEL, rad=-0.12)
    return rs


def draw_logger_dooutput_subtree(ax, source, col0, label):
    """logger.doOutput(q_potential, x) -> for each registered entry due
    this step: entry.output.doOutput(...) -> the concrete Output
    subclass's own doOutput (each just writes its own file)."""
    ld = box(ax, col0, source[1], ["Logger.doOutput", "(q_potential, x)"], COLOR_IO, FILE_LOGGER,
             h=0.95, fontsize=7.4)
    call_arrow(ax, source, ld, label, color=COLOR_IO)

    loop = box(ax, col0 + 1, source[1], ["for entry in _outputs:", "if entry.shouldBeCalled(x):",
               "entry.output.doOutput(...)"], COLOR_IO, FILE_LOGGER, h=1.15, fontsize=7.1)
    call_arrow(ax, ld, loop, "66-70", color=COLOR_IO)

    targets = [
        ("OutputTecplot.doOutput:", "writes .plt Tecplot file"),
        ("OutputRestart.doOutput:", "writes .bin restart file"),
        ("OutputForce.doOutput:", "x.computeNetForce() -> writes .force line"),
        ("OutputEnergy.doOutput:", "writes .energy line"),
    ]
    yy = source[1] + 1.35
    for t1, t2 in targets:
        tbox = box(ax, col0 + 2, yy, [t1, t2], COLOR_IO, FILE_OUTPUT, h=0.85, fontsize=6.8)
        call_arrow(ax, loop, tbox, "69", color=COLOR_IO, rad=0.0, fontsize=6.2)
        yy -= 1.05
    return ld


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(48, 56))

    # ===================== COLUMN 0: ibpm.py main(), in ACTUAL SOURCE ORDER =====================
    y = 60.0
    dy = 1.9
    n01 = box(ax, 0, y, ["ParmParser(argc, argv)", "[166]"], COLOR_DRIVER, FILE_IBPM)
    y -= dy
    n02 = box(ax, 0, y, ["Grid(nx, ny, ngrid, ...)", "[280]"], COLOR_DRIVER, FILE_IBPM)
    y -= dy
    n03 = box(ax, 0, y, ["geom.load(geomFile)", "[285]"], COLOR_DRIVER, FILE_IBPM)
    y -= dy
    n04 = box(ax, 0, y, ["BaseFlow(grid, mag, alpha)", "[298]"], COLOR_DRIVER, FILE_IBPM)
    y -= dy
    n05 = box(ax, 0, y, ["Build model + solver", "(branches on -model flag)", "[318-347]"],
              COLOR_DRIVER, FILE_IBPM, h=1.3)
    n05_y = y

    y -= dy + 1.3
    n06 = box(ax, 0, y, ["x = State(grid, ...);", "x.load(icFile)", "[361-369]"], COLOR_DRIVER,
              FILE_IBPM, h=1.3)

    y -= dy + 0.3
    n07 = box(ax, 0, y, ["geom.moveBodies(x.time)", "[394]"], COLOR_DRIVER, FILE_IBPM)
    n07_y = y
    draw_movebodies_subtree(ax, n07, 1, "394")

    y -= dy + 0.3
    n08 = box(ax, 0, y, ["model.init()", "[397]"], COLOR_DRIVER, FILE_IBPM)
    n08_y = y
    mi = box(ax, 1, y, ["NavierStokesModel.init()"], COLOR_MODEL, FILE_NSM, h=0.85, fontsize=7.4)
    call_arrow(ax, n08, mi, "397", color=COLOR_MODEL)
    reg = box(ax, 2, y, ["if not hasBeenInitialized:", "regularizer.update()"], COLOR_MODEL,
              FILE_NSM, h=0.95, fontsize=7.3)
    call_arrow(ax, mi, reg, "71", color=COLOR_MODEL)

    y -= dy + 0.3
    n09 = box(ax, 0, y, ["solver.load(...) or", "solver.init()+save(...)", "[400-405]"],
              COLOR_DRIVER, FILE_IBPM, h=1.3)
    n09_y = y
    sl = box(ax, 1, y + 0.85, ["solver.load(basename)"], COLOR_SOLVER, FILE_IBS, h=0.85, fontsize=7.3)
    call_arrow(ax, n09, sl, "400", color=COLOR_SOLVER, rad=0.12)
    si = box(ax, 1, y - 0.85, ["if load failed:", "solver.init(); solver.save(basename)"],
             COLOR_SOLVER, FILE_IBS, h=0.95, fontsize=7.1)
    call_arrow(ax, n09, si, "403-405", color=COLOR_SOLVER, rad=-0.12)
    ibs_loop = box(ax, 2, y, ["for i in range(nsteps):", "_solver[i].load/init/save(...)"],
                   COLOR_SOLVER, FILE_IBS, h=1.0, fontsize=7.2)
    call_arrow(ax, sl, ibs_loop, "122-128", color=COLOR_SOLVER, rad=0.12)
    call_arrow(ax, si, ibs_loop, "112-114,\n130-139", color=COLOR_SOLVER, rad=-0.12)
    chol = box(ax, 3, y, ["CholeskySolver.init():", "computeMatrixM()", "-> computeFactorization()",
               "ConjugateGradientSolver:", "init() is a no-op (base class)"], COLOR_LOWLEVEL,
               FILE_CHOL_CG, h=1.65, fontsize=6.9)
    call_arrow(ax, ibs_loop, chol, "110-122", color=COLOR_SOLVER)

    y -= dy + 0.4
    n10 = box(ax, 0, y, ["model.updateOperators(x.time)", "[411]"], COLOR_DRIVER, FILE_IBPM)
    n10_y = y
    draw_update_operators_subtree(ax, n10, 1, "411")

    y -= dy + 1.35
    n11 = box(ax, 0, y, ["model.refreshState(x)", "[412]"], COLOR_DRIVER, FILE_IBPM)
    n11_y = y
    draw_refresh_state_subtree(ax, n11, 1, "412")

    y -= dy + 0.4
    n12 = box(ax, 0, y, ["Setup Output objects", "(Tecplot/Restart/Force/Energy);",
              "logger.addOutput(...)", "[417-443]"], COLOR_DRIVER, FILE_IBPM, h=1.5)

    y -= dy + 0.15
    n13 = box(ax, 0, y, ["logger.init()", "[445]"], COLOR_DRIVER, FILE_IBPM)
    n13_y = y
    li = box(ax, 1, y, ["Logger.init()"], COLOR_IO, FILE_LOGGER, h=0.8, fontsize=7.4)
    call_arrow(ax, n13, li, "445", color=COLOR_IO)
    li_loop = box(ax, 2, y, ["for entry in _outputs:", "entry.output.init()"], COLOR_IO,
                  FILE_LOGGER, h=0.95, fontsize=7.3)
    call_arrow(ax, li, li_loop, "76-79", color=COLOR_IO)

    y -= dy + 0.15
    n14 = box(ax, 0, y, ["logger.doOutput(q_potential, x)", "(initial output)", "[446]"],
              COLOR_DRIVER, FILE_IBPM, h=1.15)
    n14_y = y
    draw_logger_dooutput_subtree(ax, n14, 1, "446")

    y -= dy + 1.4
    loop_top_y = y
    n15 = box(ax, 0, y, ["for i in 1..numSteps:", "[449]"], COLOR_DRIVER, FILE_IBPM,
              facecolor="#fdf2ef")

    y -= dy
    n16 = box(ax, 0, y, ["solver.advance(x)", "[452]"], COLOR_SOLVER, FILE_IBPM,
              facecolor="#f5eefb")

    # ================= BRANCH B: n16 "solver.advance(x)" (the big one) =================
    b1_y = n16[1]
    b1 = box(ax, 1, b1_y, ["IBSolver.advance(x)"], COLOR_SOLVER, FILE_IBS, h=0.9, fontsize=7.8)
    call_arrow(ax, n16, b1, "452", color=COLOR_SOLVER, lw=1.8, mutation_scale=15)

    b2 = box(ax, 2, b1_y, ["for i in range(nsteps):", "nonlinear = self.N(x)"], COLOR_SOLVER,
             FILE_IBS, h=1.0, fontsize=7.5)
    call_arrow(ax, b1, b2, "175,177", color=COLOR_SOLVER)

    b3 = box(ax, 3, b1_y + 1.5,
             ["N(x) polymorphic override:", "Nonlinear: Curl(Cross(q,ω))  [241]",
              "Linearized: Curl(Cross(q0,ω)+Cross(q,ω0))  [261]",
              "Adjoint: Lap(Curl(Cross(q0,q)))-Curl(Cross(q,ω0))  [282]",
              "Periodic: as Linearized, ω0->ω0periodic[k]  [311]",
              "SFD: as Nonlinear, - chi*(ω-ωhat)  [347]"],
             COLOR_NOTE, FILE_IBS, h=2.5, fontsize=6.7)
    call_arrow(ax, b2, b3, "177", color=COLOR_SOLVER, rad=0.12)

    b4 = box(ax, 3, b1_y - 1.5, ["advanceSubstep(x, nonlinear, i)"], COLOR_SOLVER, FILE_IBS,
             h=0.9, fontsize=7.6)
    call_arrow(ax, b2, b4, "180", color=COLOR_SOLVER, rad=-0.12)

    b5 = box(ax, 4, b1_y + 2.5, ["if isTimeDependent():", "model.updateOperators(t)"],
             COLOR_SOLVER, FILE_IBS, h=1.0, fontsize=7.3)
    call_arrow(ax, b4, b5, "187-188", color=COLOR_SOLVER, rad=0.32)
    draw_update_operators_subtree(ax, b5, 5, "188")

    b6 = box(ax, 4, b1_y + 0.55,
             ["a = Laplacian(ω)*coef", "+ coef*nonlinear (+bn*Nprev)"], COLOR_SOLVER,
             FILE_IBS, h=1.0, fontsize=7.3)
    call_arrow(ax, b4, b6, "191-203", color=COLOR_SOLVER, rad=0.14)

    b7 = box(ax, 4, b1_y - 0.75, ["b = model.getConstraints()"], COLOR_SOLVER, FILE_IBS,
             h=0.85, fontsize=7.3)
    call_arrow(ax, b4, b7, "206", color=COLOR_SOLVER, rad=-0.05)
    # placed clear ABOVE the wide ProjectionSolver.solve() breakdown (c_labels,
    # drawn below in column 5 across nearly all of branch B's height) so the
    # two column-5 subtrees don't collide
    gc = box(ax, 5, b1_y + 1.3, ["geometry.getVelocities()", "- toBoundary(baseFlow.getFlux())"],
             COLOR_GEOM, FILE_NSM, h=0.95, fontsize=7.0)
    call_arrow(ax, b7, gc, "105-107", color=COLOR_SOLVER, rad=-0.35)

    b8 = box(ax, 4, b1_y - 2.0, ["self._solver[i].solve(", "a, b, x.omega, x.f)"], COLOR_SOLVER,
             FILE_IBS, h=1.0, fontsize=7.3)
    call_arrow(ax, b4, b8, "209", color=COLOR_SOLVER, rad=-0.18)

    b9 = box(ax, 4, b1_y - 4.6, ["model.refreshState(x)"], COLOR_SOLVER, FILE_IBS, h=0.85,
             fontsize=7.3)
    call_arrow(ax, b4, b9, "212", color=COLOR_SOLVER, rad=-0.42)
    # skip column 5 entirely here (occupied top-to-bottom by the
    # ProjectionSolver.solve() breakdown below) -- start this subtree one
    # column further right instead of colliding with it
    draw_refresh_state_subtree(ax, b9, 6, "212")

    b_sfd = box(ax, 4, b1_y - 6.3, ["[SFDSolver only] also", "integrates filtered state ωhat"],
                COLOR_NOTE, FILE_IBS, h=1.0, fontsize=7.0)
    call_arrow(ax, b4, b_sfd, "354-386", color=COLOR_SOLVER, rad=-0.52)

    # ProjectionSolver.solve() breakdown -- one layer deeper than b8. All six
    # steps below are sequential lines inside ProjectionSolver.solve() itself
    # (projection_solver.py), even though several delegate out to
    # NavierStokesModel/HelmholtzSolver/CholeskySolver -- the "->" in each
    # box's text names that delegate, but the line cited on the arrow (and
    # the box's own file tag) is always projection_solver.py's.
    c_y = b8[1] + 2.05
    c_labels_lines = [
        (["Ainv(a,ω*) ->", "HelmholtzSolver.solve"], "128"),
        (["C(ω*,rhs) -> model.C:", "Poisson.solve+Curl+toBoundary"], "132"),
        (["Minv(rhs,f): Cholesky", "back-sub | CG iterate"], "134"),
        (["B(f,c) -> model.B:", "toFlux(f); Curl(q,ω)"], "138"),
        (["Ainv(c,c) ->", "HelmholtzSolver.solve"], "139"),
        (["ω.assign(ω* - c)"], "140"),
    ]
    for lines, lbl in c_labels_lines:
        cc = box(ax, 5, c_y, lines, COLOR_PROJ, FILE_PROJ, h=0.95, fontsize=6.7)
        call_arrow(ax, b8, cc, lbl, color=COLOR_PROJ, rad=0.0, fontsize=6.3, label_frac=0.4)
        c_y -= 1.35

    branchB_bottom = min(b_sfd[1] - b_sfd[3] / 2, c_y + 1.35 - 0.95 / 2) - 0.6

    b10 = box(ax, 2, branchB_bottom - 1.2, ["after all substeps:", "x.time += dt; x.timestep += 1"],
              COLOR_SOLVER, FILE_IBS, h=1.0, fontsize=7.4)
    # SAME column as b2 (both are statements in IBSolver.advance's own body,
    # b10 sequentially after the substep loop) -- straight DOWN, never left.
    seq_arrow(ax, b2, b10, COLOR_SOLVER, lw=1.0, ls="dashed")
    loop_arrow(ax, b4, b2, COLOR_SOLVER, x_offset=1.8, label="repeat per\nsubstep i", fontsize=6.5)

    connector_y = b1_y  # for the a8->b_hdr-style connector below (kept implicit via n16->b1 arrow)

    # ===================== back in column 0, below branch B entirely =====================
    y = branchB_bottom - dy - 1.0
    n17 = box(ax, 0, y, ["x.computeNetForce()", "[453]"], COLOR_DRIVER, FILE_IBPM)

    y -= dy + 0.15
    n18 = box(ax, 0, y, ["logger.doOutput(q_potential, x)", "[462]"], COLOR_DRIVER, FILE_IBPM,
              h=1.15)
    n18_y = y
    draw_logger_dooutput_subtree(ax, n18, 1, "462")
    loop_bottom_y = y

    y -= dy + 1.4
    n19 = box(ax, 0, y, ["logger.cleanup()", "[480]"], COLOR_DRIVER, FILE_IBPM)
    lc = box(ax, 1, y, ["Logger.cleanup()"], COLOR_IO, FILE_LOGGER, h=0.8, fontsize=7.4)
    call_arrow(ax, n19, lc, "480", color=COLOR_IO)
    lc_loop = box(ax, 2, y, ["for entry in _outputs:", "entry.output.cleanup()"], COLOR_IO,
                  FILE_LOGGER, h=0.95, fontsize=7.3)
    call_arrow(ax, lc, lc_loop, "84-88", color=COLOR_IO)

    # ===================== column-0 sequencing arrows (straight down only) =====================
    for p, q in [(n01, n02), (n02, n03), (n03, n04), (n04, n05), (n05, n06),
                 (n06, n07), (n07, n08), (n08, n09), (n09, n10), (n10, n11),
                 (n11, n12), (n12, n13), (n13, n14), (n14, n15), (n15, n16),
                 (n16, n17), (n17, n18), (n18, n19)]:
        seq_arrow(ax, p, q, COLOR_DRIVER)

    # ===================== BRANCH A: n05 "build model + solver" =====================
    a1 = box(ax, 1, n05_y + 0.9, ["NavierStokesModel", "(grid, geom, Re[, q_potential])"],
             COLOR_MODEL, FILE_NSM, h=1.0, fontsize=7.6)
    a2 = box(ax, 1, n05_y - 0.9, ["<Solver>(grid, model, dt,", "scheme, ...) -- Nonlinear/",
             "Linearized/Adjoint/Periodic/SFD"], COLOR_SOLVER, FILE_IBS, h=1.25, fontsize=7.3)
    call_arrow(ax, n05, a1, "318/324/\n330/340/346")
    call_arrow(ax, n05, a2, "319/325/\n331/341/347")

    a3 = box(ax, 2, a2[1], ["IBSolver.__init__", "-> self.createAllSolvers()"], COLOR_SOLVER,
             FILE_IBS, h=1.0, fontsize=7.4)
    call_arrow(ax, a2, a3, "84-104")

    a4 = box(ax, 3, a2[1], ["createAllSolvers(): loop substeps", "-> createSolver(beta)"],
             COLOR_SOLVER, FILE_IBS, h=1.0, fontsize=7.4)
    call_arrow(ax, a3, a4, "104")

    a5 = box(ax, 4, a2[1], ["CholeskySolver (stationary)", "ConjugateGradientSolver (moving)"],
             COLOR_LOWLEVEL, FILE_CHOL_CG, h=1.0, fontsize=7.4)
    call_arrow(ax, a4, a5, "141-165")

    # ===================== loop-back arrows (the only intentionally backward arrows) ============
    loop_arrow(ax, n18, n15, COLOR_DRIVER, x_offset=-2.2, label="repeat\nnumSteps\ntimes")

    ax.set_title(
        "ibpm_py -- execution flow of a simulation run (py/ibpm.py: main())\n"
        "Column 0 = ibpm.py's own statements, TOP TO BOTTOM IN ACTUAL SOURCE ORDER. Every "
        "arrow is labeled with the file:line of the CALL SITE it leaves from; call depth "
        "increases strictly left to right; the only backward (dashed) arrows are the two real "
        "loops.", fontsize=13, pad=10,
    )

    all_y = [b[1] - b[3] / 2 for b in [n19, lc_loop]] + [c_y]
    ax.set_xlim(-2.5, COL_X[-1] + COL_W[-1] / 2 + 2)
    ax.set_ylim(min(all_y) - 1.5, 61.5)
    ax.axis("off")
    fig.subplots_adjust(left=0.005, right=0.995, top=0.965, bottom=0.005)

    out_png = os.path.join(OUT_DIR, "main_execution_flowchart.png")
    fig.savefig(out_png, dpi=150)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
