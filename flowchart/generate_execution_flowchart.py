"""
generate_execution_flowchart.py

Hand-traced "what calls what" flowchart for actually running a simulation
with py/ibpm.py's main() -- i.e. the numerical algorithm, not just the
import graph (see generate_module_graph.py for the mechanical version).

Layout: a strict left-to-right layered call graph. Column 0 is every
statement of interest in ibpm.py's main(), top to bottom, in execution
order -- so "where does ibpm.py call into X" is always answered by scanning
straight down column 0. Every rightward arrow is labeled with the file:line
of the CALL SITE (i.e. the line in the box the arrow leaves, not the line
the callee is defined on) -- so hovering over any arrow answers "when is
this called". Call depth increases left to right: column 1 is what column-0
code calls directly, column 2 is what column-1 code calls, and so on.

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
COLOR_NOTE = "#d35400"        # polymorphic dispatch / SFD-only notes

# column x-centers (call depth 0 = ibpm.py, increasing = deeper calls)
COL_X = [3.0, 10.5, 18.3, 26.7, 35.2, 43.2]
COL_W = [6.0, 6.6, 7.2, 7.8, 6.6, 6.2]


def box(ax, col, y, lines, edgecolor, file, h=1.0, fontsize=8.0, facecolor="white",
        linewidth=1.4, zorder=2, style="round,pad=0.02,rounding_size=0.07"):
    """Draw one flowchart box. `file` is the py/ source file this box's code
    belongs to; it's rendered as a small tag centered directly above the
    box's top edge (drawn outside the box, so it never eats into body text
    or collides with an incoming arrow, which always lands at top-center)."""
    x, w = COL_X[col], COL_W[col]
    b = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=style,
                         linewidth=linewidth, edgecolor=edgecolor,
                         facecolor=facecolor, zorder=zorder)
    ax.add_patch(b)
    ax.text(x, y, "\n".join(lines), ha="center", va="center", fontsize=fontsize,
             color="#222222", zorder=zorder + 1, linespacing=1.25)
    # file tag: upper-left, outside the box, so it never sits in the path of
    # a straight vertical incoming arrow (which lands at top-center)
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
    """Arrow for an actual function call: src (left) -> dst (right),
    labeled with the file:line of the call site (a line inside src's own
    code) -- NOT the line dst is defined on."""
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
    """Plain top-to-bottom sequencing arrow within the same column."""
    a = FancyArrowPatch(bottom(src), top(dst), arrowstyle="-|>",
                         mutation_scale=11, color=color, linewidth=lw,
                         linestyle=ls, zorder=1)
    ax.add_patch(a)


def loop_arrow(ax, src, dst, color, x_offset, label, fontsize=7.0, ls="dashed"):
    """Dashed feedback loop routed through a vertical line to the side."""
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


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(23, 18))

    # ===================== COLUMN 0: ibpm.py main(), in order =====================
    IBPM = "ibpm.py"
    y = 45.0
    dy = 1.75
    n01 = box(ax, 0, y, ["ParmParser(argc, argv)", "[166]"], COLOR_DRIVER, IBPM)
    y -= dy
    n02 = box(ax, 0, y, ["Grid(nx, ny, ngrid, ...)", "[280]"], COLOR_DRIVER, IBPM)
    y -= dy
    n03 = box(ax, 0, y, ["geom.load(geomFile)", "[285]"], COLOR_DRIVER, IBPM)
    y -= dy
    n04 = box(ax, 0, y, ["BaseFlow(grid, mag, alpha)", "[298]"], COLOR_DRIVER, IBPM)
    y -= dy
    n05 = box(ax, 0, y, ["Build model + solver", "(branches on -model flag)", "[318-347]"],
              COLOR_DRIVER, IBPM, h=1.3)
    y -= dy - 0.1
    n06 = box(ax, 0, y, ["model.init()", "[397]"], COLOR_DRIVER, IBPM)
    y -= dy
    n07 = box(ax, 0, y, ["solver.load(...) or", "solver.init()+save(...)", "[400-405]"],
              COLOR_DRIVER, IBPM, h=1.3)
    y -= dy + 0.05
    n08 = box(ax, 0, y, ["model.updateOperators(t);", "model.refreshState(x)", "[411-412]"],
              COLOR_DRIVER, IBPM, h=1.3)
    y -= dy + 0.05
    n09 = box(ax, 0, y, ["x = State(...); x.load(icFile)", "[363-369]"], COLOR_DRIVER, IBPM)
    y -= dy
    n10 = box(ax, 0, y, ["geom.moveBodies(x.time)", "[394]"], COLOR_DRIVER, IBPM)
    y -= dy
    n11 = box(ax, 0, y, ["Register outputs w/ Logger;", "logger.init()/doOutput(...)", "[417-446]"],
              COLOR_DRIVER, IBPM, h=1.3)
    y -= dy + 0.35
    n12 = box(ax, 0, y, ["for i in 1..numSteps:", "[449]"], COLOR_DRIVER, IBPM, facecolor="#fdf2ef")
    y -= dy
    n13 = box(ax, 0, y, ["solver.advance(x)", "[452]"], COLOR_SOLVER, IBPM, facecolor="#f5eefb")
    y -= dy
    n14 = box(ax, 0, y, ["x.computeNetForce()", "[453]"], COLOR_DRIVER, IBPM)
    y -= dy
    n15 = box(ax, 0, y, ["logger.doOutput(q_potential, x)", "[462]"], COLOR_DRIVER, IBPM)
    n15_loop_bottom = y
    y -= dy + 0.4
    n16 = box(ax, 0, y, ["logger.cleanup()", "[480]"], COLOR_DRIVER, IBPM)

    for p, q in [(n01, n02), (n02, n03), (n03, n04), (n04, n05), (n05, n06),
                 (n06, n07), (n07, n08), (n08, n09), (n09, n10), (n10, n11),
                 (n11, n12), (n12, n13), (n13, n14), (n14, n15)]:
        seq_arrow(ax, p, q, COLOR_DRIVER)
    seq_arrow(ax, n15, n16, COLOR_DRIVER)
    loop_arrow(ax, n15, n12, COLOR_DRIVER, x_offset=-1.6, label="repeat\nnumSteps\ntimes")

    # ===================== BRANCH A: n05 "build model + solver" =====================
    NSM = "navier_stokes_model.py"
    IBS = "ib_solver.py"
    PROJ = "projection_solver.py"
    CHOL_CG = "cholesky_solver.py | conjugate_gradient_solver.py"

    a1 = box(ax, 1, n05[1] + 0.9, ["NavierStokesModel", "(grid, geom, Re[, q_potential])"],
             COLOR_MODEL, NSM, h=1.0, fontsize=7.6)
    a2 = box(ax, 1, n05[1] - 0.9, ["<Solver>(grid, model, dt,", "scheme, ...) -- Nonlinear/", "Linearized/Adjoint/Periodic/SFD"],
             COLOR_SOLVER, IBS, h=1.25, fontsize=7.3)
    call_arrow(ax, n05, a1, "318/324/\n330/340/346")
    call_arrow(ax, n05, a2, "319/325/\n331/341/347")

    a3 = box(ax, 2, a2[1], ["IBSolver.__init__", "-> self.createAllSolvers()"], COLOR_SOLVER,
             IBS, h=1.0, fontsize=7.4)
    call_arrow(ax, a2, a3, "84-104")

    a4 = box(ax, 3, a2[1], ["createAllSolvers(): loop substeps", "-> createSolver(beta)"], COLOR_SOLVER,
             IBS, h=1.0, fontsize=7.4)
    call_arrow(ax, a3, a4, "104")

    a5 = box(ax, 4, a2[1], ["CholeskySolver (stationary)", "ConjugateGradientSolver (moving)"],
             COLOR_LOWLEVEL, CHOL_CG, h=1.0, fontsize=7.4)
    call_arrow(ax, a4, a5, "141-165")

    # ===================== BRANCH B: n13 "solver.advance(x)" =====================
    b1_y = n13[1]
    b1 = box(ax, 1, b1_y, ["IBSolver.advance(x)"], COLOR_SOLVER, IBS, h=0.9, fontsize=7.8)
    call_arrow(ax, n13, b1, "452", color=COLOR_SOLVER, lw=1.8, mutation_scale=15)

    b2 = box(ax, 2, b1_y, ["for i in range(nsteps):", "nonlinear = self.N(x)"], COLOR_SOLVER,
             IBS, h=1.0, fontsize=7.5)
    call_arrow(ax, b1, b2, "175,177", color=COLOR_SOLVER)

    b3 = box(ax, 3, b1_y + 1.45,
             ["N(x) polymorphic override:", "Nonlinear: Curl(Cross(q,ω))  [241]",
              "Linearized: Curl(Cross(q0,ω)+Cross(q,ω0))  [261]",
              "Adjoint: Lap(Curl(Cross(q0,q)))-Curl(Cross(q,ω0))  [282]",
              "Periodic: as Linearized, ω0->ω0periodic[k]  [311]",
              "SFD: as Nonlinear, - chi*(ω-ωhat)  [347]"],
             COLOR_NOTE, IBS, h=2.5, fontsize=6.7)
    call_arrow(ax, b2, b3, "177", color=COLOR_SOLVER, rad=0.12)

    b4 = box(ax, 3, b1_y - 1.45, ["advanceSubstep(x, nonlinear, i)"], COLOR_SOLVER, IBS, h=0.9, fontsize=7.6)
    call_arrow(ax, b2, b4, "180", color=COLOR_SOLVER, rad=-0.12)

    b5 = box(ax, 4, b1_y + 1.65, ["if isTimeDependent():", "model.updateOperators(t)"], COLOR_SOLVER,
             IBS, h=1.0, fontsize=7.3)
    call_arrow(ax, b4, b5, "187-188", color=COLOR_SOLVER, rad=0.28)

    b6 = box(ax, 4, b1_y + 0.2,
             ["a = Laplacian(ω)*coef", "+ coef*nonlinear (+bn*Nprev)"], COLOR_SOLVER,
             IBS, h=1.0, fontsize=7.3)
    call_arrow(ax, b4, b6, "191-203", color=COLOR_SOLVER, rad=0.1)

    b7 = box(ax, 4, b1_y - 1.25, ["b = model.getConstraints()"], COLOR_SOLVER, IBS, h=0.85, fontsize=7.3)
    call_arrow(ax, b4, b7, "206", color=COLOR_SOLVER, rad=-0.1)

    b8 = box(ax, 4, b1_y - 2.5, ["self._solver[i].solve(", "a, b, x.omega, x.f)"], COLOR_SOLVER,
             IBS, h=1.0, fontsize=7.3)
    call_arrow(ax, b4, b8, "209", color=COLOR_SOLVER, rad=-0.28)

    b9 = box(ax, 4, b1_y - 3.9, ["model.refreshState(x)"], COLOR_SOLVER, IBS, h=0.85, fontsize=7.3)
    call_arrow(ax, b4, b9, "212", color=COLOR_SOLVER, rad=-0.36)

    b_sfd = box(ax, 4, b1_y - 5.1, ["[SFDSolver only] also", "integrates filtered state ωhat"],
                COLOR_NOTE, IBS, h=1.0, fontsize=7.0)
    call_arrow(ax, b4, b_sfd, "354-386", color=COLOR_SOLVER, rad=-0.42)

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
        cc = box(ax, 5, c_y, lines, COLOR_PROJ, PROJ, h=0.95, fontsize=6.7)
        call_arrow(ax, b8, cc, lbl, color=COLOR_PROJ, rad=0.0, fontsize=6.3, label_frac=0.4)
        c_y -= 1.35

    b10 = box(ax, 2, b1_y - 7.2, ["after all substeps:", "x.time += dt; x.timestep += 1"],
              COLOR_SOLVER, IBS, h=1.0, fontsize=7.4)
    call_arrow(ax, b2, b10, "182-183\n(after loop)", color=COLOR_SOLVER, rad=0.0, lw=1.0)
    loop_arrow(ax, b4, b2, COLOR_SOLVER, x_offset=1.55, label="repeat per\nsubstep i", fontsize=6.5)

    ax.set_title(
        "ibpm_py -- execution flow of a simulation run (py/ibpm.py: main())\n"
        "Column 0 = ibpm.py's own statements, top to bottom. Each arrow is labeled with the "
        "file:line of the CALL SITE it leaves from; call depth increases left to right.",
        fontsize=12.5, pad=8,
    )

    ax.set_xlim(-2.5, 47.5)
    ax.set_ylim(min(b_sfd[1], b10[1]) - 1.3, 46.2)
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.01)

    out_png = os.path.join(OUT_DIR, "main_execution_flowchart.png")
    fig.savefig(out_png, dpi=190)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
