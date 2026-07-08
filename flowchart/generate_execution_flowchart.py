"""
generate_execution_flowchart.py

Hand-traced "what calls what" flowchart for actually running a simulation
with py/ibpm.py's main() -- i.e. the numerical algorithm, not just the
import graph (see generate_module_graph.py for the mechanical version).

Every box cites the file:line where that call happens in THIS repo, as of
the commit this was generated against. Nothing here is inferred: each edge
was found by opening the cited file and reading the call. See
flowchart/README.md for how to re-verify each box yourself, and
flowchart/execution_flow_refs.csv for the same citations in a table.

Usage:
    python3 flowchart/generate_execution_flowchart.py

Output:
    flowchart/output/main_execution_flowchart.png
"""

from __future__ import annotations

import csv
import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

COLOR_DRIVER = "#c0392b"      # ibpm.py:main()
COLOR_SOLVER = "#8e44ad"      # ib_solver.py (IBSolver.advance/advanceSubstep)
COLOR_PROJ = "#2980b9"        # projection_solver.py
COLOR_MODEL = "#16a085"       # navier_stokes_model.py
COLOR_LOWLEVEL = "#7f8c8d"    # elliptic_solver.py / regularizer.py / vector_operations.py
COLOR_NOTE = "#d35400"        # side-notes / polymorphic dispatch


def draw_box(ax, x, y, w, h, lines, edgecolor, fontsize=8.2, facecolor="white",
             linewidth=1.5, zorder=2, style="round,pad=0.02,rounding_size=0.08"):
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=style, linewidth=linewidth,
        edgecolor=edgecolor, facecolor=facecolor, zorder=zorder,
    )
    ax.add_patch(box)
    text = "\n".join(lines)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
             color="#222222", zorder=zorder + 1, linespacing=1.35)
    return (x, y, w, h)


def arrow(ax, p1, p2, color="#333333", rad=0.0, lw=1.4, style="-|>",
          ls="solid", mutation_scale=13, zorder=1):
    a = FancyArrowPatch(p1, p2, connectionstyle=f"arc3,rad={rad}",
                          arrowstyle=style, mutation_scale=mutation_scale,
                          color=color, linewidth=lw, linestyle=ls, zorder=zorder)
    ax.add_patch(a)


def bottom(box):
    x, y, w, h = box
    return (x, y - h / 2)


def top(box):
    x, y, w, h = box
    return (x, y + h / 2)


def left(box):
    x, y, w, h = box
    return (x - w / 2, y)


def right(box):
    x, y, w, h = box
    return (x + w / 2, y)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(30, 34))

    # ================= LEFT COLUMN: ibpm.py main() driver =================
    LX = 5.0
    LW = 9.4
    y = 33.0

    a1 = draw_box(ax, LX, y, LW, 1.5,
                  ["1. Parse CLI arguments", "ParmParser(argc, argv)  [ibpm.py:166]"],
                  COLOR_DRIVER)
    y -= 2.1
    a2 = draw_box(ax, LX, y, LW, 1.9,
                  ["2. Build Grid / Geometry / BaseFlow",
                   "Grid(...)  [ibpm.py:280]",
                   "Geometry().load(geomFile)  [ibpm.py:283-288]",
                   "BaseFlow(grid, magnitude, alpha)  [ibpm.py:298]"],
                  COLOR_DRIVER)
    y -= 2.4
    a3 = draw_box(ax, LX, y, LW, 2.6,
                  ["3. Build model + solver, per -model flag",
                   "NavierStokesModel(grid, geom, Re[, q_potential])",
                   "  [ibpm.py:318/324/330/340/346]",
                   "solver = Nonlinear-/Linearized-/Adjoint-/",
                   "  LinearizedPeriodic-/SFDSolver(grid, model, dt, scheme,...)",
                   "  [ibpm.py:319/325/331/341/347]"],
                  COLOR_DRIVER)
    y -= 2.55
    a3note = draw_box(ax, LX, y, LW, 1.9,
                       ["IBSolver.__init__ eagerly calls createAllSolvers()",
                        "-> createSolver(beta) per RK/AB2 substep, picks:",
                        "CholeskySolver  (stationary body)  [ib_solver.py:153-165]",
                        "ConjugateGradientSolver  (moving body)"],
                       COLOR_NOTE, fontsize=7.6)
    y -= 2.35
    a4 = draw_box(ax, LX, y, LW, 2.5,
                  ["4. Initialize model + solver",
                   "model.init() -> regularizer.update()  [ibpm.py:397]",
                   "solver.load(...) or solver.init()+save(...)",
                   "  -> ProjectionSolver.init() per substep  [ibpm.py:400-405]",
                   "model.updateOperators(t); model.refreshState(x)  [ibpm.py:411-412]"],
                  COLOR_DRIVER)
    y -= 2.5
    a5 = draw_box(ax, LX, y, LW, 1.9,
                  ["5. Load initial condition",
                   "x = State(grid, geom.getNumPoints())  [ibpm.py:363]",
                   "x.load(icFile)  [ibpm.py:369]",
                   "geom.moveBodies(x.time)  [ibpm.py:394]"],
                  COLOR_DRIVER)
    y -= 2.4
    a6 = draw_box(ax, LX, y, LW, 2.3,
                  ["6. Register outputs with Logger",
                   "OutputTecplot/Restart/Force/Energy(...)  [ibpm.py:417-426]",
                   "logger.addOutput(output, everyNsteps)  [ibpm.py:432-443]",
                   "logger.init(); logger.doOutput(q_potential, x)  [ibpm.py:445-446]"],
                  COLOR_DRIVER)

    y -= 2.55
    loop_top_y = y
    a7 = draw_box(ax, LX, y, LW, 1.1,
                  ["7. for i in 1..numSteps:   [ibpm.py:449]"],
                  COLOR_DRIVER, facecolor="#fdf2ef")
    y -= 1.75
    a8 = draw_box(ax, LX, y, LW, 1.1,
                  ["solver.advance(x)   [ibpm.py:452]", "-- see right panel -->"],
                  COLOR_SOLVER, facecolor="#f5eefb")
    y -= 1.9
    a9 = draw_box(ax, LX, y, LW, 2.15,
                  ["x.computeNetForce()  [ibpm.py:453]",
                   "logger.doOutput(q_potential, x)  [ibpm.py:462]",
                   "-> Output.doOutput(...) for each registered",
                   "   OutputTecplot / OutputRestart / OutputForce / OutputEnergy"],
                  COLOR_DRIVER, facecolor="#fdf2ef")
    loop_bottom_y = y

    y -= 2.35
    a10 = draw_box(ax, LX, y, LW, 1.1,
                   ["8. logger.cleanup()   [ibpm.py:480]"],
                   COLOR_DRIVER)

    for p, q in [(a1, a2), (a2, a3), (a3, a3note), (a3note, a4), (a4, a5),
                 (a5, a6), (a6, a7), (a7, a8), (a8, a9)]:
        arrow(ax, bottom(p), top(q), color=COLOR_DRIVER)
    arrow(ax, bottom(a9), top(a10), color=COLOR_DRIVER)
    # loop-back arrow from a9 to a7
    loop_x = LX - LW / 2 - 0.9
    arrow(ax, (loop_x, loop_bottom_y), (loop_x, loop_top_y), color=COLOR_DRIVER,
          rad=0.0, ls="dashed")
    arrow(ax, left(a9), (loop_x, loop_bottom_y), color=COLOR_DRIVER, ls="dashed", mutation_scale=1)
    arrow(ax, (loop_x, loop_top_y), left(a7), color=COLOR_DRIVER, ls="dashed")
    ax.text(loop_x - 0.35, (loop_top_y + loop_bottom_y) / 2, "repeat\nnumSteps\ntimes",
             ha="center", va="center", fontsize=7.5, color=COLOR_DRIVER, rotation=90)

    # ================= RIGHT PANEL: IBSolver.advance(x) detail =================
    RX = 21.5
    RW = 17.6
    ry = 30.0

    panel_top = 31.4
    b_hdr = draw_box(ax, RX, ry, RW, 1.35,
                      ["IBSolver.advance(x)  --  ib_solver.py:171",
                       "one call per outer-loop iteration (step 7-8, left)"],
                      COLOR_SOLVER, facecolor="#f5eefb", fontsize=9)
    ry -= 1.95
    b1 = draw_box(ax, RX, ry, RW, 1.0,
                  ["for i in range(scheme.nsteps()):   [ib_solver.py:175]",
                   "(nsteps = 1 for Euler, 2 for AB2, 3 for RK3/RK3b)"],
                  COLOR_SOLVER, fontsize=7.6)
    ry -= 1.85
    b2 = draw_box(ax, RX, ry, RW, 1.1,
                  ["nonlinear = self.N(x)   [ib_solver.py:177]",
                   "polymorphic -- see side-note below"],
                  COLOR_SOLVER)
    ry -= 2.3
    b2note = draw_box(ax, RX, ry, RW, 2.9,
                       ["N(x) override per solver subclass:",
                        "NonlinearIBSolver:  Curl(CrossProduct(x.q,x.omega))  [241]",
                        "LinearizedIBSolver: Curl(CrossProduct(x0.q,x.omega)",
                        "                    + CrossProduct(x.q,x0.omega))  [261]",
                        "AdjointIBSolver:    Laplacian(Curl(CrossProduct(x0.q,x.q)))",
                        "                    - Curl(CrossProduct(x.q,x0.omega))  [282]",
                        "LinearizedPeriodicIBSolver: as Linearized, x0 -> x0periodic[k]  [311]",
                        "SFDSolver: as Nonlinear, minus chi*(x.omega - xhat.omega)  [347]"],
                       COLOR_NOTE, fontsize=7.3)
    ry -= 3.15
    b3 = draw_box(ax, RX, ry, RW, 1.0,
                  ["advanceSubstep(x, nonlinear, i)   [ib_solver.py:185]"],
                  COLOR_SOLVER)
    ry -= 1.85
    b3a = draw_box(ax, RX, ry, RW, 1.55,
                   ["if model.isTimeDependent():",
                    "  model.updateOperators(t)  [ib_solver.py:187-188]",
                    "  -> geometry.moveBodies(t), regularizer.update()  [navier_stokes_model.py:110-116]"],
                   COLOR_MODEL, fontsize=7.4)
    ry -= 1.95
    b3b = draw_box(ax, RX, ry, RW, 1.4,
                   ["a = Laplacian(x.omega)*coeffs + coeffs*nonlinear",
                    "    (+ bn(i)*Nprev when bn(i) != 0, i.e. AB2)  [ib_solver.py:191-203]",
                    "b = model.getConstraints()  [ib_solver.py:206]  -> geometry.getVelocities()",
                    "    - regularizer.toBoundary(baseFlow.getFlux())  [navier_stokes_model.py:99-108]"],
                   COLOR_SOLVER, fontsize=7.4)

    ry -= 2.05
    b4 = draw_box(ax, RX, ry, RW, 1.15,
                  ["self._solver[i].solve(a, b, x.omega, x.f)",
                   "-> ProjectionSolver.solve(...)   [projection_solver.py:108]"],
                  COLOR_PROJ)

    # nested sub-steps of ProjectionSolver.solve, indented
    ry -= 1.7
    sub_w = RW - 1.6
    b4a = draw_box(ax, RX, ry, sub_w, 1.0,
                   ["Ainv(a, omegaStar) -> HelmholtzSolver.solve(a, omegaStar)",
                    "[projection_solver.py:128,144-146]"],
                   COLOR_LOWLEVEL, fontsize=7.3)
    ry -= 1.55
    b4b = draw_box(ax, RX, ry, sub_w, 1.55,
                   ["C(omegaStar, rhs) -> model.C(omega, f)  [projection_solver.py:132]",
                    "-> computeFluxWithoutBaseFlow: vorticityToStreamfunction",
                    "   (PoissonSolver.solve) then Curl  [navier_stokes_model.py:132-148]",
                    "-> regularizer.toBoundary(q)"],
                   COLOR_MODEL, fontsize=7.3)
    ry -= 1.9
    b4c = draw_box(ax, RX, ry, sub_w, 1.55,
                   ["Minv(rhs, f)  [projection_solver.py:134]",
                    "CholeskySolver.Minv: back-substitute precomputed factorization",
                    "  of M = C.Ainv.B  [cholesky_solver.py:270]",
                    "ConjugateGradientSolver.Minv: iterate, calling M()->B,Ainv,C each step"],
                   COLOR_LOWLEVEL, fontsize=7.3)
    ry -= 1.9
    b4d = draw_box(ax, RX, ry, sub_w, 1.3,
                   ["B(f, c) -> model.B(f, omega)  [projection_solver.py:138]",
                    "-> regularizer.toFlux(f); Curl(q, omega)  [navier_stokes_model.py:122-130]"],
                   COLOR_MODEL, fontsize=7.3)
    ry -= 1.55
    b4e = draw_box(ax, RX, ry, sub_w, 0.85,
                   ["Ainv(c, c) -> HelmholtzSolver.solve(c, c)  [projection_solver.py:139]"],
                   COLOR_LOWLEVEL, fontsize=7.3)
    ry -= 1.35
    b4f = draw_box(ax, RX, ry, sub_w, 0.85,
                   ["omega.assign(omegaStar - c)  [projection_solver.py:140]"],
                   COLOR_PROJ, fontsize=7.3)

    ry -= 1.7
    b5 = draw_box(ax, RX, ry, RW, 1.5,
                  ["model.refreshState(x)  [ib_solver.py:212]",
                   "-> model.computeFlux(x.omega, x.q)  [navier_stokes_model.py:158-166]",
                   "   = computeFluxWithoutBaseFlow(omega, q) + baseFlow.getFlux()"],
                  COLOR_MODEL, fontsize=7.6)
    ry -= 1.85
    b6 = draw_box(ax, RX, ry, RW, 1.3,
                  ["[SFDSolver only] advanceSubstep also integrates the",
                   "filtered state _xhat, used by N() above  [ib_solver.py:354-386]"],
                  COLOR_NOTE, fontsize=7.4, style="round,pad=0.02,rounding_size=0.08")
    ry -= 1.75
    b7 = draw_box(ax, RX, ry, RW, 1.0,
                  ["after all substeps: x.time += dt; x.timestep += 1   [ib_solver.py:182-183]"],
                  COLOR_SOLVER, fontsize=7.8)

    for p, q in [(b_hdr, b1), (b1, b2), (b2, b2note), (b2note, b3), (b3, b3a),
                 (b3a, b3b), (b3b, b4), (b4, b4a), (b4a, b4b), (b4b, b4c),
                 (b4c, b4d), (b4d, b4e), (b4e, b4f), (b4f, b5), (b5, b6), (b6, b7)]:
        arrow(ax, bottom(p), top(q), color=COLOR_SOLVER if p in (b_hdr, b1, b2, b2note, b3, b6) else "#555555")

    # dashed loop arrow: after b7, back up to b1 for the next substep
    loop_x2 = RX + RW / 2 + 0.9
    arrow(ax, right(b7), (loop_x2, ry), color=COLOR_SOLVER, mutation_scale=1)
    arrow(ax, (loop_x2, ry), (loop_x2, top(b1)[1]), color=COLOR_SOLVER, ls="dashed")
    arrow(ax, (loop_x2, top(b1)[1]), right(b1), color=COLOR_SOLVER, ls="dashed")
    ax.text(loop_x2 + 0.35, (ry + top(b1)[1]) / 2, "repeat per\nsubstep i",
             ha="center", va="center", fontsize=7.5, color=COLOR_SOLVER, rotation=90)

    # connector from left panel (a8) into right panel (b_hdr)
    arrow(ax, right(a8), left(b_hdr), color=COLOR_SOLVER, rad=-0.15, lw=2.0, mutation_scale=18)

    # panel border for right side
    ax.add_patch(Rectangle((RX - RW / 2 - 1.3, ry - 0.8), RW + 2.6, panel_top - (ry - 0.8),
                             fill=False, linestyle="dashed", edgecolor="#999999", linewidth=1.2, zorder=0))

    ax.set_title(
        "ibpm_py -- hand-traced execution flow of a simulation run (py/ibpm.py: main())\n"
        "Every box cites its exact file:line; cross-check against the source (see flowchart/README.md)",
        fontsize=14, pad=14,
    )

    ax.set_xlim(-2, 33)
    ax.set_ylim(ry - 2, 34.2)
    ax.axis("off")

    out_png = os.path.join(OUT_DIR, "main_execution_flowchart.png")
    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
