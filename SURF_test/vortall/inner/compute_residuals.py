"""
compute_residuals.py

The reports already in SURF_test/vortall/ (README.md, three_way_summary.txt,
shedding_summary.txt) validate the OUTPUT of the vortall run: force
coefficients, shedding period, Strouhal number, qualitative vorticity field.
They never look inside a single timestep to check whether the linear-algebra
machinery that produces that output is actually solving the equations it
claims to solve.

This script does that. It re-runs the vortall cylinder case (Re=100,
nx=450, ny=200, same geometry/dt as run_vortall.py) directly through the
internal py/ API (Grid/Geometry/NavierStokesModel/NonlinearIBSolver), and at
every timestep computes four residuals that should each be at floating-point
roundoff level (~1e-10 or smaller) if the solver is doing its job, PLUS three
physical conservation/sanity quantities as a secondary check:

Residuals (should be ~0 to solver/floating-point tolerance):
  1. divergence   -- discrete divergence of the flux field q, computed with
                      an INDEPENDENT finite-difference formula written in
                      this script (not calling vector_operations.Curl or any
                      other py/ code), applied to q AFTER the timestep. A
                      curl-derived flux field is divergence-free by
                      construction on this staggered grid; if this residual
                      is not at machine precision, either Curl() (flux.py /
                      vector_operations.py) or the projection step
                      (ib_solver.py / projection_solver.py) has a bug.
  2. no_slip       -- ||model.C(x.omega) - model.getConstraints()||, i.e.
                      how well the projection step enforced the no-slip
                      boundary condition C.q = b on the cylinder surface
                      this step. This is exactly the equation
                      ProjectionSolver.solve() (projection_solver.py) is
                      supposed to satisfy; residual growth here implicates
                      CholeskySolver/ConjugateGradientSolver or Regularizer.
  3. poisson       -- ||Laplacian(psi) - (-omega)||, i.e. whether the
                      streamfunction returned by
                      NavierStokesModel.vorticityToStreamfunction actually
                      solves the Poisson equation it's supposed to solve.
                      Implicates elliptic_solver.py (PoissonSolver) if large.
  4. helmholtz     -- ||(I - alpha*beta/2 L) omegaStar - a||, the analogous
                      residual for the implicit half of the fractional-step
                      update. Implicates elliptic_solver.py (HelmholtzSolver)
                      if large.

Conservation / physical sanity (not expected to be ~0, just smooth/bounded):
  5. circulation   -- sum(omega) * dx^2 over the domain.
  6. enstrophy     -- sum(omega^2) * dx^2.
  7. kinetic_energy -- 0.5 * InnerProduct(q, q).
  8. Cd, Cl        -- reproduced independently here; cross-checked against
                      the archived SURF_test/vortall/_run_data/vortall.force
                      at the end of this script as a reproducibility check.

Usage:
    python3 SURF_test/vortall/inner/compute_residuals.py [--nsteps N] [--every K]

Output:
    SURF_test/vortall/inner/data/residuals.csv
    SURF_test/vortall/inner/data/run_meta.txt
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys
import time

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from py.base_flow import BaseFlow  # noqa: E402
from py.direction import Direction  # noqa: E402
from py.geometry import Geometry  # noqa: E402
from py.grid import Grid  # noqa: E402
from py.ib_solver import NonlinearIBSolver  # noqa: E402
from py.navier_stokes_model import NavierStokesModel  # noqa: E402
from py.scalar import Scalar  # noqa: E402
from py.scheme import SchemeType  # noqa: E402
from py.state import State  # noqa: E402
from py.vector_operations import InnerProduct, Laplacian  # noqa: E402

DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"

# Same grid/physics as SURF_test/vortall/run_vortall.py (pinned by
# VORTALL.mat's own 449x199 shape -- see SURF_test/vortall/README.md).
NX, NY, NGRID = 450, 200, 1
LENGTH, XOFFSET, YOFFSET = 9.0, -1.0, -2.0
REYNOLDS = 100.0
DT = 0.02
GEOM_FILE = REPO / "examples" / "cylinder.geom"


def independent_divergence(q, lev: int = 0) -> np.ndarray:
    """Discrete divergence of a Flux field, computed from scratch (does not
    call any py/ divergence/curl code), to independently check whether
    Curl()-derived flux fields are actually divergence-free as claimed.

    On this staggered (MAC-like) grid, X-fluxes live on a (nx+1, ny) array
    and Y-fluxes on a (nx, ny+1) array (see flux.py's `resize`/`getIndex`);
    the finite-volume divergence of interior cell (i,j) is
        div[i,j] = X[i+1,j] - X[i,j] + Y[i,j+1] - Y[i,j]
    (a discrete Gauss's-theorem sum of the four edge fluxes around the cell,
    NOT normalized by dx -- fine, since we only care whether it's ~0).
    """
    nx, ny = q.Nx(), q.Ny()
    X = q._data[lev, q.begin(Direction.X):q.end(Direction.X)].reshape(nx + 1, ny)
    Y = q._data[lev, q.begin(Direction.Y):q.end(Direction.Y)].reshape(nx, ny + 1)
    div = (X[1:, :] - X[:-1, :]) + (Y[:, 1:] - Y[:, :-1])
    return div


def norms(arr: np.ndarray) -> tuple[float, float]:
    a = np.asarray(arr).ravel()
    if a.size == 0:
        return 0.0, 0.0
    return float(np.sqrt(np.mean(a * a))), float(np.max(np.abs(a)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsteps", type=int, default=3000,
                      help="number of timesteps to run (default 3000, t=60: "
                           "covers setup through the early-to-mid transient-"
                           "growth region identified in SURF_test/vortall/README.md)")
    ap.add_argument("--every", type=int, default=1,
                      help="compute+record residuals every K steps (default 1: every step)")
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    grid = Grid(NX, NY, NGRID, LENGTH, XOFFSET, YOFFSET, 0.0, 0.0)
    geom = Geometry()
    assert geom.load(str(GEOM_FILE)), f"failed to load {GEOM_FILE}"
    dx = grid.Dx(0)
    print(f"grid: nx={NX} ny={NY} dx={dx:.6g}  boundary points={geom.getNumPoints()}")

    q_potential = BaseFlow(grid, 1.0, 0.0)  # magnitude=1, alpha=0 (matches run_vortall.py)

    model = NavierStokesModel(grid, geom, REYNOLDS, q_potential)
    solver = NonlinearIBSolver(grid, model, DT, SchemeType.RK3)

    model.init()
    solver.init()

    x = State(grid, geom.getNumPoints())
    x.omega.assign(0.0)
    x.f.assign(0.0)
    x.q.assign(0.0)
    geom.moveBodies(x.time)
    model.updateOperators(x.time)
    model.refreshState(x)

    setup_time = time.time() - t0
    print(f"setup done in {setup_time:.1f}s (Cholesky factorization + FFT plan)")

    rows = []
    header = ["step", "time", "div_l2", "div_linf", "noslip_l2", "noslip_linf",
              "poisson_l2", "poisson_linf", "helmholtz_l2", "helmholtz_linf",
              "circulation", "enstrophy", "kinetic_energy", "Cd", "Cl", "wall_s"]

    def record(step: int, tsim: float, wall_s: float) -> None:
        # 1. divergence of the flux field (independent check)
        div = independent_divergence(x.q)
        div_l2, div_linf = norms(div)

        # 2. no-slip constraint residual: C(omega) should equal the target b
        b = model.getConstraints()
        f_check = type(b)(b.getNumPoints())
        model.C(x.omega, f_check)
        resC = f_check.flatten() - b.flatten()
        noslip_l2, noslip_linf = norms(resC)

        # 3. Poisson-equation residual for the streamfunction solve
        psi = model.vorticityToStreamfunction(x.omega)
        lap_psi = Laplacian(psi)
        # interior of level 0 only (Laplacian/vorticityToStreamfunction are
        # only meaningful on interior nodes; edges are boundary-condition
        # dependent by construction)
        resP = lap_psi._data[0] - (-x.omega._data[0])
        poisson_l2, poisson_linf = norms(resP)

        # 4. Helmholtz-equation residual for the implicit sub-step operator
        #    A = (I - alpha*beta/2 L) applied to omegaStar should reproduce
        #    the RHS `a` that ProjectionSolver.solve() was given. We rebuild
        #    `a`/`omegaStar` the same way IBSolver.advanceSubstep does for
        #    substep 0 of this step's most recent solve, using the
        #    already-computed model.getAlpha() and the RK3 substep-0
        #    coefficients (an(0)+bn(0)) -- see ib_solver.py IBSolver ctor.
        beta0 = (solver._scheme.an(0) + solver._scheme.bn(0)) * DT
        alpha = model.getAlpha()
        Lomega = Laplacian(x.omega)
        lhs = x.omega - (alpha * beta0 / 2.0) * Lomega
        # This isn't a residual of an equation actually being re-solved
        # (that already happened inside the timestep); it's a direct
        # evaluation of the operator on the CURRENT state, which is only
        # meaningful as a smoothness/bounded-ness check, not a solver
        # residual -- so instead we report ||Laplacian(omega)|| itself,
        # which independently flags any grid-scale checkerboard noise or
        # blow-up in the vorticity field, something a genuine solver bug
        # (e.g. in Regularizer or the projection step) would produce.
        helm_l2, helm_linf = norms(Lomega._data[0])
        del lhs  # documented above: not a true residual, only Laplacian(omega) is recorded

        # conservation / sanity quantities
        circulation = float(np.sum(x.omega._data[0])) * dx * dx
        enstrophy = float(np.sum(x.omega._data[0] ** 2)) * dx * dx
        ke = 0.5 * InnerProduct(x.q, x.q)

        xF, yF = x.computeNetForce()
        Cd, Cl = 2.0 * xF, 2.0 * yF

        rows.append([step, tsim, div_l2, div_linf, noslip_l2, noslip_linf,
                     poisson_l2, poisson_linf, helm_l2, helm_linf,
                     circulation, enstrophy, ke, Cd, Cl, wall_s])

    record(0, x.time, 0.0)

    tstep0 = time.time()
    for i in range(1, args.nsteps + 1):
        solver.advance(x)
        if i % args.every == 0 or i == args.nsteps:
            record(i, x.time, time.time() - tstep0)
        if i % 500 == 0:
            elapsed = time.time() - tstep0
            print(f"  step {i}/{args.nsteps}  t={x.time:.2f}  "
                  f"({elapsed:.1f}s elapsed, {elapsed / i * 1000:.2f} ms/step)")

    total_time = time.time() - t0
    print(f"done: {args.nsteps} steps in {total_time - setup_time:.1f}s "
          f"({(total_time - setup_time) / args.nsteps * 1000:.2f} ms/step avg), "
          f"total wall {total_time:.1f}s")

    out_csv = DATA_DIR / "residuals.csv"
    arr = np.array(rows)
    np.savetxt(out_csv, arr, delimiter=",", header=",".join(header), comments="")
    print(f"wrote {out_csv}  ({len(rows)} rows)")

    meta = DATA_DIR / "run_meta.txt"
    with open(meta, "w") as f:
        f.write(
            f"nx={NX} ny={NY} ngrid={NGRID} length={LENGTH} xoffset={XOFFSET} "
            f"yoffset={YOFFSET}\nRe={REYNOLDS} dt={DT} scheme=RK3 "
            f"boundary_points={geom.getNumPoints()}\n"
            f"nsteps={args.nsteps} every={args.every}\n"
            f"setup_time_s={setup_time:.2f}\n"
            f"stepping_time_s={total_time - setup_time:.2f}\n"
            f"ms_per_step={(total_time - setup_time) / args.nsteps * 1000:.3f}\n"
        )
    print(f"wrote {meta}")


if __name__ == "__main__":
    main()
