# conjugate_gradient_solver.py
#
# Python port of src/ConjugateGradientSolver.h / src/ConjugateGradientSolver.cc
#
# Subclass of ProjectionSolver in which the system M f = b is solved
# iteratively, using a conjugate-gradient method. Assumes M is symmetric, and
# iterates until a specified tolerance has been reached.
#
# Original author: Clancy Rowley (8 Aug 2008)
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES (see also inline NOTE(port) comments):
#
#   1. `numIterations` is never incremented -- FAITHFUL REPRODUCTION OF THE
#      ORIGINAL.  In the C++ Minv(), `numIterations` is initialized to 0 and
#      the loop condition tests `numIterations < MAX_ITERATIONS`, but the loop
#      body never increments it. The iteration cap is therefore dead code in
#      the original: the loop runs purely until `delta <= _toleranceSquared`.
#      This looks like a bug, but the task is a faithful port, so it is
#      reproduced exactly (numIterations stays 0, MAX_ITERATIONS never trips).
#      Flagged so it is a conscious decision rather than a silent one.
#
#   2. Out-parameter `x`.  As elsewhere in this port, `x` is mutated in place;
#      `x += alpha * d` uses BoundaryVector.__iadd__, which mutates the
#      caller's object (matching the C++ non-const reference).
#
#   3. Overloaded M().  The C++ calls the convenience overload `M(x)` /
#      `M(d)` that returns a new BoundaryVector; the Python ProjectionSolver.M
#      supplies the same convenience form when called with a single argument.
#
#   4. setTolerance/getTolerance were `inline` in the header; here they are
#      ordinary methods with identical behavior.
#
#   5. JAX-readiness.  The per-iteration work is BoundaryVector arithmetic and
#      InnerProduct (numpy dot); the Python-level loop is the CG iteration
#      itself, whose length is data-dependent. A later JAX port would express
#      this with jax.lax.while_loop; the body already maps cleanly onto array
#      ops.
# ---------------------------------------------------------------------------

from __future__ import annotations

import math

from .boundary_vector import BoundaryVector, InnerProduct
from .grid import Grid
from .navier_stokes_model import NavierStokesModel
from .projection_solver import ProjectionSolver

MAX_ITERATIONS = 3000


class ConjugateGradientSolver(ProjectionSolver):
    """Solve M f = b iteratively via a conjugate-gradient method."""

    def __init__(
        self,
        grid: Grid,
        model: NavierStokesModel,
        beta: float,
        tolerance: float,
    ) -> None:
        # Constructor. Store the tolerance (squared) as private data.
        super().__init__(grid, model, beta)
        self._toleranceSquared: float = tolerance * tolerance

    def setTolerance(self, tolerance: float) -> None:
        self._toleranceSquared = tolerance * tolerance

    def getTolerance(self) -> float:
        return math.sqrt(self._toleranceSquared)

    def Minv(self, b: BoundaryVector, x: BoundaryVector) -> None:
        """Solve M f = b for f iteratively, using a conjugate-gradient method.
        Assumes M is symmetric.

        NOTE(port): `x` is an out-parameter, updated in place.
        """
        # NOTE(port): numIterations is never incremented below -- see port
        # note 1.
        numIterations = 0

        # r = b - M(x)
        r = BoundaryVector(b)
        q = self.M(x)
        r -= q
        d = BoundaryVector(r)
        delta = InnerProduct(r, r)

        # while error is greater than tolerance
        while (delta > self._toleranceSquared) and (numIterations < MAX_ITERATIONS):
            # alpha = r^2 / <d, Md>
            q = self.M(d)
            alpha = delta / InnerProduct(d, q)
            x += alpha * d
            r -= alpha * q
            delta_old = delta
            delta = InnerProduct(r, r)
            beta = delta / delta_old
            d = r + beta * d
