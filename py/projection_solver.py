# projection_solver.py
#
# Python port of src/ProjectionSolver.h / src/ProjectionSolver.cc
#
# Solve systems arising in a fractional step method.
#
# Solve a system of the form
#     (1 - alpha*beta/2 L) x + beta B f = a
#     C x = b
# for f and x, using a fractional step method:
#     A x^*        = a
#     C A^{-1} B_1 f = C x^* - b
#     x            = x^* - A^{-1} B_1 f
# where
#     A   = 1 - alpha*beta/2 L
#     B_1 = beta B
#
# This is an abstract base class, and may not be instantiated directly.
#
# Original author: Clancy Rowley (3 Jul 2008)
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES (see also inline NOTE(port) comments):
#
#   1. Abstract base class.  C++ makes ProjectionSolver abstract through the
#      pure virtual `Minv(...) = 0`. Here the class derives from abc.ABC and
#      Minv is an @abstractmethod, so ProjectionSolver cannot be instantiated
#      directly -- only its subclasses (CholeskySolver,
#      ConjugateGradientSolver) can be.
#
#   2. Destructor dropped.  The C++ `~ProjectionSolver()` is empty; Python is
#      garbage collected, so there is no counterpart.
#
#   3. `const Grid _grid` is a *value* member in C++ (the constructor copies
#      the caller's Grid), whereas `const NavierStokesModel& _model` is a
#      reference. We store a plain reference to `grid` here rather than
#      deep-copying it: Grid is treated as immutable throughout this codebase
#      (it is only ever read, e.g. to size Scalars in M()), so a copy vs. a
#      reference is behaviorally identical. Flagged because it is a
#      deliberate departure from the literal C++ storage class.
#
#   4. Out-parameters.  Several methods take an output object that C++ mutates
#      through a non-const reference (`Scalar& y`, `BoundaryVector& y`, etc.).
#      Python cannot rebind a caller's name, so we mutate the passed-in object
#      in place: `.assign(...)` for a whole-object copy (C++ `operator=`), or
#      the object's own in-place operators. This matches the convention
#      already used in the other ported files (see navier_stokes_model.py).
#
#   5. Overload collapse on M().  C++ overloads M():
#         void          M(const BoundaryVector& f, BoundaryVector& y)  -- in-place
#         BoundaryVector M(const BoundaryVector& f)                    -- convenience
#      Collapsed into a single `M(f, y=None)`: when y is omitted the
#      convenience form runs (allocates y, computes into it, returns it); when
#      y is supplied the in-place form runs (writes into y, returns None to
#      match the C++ `void`).
#
#   6. save()/load() take the *basename* (a str) and return False in the base
#      class, exactly as C++ ("not implemented: return false"). Subclasses
#      override.
#
#   7. JAX-readiness.  This class holds no array math of its own; it only
#      orchestrates Scalar / BoundaryVector / HelmholtzSolver / model calls.
#      Nothing here needs changing for a later JAX port beyond whatever those
#      collaborators require.
# ---------------------------------------------------------------------------

from __future__ import annotations

import abc
from typing import Optional

from .boundary_vector import BoundaryVector
from .elliptic_solver import HelmholtzSolver
from .grid import Grid
from .navier_stokes_model import NavierStokesModel
from .scalar import Scalar


class ProjectionSolver(abc.ABC):
    """Abstract base class: solve systems arising in a fractional step
    method (see module docstring)."""

    def __init__(self, grid: Grid, model: NavierStokesModel, beta: float) -> None:
        # Constructor: initialize Helmholtz solver to solve
        # (1 - beta/2 L) u = f
        self._beta: float = beta
        # NOTE(port): reference stored rather than a copy -- see port note 3.
        self._grid: Grid = grid
        self._model: NavierStokesModel = model
        self._helmholtz: HelmholtzSolver = HelmholtzSolver(
            grid, -beta / 2.0 * model.getAlpha()
        )

    # Initialization for this base class done in constructor.
    # Subclasses use this for their own initialization, if needed.
    def init(self) -> None:
        pass

    # Subclasses can override these to save and load their own state
    # (e.g. a Cholesky factorization). By default, not implemented: return
    # False.
    def save(self, basename: str) -> bool:
        return False

    def load(self, basename: str) -> bool:
        return False

    def solve(
        self,
        a: Scalar,
        b: BoundaryVector,
        omega: Scalar,
        f: BoundaryVector,
    ) -> None:
        """Solve for omega and f for a system of the form
            A omega + B f = a
            C omega       = b
        using a fractional step method:
            A omega^*        = a
            C A^{-1} B f     = C omega^* - b
            omega            = omega^* - A^{-1} B f

        NOTE(port): `omega` and `f` are out-parameters, mutated in place
        (port note 4).
        """
        # A omega^* = a
        omegaStar = Scalar(a.getGrid())
        self.Ainv(a, omegaStar)

        # C A^{-1} B f = C omega^* - b
        rhs = BoundaryVector(f.getNumPoints())
        self.C(omegaStar, rhs)
        rhs -= b                    # rhs = C omega^* - b
        self.Minv(rhs, f)           # f = Minv( rhs )

        # omega = omega^* - A^{-1} B f
        c = Scalar(a.getGrid())
        self.B(f, c)                # c = Bf
        self.Ainv(c, c)             # c = Ainv(Bf)
        omega.assign(omegaStar - c)

    # -- protected methods ----------------------------------------------

    def Ainv(self, b: Scalar, y: Scalar) -> None:
        """Solve y = A^{-1} b."""
        self._helmholtz.solve(b, y)

    def B(self, f: BoundaryVector, y: Scalar) -> None:
        """Compute y = B(f)."""
        self._model.B(f, y)
        y *= self._beta

    def C(self, x: Scalar, y: BoundaryVector) -> None:
        """Compute y = C(x)."""
        self._model.C(x, y)

    def M(
        self, f: BoundaryVector, y: Optional[BoundaryVector] = None
    ) -> Optional[BoundaryVector]:
        """Compute y = M(f), where M = C A^{-1} B.

        NOTE(port): overload collapse (port note 5).
            M(f)     -> convenience form: allocate y, compute, return it
            M(f, y)  -> in-place form: compute into y, return None
        """
        if y is None:
            y = BoundaryVector(f.getNumPoints())
            self.M(f, y)
            return y

        u = Scalar(self._grid)
        self.B(f, u)                # u = B f
        self.Ainv(u, u)             # u = Ainv B f
        self.C(u, y)                # y = C Ainv B f
        return None

    @abc.abstractmethod
    def Minv(self, b: BoundaryVector, x: BoundaryVector) -> None:
        """Compute x = M^{-1} b."""
        ...
