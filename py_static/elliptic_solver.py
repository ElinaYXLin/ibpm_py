# elliptic_solver.py
#
# Python port of src/EllipticSolver.h / src/EllipticSolver.cc
#
# Solve Poisson and Helmholtz equations on a multi-domain grid (Scalar).
#
# EllipticSolver is an abstract base class that drives a nested family of
# single-domain solvers (EllipticSolver2d, one per grid level). PoissonSolver
# and HelmholtzSolver are the concrete subclasses.
#
# Original author: Clancy Rowley (30 Sep 2008)
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES (see also inline NOTE(port) comments):
#
#   1. Abstract base class.  C++ makes EllipticSolver abstract via a pure
#      virtual destructor and a pure virtual create2dSolver(). Here it derives
#      from abc.ABC and create2dSolver is an @abstractmethod, so EllipticSolver
#      cannot be instantiated directly.
#
#   2. Destructor dropped.  The C++ destructor only `delete`s the per-level
#      EllipticSolver2d objects held in `_solvers`. Python is garbage
#      collected, so there is no counterpart and it is omitted (see also
#      port note 3 in elliptic_solver_2d.py, where the 2d solver's own
#      destructor was likewise dropped).
#
#   3. Overload collapse on solve().  C++ overloads solve():
#         void   solve(const Scalar& f, Scalar& u)   -- in-place form
#         Scalar solve(const Scalar& f)              -- convenience form
#      Python has no overloading, so both are collapsed into a single
#      solve(f, u=None): when u is omitted the convenience form runs
#      (allocates a new Scalar, solves into it, and returns it); when u is
#      supplied the in-place form runs (writes into u, returns None to match
#      the C++ `void`). This mirrors the collapse pattern used in the other
#      ported files.
#
#   4. `vector<EllipticSolver2d*> _solvers`.  Becomes a Python list. The C++
#      constructor sizes it to Ngrid null pointers; here it is initialized to
#      a list of Ngrid `None` entries, which init() then fills. The stored
#      elements are concrete EllipticSolver2d instances (no raw pointers).
#
#   5. Array2d slices.  `rhs[lev]` / `u[lev]` return numpy views of shape
#      (Nx-1, Ny-1) from Scalar.__getitem__ (see scalar.py). The 2d solver
#      writes its result through those views in place, exactly as the C++
#      `Array::Array2<double>` slices do, so the multi-grid result lands back
#      in the Scalar `u` with no extra copy.
#
#   6. JAX-readiness.  This class contains no array math of its own; it only
#      orchestrates Scalar/EllipticSolver2d calls and does integer/bit
#      arithmetic. Nothing here needs changing for a later JAX port beyond
#      whatever the Scalar / EllipticSolver2d ports require.
# ---------------------------------------------------------------------------

from __future__ import annotations

import abc
from typing import List, Optional

import numpy as np

from .bc import BC
from .elliptic_solver_2d import (
    EllipticSolver2d,
    HelmholtzSolver2d,
    PoissonSolver2d,
)
from .grid import Grid
from .scalar import Scalar


class EllipticSolver(abc.ABC):
    """Abstract base class: solve Poisson/Helmholtz equations on a
    multi-domain grid, using one EllipticSolver2d per grid level."""

    def __init__(self, grid: Grid) -> None:
        self._ngrid: int = grid.Ngrid()
        self._dx: float = grid.Dx()
        # NOTE(port): C++ `_solvers( grid.Ngrid() )` sizes the vector to
        # Ngrid null pointers; init() fills them in. See port note 4.
        self._solvers: List[Optional[EllipticSolver2d]] = [None] * grid.Ngrid()

    # Create 2d solvers
    def init(self) -> None:
        for lev in range(self._ngrid):
            # calculate grid spacing on this grid level
            dx = self._dx * (1 << lev)
            self._solvers[lev] = self.create2dSolver(dx)

    def solve(self, f: Scalar, u: Optional[Scalar] = None) -> Optional[Scalar]:
        # NOTE(port): overload collapse (port note 3).
        #   solve(f)     -> convenience form: allocate u, solve, return it
        #   solve(f, u)  -> in-place form: solve into u, return None
        if u is None:
            u = Scalar(f.getGrid())
            self.solve(f, u)
            return u

        # Multi-domain elliptic solver
        assert f.Ngrid() == self._ngrid
        assert f.Ngrid() == u.Ngrid()
        assert f.Nx() == u.Nx()
        assert f.Ny() == u.Ny()

        # First "coarsify" right-hand side (f) to coarse grids
        rhs = Scalar(f)
        rhs.coarsify()
        # Solve coarsest grid first, then finer grids
        for lev in range(f.Ngrid() - 1, -1, -1):
            # Get slices of input and output data at current grid level
            rhs1: np.ndarray = rhs[lev]
            u1: np.ndarray = u[lev]
            # if on the coarsest grid, solve with zero bcs
            if lev == f.Ngrid() - 1:
                self._solvers[lev].solve(rhs1, u1)
            else:
                # Get boundary condition from next coarser grid
                bc = BC(f.Nx(), f.Ny())
                u.getBC(lev, bc)
                # solve with specified boundary conditions
                self._solvers[lev].solve(rhs1, bc, u1)
        return None

    @abc.abstractmethod
    def create2dSolver(self, dx: float) -> EllipticSolver2d:
        ...


# -----------------------------------------------------------------------------

class PoissonSolver(EllipticSolver):
    """Solve a Poisson equation L u = f on a multi-domain grid, with zero
    boundary conditions on the outer domain of u, where L is the Laplacian."""

    def __init__(self, grid: Grid) -> None:
        super().__init__(grid)
        self._nx: int = grid.Nx()
        self._ny: int = grid.Ny()
        self.init()

    def create2dSolver(self, dx: float) -> EllipticSolver2d:
        return PoissonSolver2d(self._nx, self._ny, dx)


# -----------------------------------------------------------------------------

class HelmholtzSolver(EllipticSolver):
    """Solve a Helmholtz equation (1 + alpha * L) u = f on a multi-domain
    grid, with zero boundary conditions on u, where L is the Laplacian."""

    def __init__(self, grid: Grid, alpha: float) -> None:
        super().__init__(grid)
        self._nx: int = grid.Nx()
        self._ny: int = grid.Ny()
        self._alpha: float = alpha
        self.init()

    def create2dSolver(self, dx: float) -> EllipticSolver2d:
        return HelmholtzSolver2d(self._nx, self._ny, dx, self._alpha)
