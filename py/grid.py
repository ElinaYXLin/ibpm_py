# grid.py
#
# Python port of src/Grid.h / src/Grid.cc
#
# Faithful, line-by-line port of the ibpm::Grid class. This class holds only
# scalar parameters (no array data), so there is no numpy vectorization to
# apply here; the "vectorize, no manual loops" requirement is a no-op for
# this particular file. The class is kept as a plain Python object (no numpy
# arrays, no JAX pytree registration) so that it can later be wrapped as a
# static/hashable config object for JAX (e.g. as a NamedTuple or a
# dataclass registered as a JAX pytree with static fields), without needing
# to change its public interface now.

from __future__ import annotations

import math
from typing import Tuple


class Grid:
    """Define parameters associated with a uniform staggered grid.

    Uses a staggered grid, suitable for a finite-volume method, in which
    scalars (such as streamfunction and vorticity) are defined at cell nodes
    and fluxes are defined at cell edges.
    """

    def __init__(
        self,
        nx: int = 0,
        ny: int = 0,
        ngrid: int = 0,
        length: float = 0.0,
        xOffset: float = 0.0,
        yOffset: float = 0.0,
        xShift: float = 0.0,
        yShift: float = 0.0,
    ) -> None:
        # NOTE(port): The C++ class has three overloaded constructors:
        #   Grid(nx, ny, ngrid, length, xOffset, yOffset)
        #   Grid(nx, ny, ngrid, length, xOffset, yOffset, xShift, yShift)
        #   Grid()  (default, all zeros)
        # Python has no overloading, so these are collapsed into a single
        # constructor with defaults reproducing Grid() when called with no
        # arguments, and reproducing the 6-arg constructor when xShift/
        # yShift are left at their default of 0.0 (matching the 6-arg
        # C++ constructor, which also sets _xShift = _yShift = 0).
        if nx == 0 and ny == 0 and ngrid == 0 and length == 0.0 and xOffset == 0.0 and yOffset == 0.0 and xShift == 0.0 and yShift == 0.0:
            # Default constructor: set all parameters to zero
            self._nx: int = 0
            self._ny: int = 0
            self._ngrid: int = 0
            self._xOffset: float = 0.0
            self._yOffset: float = 0.0
            self._dx: float = 0.0
            self._xShift: float = 0.0
            self._yShift: float = 0.0
            return

        self._xShift = 0.0
        self._yShift = 0.0
        self.resize(nx, ny, ngrid, length, xOffset, yOffset)
        if xShift != 0.0 or yShift != 0.0:
            self.setXShift(xShift)
            self.setYShift(yShift)

    def resize(
        self,
        nx: int,
        ny: int,
        ngrid: int,
        length: float,
        xOffset: float,
        yOffset: float,
        xShift: float = None,
        yShift: float = None,
    ) -> None:
        """Set all grid parameters, optionally including x- and y- shifts."""
        assert ngrid == 1 or (nx % 4 == 0 and ny % 4 == 0)
        self._nx = nx
        self._ny = ny
        self._ngrid = ngrid
        self._xOffset = xOffset
        self._yOffset = yOffset
        if xShift is not None:
            self._xShift = xShift
        if yShift is not None:
            self._yShift = yShift
        self._dx = length / nx

    def Nx(self) -> int:
        """Return number of cells in x-direction."""
        return self._nx

    def NxExt(self) -> int:
        """Return number of coarse cells outside each fine domain, in x-direction.

        This counts the number of cells to the left of the fine domain.
        This distinction is key when the grid is shifted.
        To round to the nearest integer we add 0.5 and use the floor function.
        However, due to asserts enforced in Grid.cc, this is nominally an
        integer already.
        """
        # NOTE(port): C++ `_nx / 4` is integer division (both operands are
        # int), truncating toward zero, and then the *result* is used in
        # floating point arithmetic with `1 - _xShift`. Python 3's `/`
        # always produces a float, so we must use `//` (floor division) to
        # reproduce C++ integer division. Since _nx is asserted to be
        # divisible by 4 elsewhere, `//` and `/` truncation agree here for
        # nonnegative _nx.
        return int(math.floor(self._nx // 4 * (1 - self._xShift) + 0.5))

    def Ny(self) -> int:
        """Return number of cells in y-direction."""
        return self._ny

    def NyExt(self) -> int:
        """Return number of coarse cells outside each fine domain, in y-direction."""
        return int(math.floor(self._ny // 4 * (1 - self._yShift) + 0.5))

    def c2f(self, i: int, j: int) -> Tuple[int, int]:
        """Given indices (i,j) on a coarse grid, return corresponding indices
        (ii,jj) on the fine grid."""
        # NOTE(port): C++ signature returns via reference out-params
        # (int& ii, int& jj); Python has no out-params, so this returns a
        # tuple (ii, jj) instead. Callers must be updated accordingly.
        ii = (i - self.NxExt()) * 2
        jj = (j - self.NyExt()) * 2
        return ii, jj

    def f2c(self, ii: int, jj: int) -> Tuple[int, int]:
        """Given indices (ii,jj) on a fine grid, return indices (i,j) of the
        corresponding point on the coarse grid, or of the nearest point
        below and to the left."""
        # NOTE(port): C++ `ii/2` and `jj/2` are integer division (ii, jj are
        # int); reproduced here with `//`. Also returns a tuple instead of
        # using reference out-params (see c2f above).
        i = ii // 2 + self.NxExt()
        j = jj // 2 + self.NyExt()
        return i, j

    def Ngrid(self) -> int:
        """Return number of grid levels for multi-domain solution."""
        return self._ngrid

    def Dx(self, lev: int = None) -> float:
        """Return grid spacing on finest level (same in x- and y-directions),
        or at the specified grid level if `lev` is given."""
        # NOTE(port): C++ overloads `Dx()` (finest level) and `Dx(int lev)`
        # (specified level) as two separate methods. Collapsed here into one
        # method with `lev=None` meaning "finest level", matching `Dx()`.
        if lev is None:
            return self._dx
        assert lev >= 0 and lev < self._ngrid
        return self._dx * (1 << lev)

    def getXCenter(self, lev: int, i: int) -> float:
        """Return the x-coordinate of the center of cell i  (i in 0..m-1)."""
        assert lev >= 0 and lev < self._ngrid
        assert i >= 0 and i <= self._nx
        return self._getXOffset(lev) + (i + 0.5) * self.Dx(lev)

    def getYCenter(self, lev: int, j: int) -> float:
        """Return the y-coordinate of the center of cell j  (j in 0..n-1)."""
        assert lev >= 0 and lev < self._ngrid
        assert j >= 0 and j <= self._ny
        return self._getYOffset(lev) + (j + 0.5) * self.Dx(lev)

    def getXEdge(self, lev: int, i: int) -> float:
        """Return the x-coordinate of the left edge of cell i  (i in 0..m)."""
        assert lev >= 0 and lev < self._ngrid
        assert i >= 0 and i <= self._nx
        return self._getXOffset(lev) + i * self.Dx(lev)

    def getYEdge(self, lev: int, j: int) -> float:
        """Return the y-coordinate of the bottom edge of cell j  (j in 0..n)."""
        assert lev >= 0 and lev < self._ngrid
        assert j >= 0 and j <= self._ny
        return self._getYOffset(lev) + j * self.Dx(lev)

    def getXGridIndex(self, x: float) -> int:
        """Return the grid index i corresponding to the given x-coordinate.
        Currently, only works for the finest grid level."""
        xpos = x - self._xOffset
        assert xpos <= self._dx * self._nx
        i = int(math.floor(xpos / self._dx))
        if (xpos - i * self._dx) >= (self._dx / 2):
            i = int(math.ceil(xpos / self._dx))
        return i

    def getYGridIndex(self, y: float) -> int:
        """Return the grid index j corresponding to the given y-coordinate.
        Currently, only works for the finest grid level."""
        ypos = y - self._yOffset
        assert ypos <= self._dx * self._ny
        j = int(math.floor(ypos / self._dx))
        if (ypos - j * self._dx) >= (self._dx / 2):
            j = int(math.ceil(ypos / self._dx))
        return j

    def setXShift(self, xShift: float) -> None:
        """Set shift parameter in x."""
        assert abs(xShift) <= 1
        assert math.fmod(xShift * self._nx, 4) == 0
        self._xShift = xShift

    def getXShift(self) -> float:
        """Get the current x-shift parameter."""
        return self._xShift

    def setYShift(self, yShift: float) -> None:
        """Set shift parameter in y."""
        assert abs(yShift) <= 1
        assert math.fmod(yShift * self._ny, 4) == 0
        self._yShift = yShift

    def getYShift(self) -> float:
        """Get the current y-shift parameter."""
        return self._yShift

    def isEqualTo(self, grid2: "Grid") -> bool:
        """Compare two grids."""
        # NOTE(port): The C++ implementation multiplies booleans together
        # (bool * bool * ... in C++ promotes bools to int and multiplies,
        # returning nonzero/zero which converts back to bool) instead of
        # short-circuit `&&`. This is reproduced faithfully with `and`,
        # which is semantically equivalent here since none of the operands
        # have side effects (all are pure comparisons).
        nx_eq = self._nx == grid2.Nx()
        ny_eq = self._ny == grid2.Ny()
        ngrid_eq = self._ngrid == grid2.Ngrid()
        dx_eq = self._dx == grid2.Dx()
        xOffset_eq = self._xOffset == grid2.getXEdge(0, 0)
        yOffset_eq = self._yOffset == grid2.getYEdge(0, 0)
        xShift_eq = self._xShift == grid2.getXShift()
        yShift_eq = self._yShift == grid2.getYShift()
        return bool(
            nx_eq and ny_eq and ngrid_eq and dx_eq
            and xOffset_eq and yOffset_eq and xShift_eq and yShift_eq
        )

    def _getXOffset(self, lev: int) -> float:
        """Return the x-coordinate of the left-most gridpoint of level lev."""
        return self._xOffset + 0.5 * (self._xShift - 1) * ((1 << lev) - 1) * (self._nx * self._dx)

    def _getYOffset(self, lev: int) -> float:
        """Return the y-coordinate of the bottom gridpoint of level lev."""
        return self._yOffset + 0.5 * (self._yShift - 1) * ((1 << lev) - 1) * (self._ny * self._dx)
