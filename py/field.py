# field.py
#
# Python port of src/Field.h / src/Field.cc
#
# Abstract base class for a field defined on a grid: scalar fields (Scalar)
# and vector fields described by Flux objects are subclasses of a Field.
#
# This class holds only a Grid (scalar parameters), so there is no array
# data and no numpy vectorization applies here.

from __future__ import annotations

from typing import Optional

from .grid import Grid


class Field:
    """Base class for a field defined on a grid."""

    def __init__(self, grid: Optional[Grid] = None) -> None:
        # NOTE(port): C++ has three constructors: Field() (default, builds
        # a Grid with nx=ny=-1 to avoid divide-by-zero while marking "no
        # grid defined"), Field(const Grid&), and the copy constructor
        # Field(const Field&). Python collapses these into one constructor:
        # `grid=None` reproduces Field(), and `grid=<Grid or Field>`
        # reproduces both Field(const Grid&) (copies the Grid) and
        # Field(const Field&) (copies the other Field's Grid), since a
        # Grid copy is a plain attribute copy either way.
        if grid is None:
            # Note: cannot set nx to zero or computation of dx will divide
            # by zero. Set to -1 to indicate no grid defined.
            nx = -1
            ny = -1
            ngrid = 1
            length = 0.0
            xOffset = 0.0
            yOffset = 0.0
            self._grid = Grid()
            self._grid.resize(nx, ny, ngrid, length, xOffset, yOffset)
        elif isinstance(grid, Field):
            self._grid = grid._grid
        else:
            self._grid = grid

    def Nx(self) -> int:
        return self._grid.Nx()

    def NxExt(self) -> int:
        return self._grid.NxExt()

    def Ny(self) -> int:
        return self._grid.Ny()

    def NyExt(self) -> int:
        return self._grid.NyExt()

    def Ngrid(self) -> int:
        return self._grid.Ngrid()

    def Dx(self, lev: Optional[int] = None) -> float:
        # NOTE(port): collapses the two C++ overloads Dx() and Dx(int lev)
        # into one method, same as Grid.Dx() (see grid.py).
        return self._grid.Dx(lev)

    def getXShift(self) -> float:
        return self._grid.getXShift()

    def getYShift(self) -> float:
        return self._grid.getYShift()

    def getGrid(self) -> Grid:
        return self._grid

    def setGrid(self, grid: Grid) -> None:
        self._grid = grid

    def getXCenter(self, lev: int, i: int) -> float:
        """Return the x-coordinate of the center of cell i  (i in 0..m-1)."""
        return self._grid.getXCenter(lev, i)

    def getYCenter(self, lev: int, j: int) -> float:
        """Return the y-coordinate of the center of cell j  (j in 0..n-1)."""
        return self._grid.getYCenter(lev, j)

    def getXEdge(self, lev: int, i: int) -> float:
        """Return the x-coordinate of the left edge of cell i  (i in 0..m)."""
        return self._grid.getXEdge(lev, i)

    def getYEdge(self, lev: int, j: int) -> float:
        """Return the y-coordinate of the bottom edge of cell j  (j in 0..n)."""
        return self._grid.getYEdge(lev, j)
