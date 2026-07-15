# bc.py
#
# Python port of src/BC.h / src/BC.cc
#
# Class for storing and accessing boundary conditions for a 2d scalar field.
#
# Boundary data is stored in a 1d array. For an (8 x 4) grid, the data is
# arranged as follows:
#
#  4  5  6  7  8  9 10 11 12
#  3                      13
#  2                      14
#  1                      15
#  0 23 22 21 20 19 18 17 16
#
# Note that the total number of boundary points is
#   (nx+1) * 2 + (ny-1) * 2 = 2 * (nx + ny)

from __future__ import annotations

import numpy as np


class BC:
    """Boundary conditions for a 2d scalar field, stored in a 1d array."""

    def __init__(self, nx: int, ny: int) -> None:
        self._nx = nx
        self._ny = ny
        # NOTE(port): C++ Array::Array1<double> is 0-based here (no offset
        # passed to the constructor), so this maps directly onto a
        # zero-indexed numpy array of the same size.
        self._data = np.zeros(2 * (nx + ny), dtype=np.float64)

    def left(self, j: int) -> float:
        """Return the value on the left boundary, at index j (0..ny)."""
        assert j >= 0 and j <= self._ny
        return self._data[j]

    def setLeft(self, j: int, value: float) -> None:
        # NOTE(port): C++ `left(j)` returns `double&`, usable as an lvalue
        # (e.g. `bc.left(j) = value`). Python has no reference return, so
        # writes are done via this separate `set*` method instead of via
        # the accessor. See also right/top/bottom below.
        assert j >= 0 and j <= self._ny
        self._data[j] = value

    def right(self, j: int) -> float:
        """Return the value on the right boundary, at index j (0..ny)."""
        assert j >= 0 and j <= self._ny
        return self._data[2 * self._ny + self._nx - j]

    def setRight(self, j: int, value: float) -> None:
        assert j >= 0 and j <= self._ny
        self._data[2 * self._ny + self._nx - j] = value

    def top(self, i: int) -> float:
        """Return the value on the top boundary, at index i (0..nx)."""
        assert i >= 0 and i <= self._nx
        return self._data[self._ny + i]

    def setTop(self, i: int, value: float) -> None:
        assert i >= 0 and i <= self._nx
        self._data[self._ny + i] = value

    def bottom(self, i: int) -> float:
        """Return the value on the bottom boundary, at index i (0..nx)."""
        assert i >= 0 and i <= self._nx
        if i == 0:
            return self._data[0]
        else:
            return self._data[2 * (self._nx + self._ny) - i]

    def setBottom(self, i: int, value: float) -> None:
        assert i >= 0 and i <= self._nx
        if i == 0:
            self._data[0] = value
        else:
            self._data[2 * (self._nx + self._ny) - i] = value

    def assign(self, a: float) -> "BC":
        """Set every boundary value to `a` (C++ operator=(double))."""
        # NOTE(port): C++ overloads `operator=(double a)` so `bc = a;`
        # zero-fills (or fills with any scalar) the whole array. Python
        # cannot overload plain assignment, so this is exposed as an
        # explicit `assign()` method instead.
        self._data[:] = a
        return self

    def Nx(self) -> int:
        return self._nx

    def Ny(self) -> int:
        return self._ny
