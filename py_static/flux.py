# flux.py
#
# Python port of src/Flux.h / src/Flux.cc
#
# Store a 2D array of fluxes, located at cell edges.
#
# For a grid with nx cells in the x-direction and ny cells in the
# y-direction, there are (nx+1,ny) fluxes in the x-direction, and (nx,ny+1)
# fluxes in the y-direction. These are accessible via q.x(...) and q.y(...).

from __future__ import annotations

import math
from typing import Optional, Union

import numpy as np

from .direction import Direction
from .field import Field
from .grid import Grid
from .tangent_se2 import TangentSE2

FloatArray = np.ndarray


class Flux(Field):
    """A 2D array of fluxes, located at cell edges."""

    # Type used for referencing elements (C++ `typedef int index`)
    index = int

    def __init__(self, arg: Optional[Union[Grid, "Flux"]] = None) -> None:
        # NOTE(port): as in scalar.py, the three C++ constructors
        # (Flux(const Grid&), Flux() default, Flux(const Flux&) copy) are
        # collapsed into a single constructor dispatching on the runtime
        # type of `arg`.
        if arg is None:
            super().__init__()
            self._numXFluxes: int = 0
            self._numFluxes: int = 0
            self._data: Optional[FloatArray] = None
            return
        if isinstance(arg, Flux):
            super().__init__(arg.getGrid())
            self.resize(arg.getGrid())
            self._data[...] = arg._data
            return
        super().__init__(arg)
        self.resize(arg)

    def resize(self, grid: Grid) -> None:
        """Set all parameters and reallocate arrays based on the Grid
        dimensions."""
        self.setGrid(grid)
        nx = self.Nx()
        ny = self.Ny()
        self._numXFluxes = nx * ny + ny
        self._numFluxes = 2 * nx * ny + nx + ny
        # NOTE(port): see the corresponding note in scalar.py `resize` --
        # C++ `_data.Allocate(...)` leaves memory uninitialized; we zero-
        # fill with np.zeros for deterministic behavior.
        self._data = np.zeros((self.Ngrid(), self._numFluxes), dtype=np.float64)

    def print(self) -> None:
        """Print the X and Y components to standard out (for debugging)."""
        print("X:")
        for lev in range(self.Ngrid()):
            for j in range(self.Ny() - 1, -1, -1):
                row = " ".join(str(self(lev, Direction.X, i, j)) for i in range(self.Nx() + 1))
                print(row)
            print()
        print("Y:")
        for lev in range(self.Ngrid()):
            for j in range(self.Ny(), -1, -1):
                row = " ".join(str(self(lev, Direction.Y, i, j)) for i in range(self.Nx()))
                print(row)
            print()

    def assign(self, other: Union["Flux", float]) -> "Flux":
        """Copy assignment, from another Flux or from a scalar value.

        NOTE(port): see Scalar.assign() -- Python cannot overload plain
        `=`, so both C++ `operator=(const Flux&)` and
        `operator=(double)` are exposed through this one method.
        """
        if isinstance(other, Flux):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data[...] = other._data
        else:
            self._data[...] = other
        return self

    def __call__(self, lev: int, dir_or_ind, i: Optional[int] = None, j: Optional[int] = None) -> float:
        """q(lev,dir,i,j) refers to the flux in direction dir (X or Y) at
        edge (i,j). q(lev,ind) refers to the value at the given flat index.

        NOTE(port): C++ overloads `operator()` for both (lev,dir,i,j) and
        (lev, index) signatures, each with `double&` and `const double`
        variants (read/write vs read-only). Python collapses these into one
        `__call__` that dispatches on argument count; as in scalar.py, this
        only supports reads (`double&` write-through has no Python
        equivalent) -- use `set(...)` for writes.
        """
        assert lev >= 0 and lev < self.Ngrid()
        if i is None and j is None:
            ind = dir_or_ind
            assert ind >= 0 and ind < self._numFluxes
            return float(self._data[lev, ind])
        dir_ = dir_or_ind
        return float(self._data[lev, self.getIndex(dir_, i, j)])

    def set(self, lev: int, dir_or_ind, i: Optional[int] = None, j: Optional[int] = None, value: float = None) -> None:
        """Write counterpart to __call__ -- see the NOTE(port) there."""
        assert lev >= 0 and lev < self.Ngrid()
        if i is None and j is None:
            ind = dir_or_ind
            assert ind >= 0 and ind < self._numFluxes
            self._data[lev, ind] = value
            return
        dir_ = dir_or_ind
        self._data[lev, self.getIndex(dir_, i, j)] = value

    def x(self, lev: int, dir_or_ind, i: Optional[int] = None) -> float:
        """q.x(lev,dir,i) refers to the x-coordinate of the flux (dir,i,j),
        or q.x(lev,ind) returns the x-coordinate of the flux specified by
        the flat index ind."""
        assert lev >= 0 and lev < self.Ngrid()
        if i is None:
            ind = dir_or_ind
            dir_ = 0 if ind < self._numXFluxes else 1
            i_ = (ind - dir_ * self._numXFluxes) // (self.Ny() + dir_)
            return self.x(lev, dir_, i_)
        dir_ = dir_or_ind
        assert dir_ >= Direction.X and dir_ <= Direction.Y
        assert i >= 0
        if dir_ == Direction.X:
            assert i < self.Nx() + 1
            return self.getXEdge(lev, i)
        else:
            assert i < self.Nx()
            return self.getXCenter(lev, i)

    def y(self, lev: int, dir_or_ind, j: Optional[int] = None) -> float:
        """q.y(lev,dir,j) refers to the y-coordinate of the flux (dir,i,j),
        or q.y(lev,ind) returns the y-coordinate of the flux specified by
        the flat index ind."""
        assert lev >= 0 and lev < self.Ngrid()
        if j is None:
            ind = dir_or_ind
            dir_ = 0 if ind < self._numXFluxes else 1
            i_ = (ind - dir_ * self._numXFluxes) // (self.Ny() + dir_)
            j_ = ind - dir_ * self._numXFluxes - i_ * (self.Ny() + dir_)
            return self.y(lev, dir_, j_)
        dir_ = dir_or_ind
        assert dir_ >= Direction.X and dir_ <= Direction.Y
        assert j >= 0
        if dir_ == Direction.X:
            assert j < self.Ny()
            return self.getYCenter(lev, j)
        else:
            assert j < self.Ny() + 1
            return self.getYEdge(lev, j)

    def begin(self, dir: Optional[int] = None) -> int:
        """Returns an index that refers to the first element (overall, or
        the first element in direction dir if given)."""
        if dir is None:
            return 0
        assert dir >= Direction.X and dir <= Direction.Y
        return dir * self._numXFluxes

    def end(self, dir: Optional[int] = None) -> int:
        """Returns an index one past the last element (overall, or one past
        the last element in direction dir if given)."""
        if dir is None:
            return self._numFluxes
        assert dir >= Direction.X and dir <= Direction.Y
        if dir == Direction.X:
            return self._numXFluxes
        else:
            return self._numFluxes

    def getIndex(self, dir: int, i: int, j: int) -> int:
        """Returns an index for the value in direction dir at point (i,j)."""
        assert dir >= Direction.X and dir <= Direction.Y
        assert i >= 0 and j >= 0
        assert (i < self.Nx() + 1) if (dir == Direction.X) else (i < self.Nx())
        assert (j < self.Ny() + 1) if (dir == Direction.Y) else (j < self.Ny())
        # Tricky expression:
        #   j in [0..ny-1] for X fluxes (dir = X)
        #   j in [0..ny] for Y fluxes   (dir = Y)
        return dir * self._numXFluxes + i * (self.Ny() + dir) + j

    # -- arithmetic operators -------------------------------------------
    #
    # NOTE(port): as in scalar.py, elementwise C++ loops over `_data.Size()`
    # are replaced by direct numpy operations on the whole `_data` array.

    def __iadd__(self, other: Union["Flux", float]) -> "Flux":
        if isinstance(other, Flux):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data += other._data
        else:
            self._data += other
        return self

    def __isub__(self, other: Union["Flux", float]) -> "Flux":
        if isinstance(other, Flux):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data -= other._data
        else:
            self._data -= other
        return self

    def __imul__(self, a: float) -> "Flux":
        self._data *= a
        return self

    def __itruediv__(self, a: float) -> "Flux":
        self._data /= a
        return self

    def __add__(self, other: Union["Flux", float]) -> "Flux":
        g = Flux(self)
        g += other
        return g

    def __sub__(self, other: Union["Flux", float]) -> "Flux":
        g = Flux(self)
        g -= other
        return g

    def __mul__(self, a: float) -> "Flux":
        g = Flux(self)
        g *= a
        return g

    def __truediv__(self, a: float) -> "Flux":
        g = Flux(self)
        g /= a
        return g

    def __neg__(self) -> "Flux":
        g = Flux(self)
        g *= -1
        return g

    def __radd__(self, a: float) -> "Flux":
        """a + f"""
        g = Flux(self)
        g += a
        return g

    def __rsub__(self, a: float) -> "Flux":
        """a - f"""
        g = -self
        g += a
        return g

    def __rmul__(self, a: float) -> "Flux":
        """a * f"""
        g = Flux(self)
        g *= a
        return g

    @staticmethod
    def UniformFlow(grid: Grid, magnitude: float, angle: float) -> "Flux":
        """Return Flux for a uniform flow with the specified magnitude and
        direction (angle, in radians)."""
        u = magnitude * math.cos(angle)
        v = magnitude * math.sin(angle)

        q = Flux(grid)
        # NOTE(port): the C++ loop `for (ind = q.begin(X); ind != q.end(X);
        # ++ind) q(lev,ind) = u*dx;` assigns a single constant to a
        # contiguous range of flat indices -- vectorized here as a slice
        # assignment instead of a manual per-index Python loop. The level
        # loop is kept (Ngrid() is small, and Dx(lev) differs per level).
        for lev in range(grid.Ngrid()):
            dx = grid.Dx(lev)
            q._data[lev, q.begin(Direction.X):q.end(Direction.X)] = u * dx
            q._data[lev, q.begin(Direction.Y):q.end(Direction.Y)] = v * dx
        return q

    def setFlow(self, g: TangentSE2, xCenter: float, yCenter: float) -> None:
        """Compute the unsteady base flow motion corresponding to a moving
        body with motion g in TSE(2), about center of rotation
        (xCenter, yCenter)."""
        xdot, ydot, thetadot = g.getVelocity()
        nx = self.Nx()
        ny = self.Ny()

        # NOTE(port): The C++ loop calls x(lev,ind)/y(lev,ind) once per
        # flat index, which internally decomposes ind -> (dir,i,j) and then
        # calls getXEdge/getXCenter/getYEdge/getYCenter (each an O(1)
        # affine function of i or j alone). Rather than replicate that with
        # a Python loop calling x()/y() per index (a manual loop where
        # vectorization is equivalent), we directly build the 2D grids of
        # x- and y-coordinates for the X-flux and Y-flux sub-arrays using
        # the same affine formulas (offset + index * Dx(lev)), then
        # evaluate the (xdot - thetadot*ydiff)*dx / (ydot + thetadot*xdiff)*dx
        # formulas over the whole sub-array at once, and finally flatten
        # back into the flat (dir,i,j) index layout used by `_data`, which
        # is index = i*(Ny()+dir) + j -- i.e. row-major over (i,j) exactly
        # like `_data[lev, begin(dir):end(dir)].reshape(nx_dir, ny_dir)`.
        for lev in range(self.Ngrid()):
            dx = self.Dx(lev)

            # X fluxes: i in 0..nx, j in 0..ny-1
            i_x = np.arange(0, nx + 1)
            j_x = np.arange(0, ny)
            x_x = self.getXEdge(lev, 0) + i_x * dx
            y_x = self.getYCenter(lev, 0) + j_x * dx
            X_i, Y_j = np.meshgrid(x_x - xCenter, y_x - yCenter, indexing="ij")
            xflux_vals = (xdot - thetadot * Y_j) * dx
            self._data[lev, self.begin(Direction.X):self.end(Direction.X)] = xflux_vals.reshape(-1)

            # Y fluxes: i in 0..nx-1, j in 0..ny
            i_y = np.arange(0, nx)
            j_y = np.arange(0, ny + 1)
            x_y = self.getXCenter(lev, 0) + i_y * dx
            y_y = self.getYEdge(lev, 0) + j_y * dx
            X_i2, Y_j2 = np.meshgrid(x_y - xCenter, y_y - yCenter, indexing="ij")
            yflux_vals = (ydot + thetadot * X_i2) * dx
            self._data[lev, self.begin(Direction.Y):self.end(Direction.Y)] = yflux_vals.reshape(-1)
