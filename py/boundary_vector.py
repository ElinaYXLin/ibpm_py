# boundary_vector.py
#
# Python port of src/BoundaryVector.h / src/BoundaryVector.cc
#
# Store vectors located at boundary points.
#
# Examples are forces, velocities, and coordinates of points on a body.

from __future__ import annotations

from typing import Optional, Union

import numpy as np

from .direction import Direction

FloatArray = np.ndarray


class BoundaryVector:
    """Vectors located at boundary points (e.g. forces, velocities,
    coordinates), stored as a flat array of length 2*numPoints, with the
    X-components first, then the Y-components."""

    def __init__(self, arg: Optional[Union[int, "BoundaryVector"]] = None) -> None:
        # NOTE(port): C++ has three constructors: BoundaryVector()
        # (default, no allocation), BoundaryVector(int numPoints)
        # (allocate), and BoundaryVector(const BoundaryVector&) (copy).
        # Collapsed here into one constructor dispatching on the type of
        # `arg`: None reproduces the default constructor, an int
        # reproduces BoundaryVector(int), and a BoundaryVector reproduces
        # the copy constructor.
        if arg is None:
            self._numPoints: int = 0
            self._data: Optional[FloatArray] = None
            return
        if isinstance(arg, BoundaryVector):
            self.resize(arg._numPoints)
            self._data[...] = arg._data
            return
        self.resize(arg)

    def resize(self, numPoints: int) -> None:
        """Reallocate memory for the given number of points."""
        self._numPoints = numPoints
        # NOTE(port): C++ `_data.Allocate(...)` leaves memory uninitialized
        # (see the equivalent note in scalar.py); zero-filled here instead
        # for deterministic behavior.
        self._data = np.zeros(numPoints * int(Direction.XY), dtype=np.float64)

    def getNumPoints(self) -> int:
        """Return the number of boundary points."""
        return self._numPoints

    def getSize(self) -> int:
        """Return the number of elements in the array (twice the number of
        boundary points)."""
        return int(Direction.XY) * self._numPoints

    def print(self) -> None:
        """Print the contents to standard output, for debugging."""
        print(self._data)

    def __str__(self) -> str:
        # NOTE(port): reproduces C++ `operator<<(ostream&, const
        # BoundaryVector&)`, which Python spells as `__str__` /
        # `print(bv)` rather than `cout << bv`.
        lines = [f"  ({self(Direction.X, i)}, {self(Direction.Y, i)})" for i in range(self._numPoints)]
        return "\n".join(lines) + ("\n" if lines else "")

    def __call__(self, dir_or_ind, i: Optional[int] = None) -> float:
        """f(dir,i) refers to the value in direction dir (X or Y) at point
        i, or f(ind) refers to the value at the flat index ind.

        NOTE(port): as in scalar.py/flux.py, C++ `operator()` returning
        `double&` (read-write) is split here into `__call__` (read-only)
        and `set(...)` (write) since Python cannot return a writable
        reference from a function call.
        """
        if i is None:
            ind = dir_or_ind
            assert ind >= 0 and ind < self._numPoints * int(Direction.XY)
            return float(self._data[ind])
        dir_ = dir_or_ind
        assert dir_ >= Direction.X and dir_ <= Direction.Y
        assert i >= 0 and i < self._numPoints
        return float(self._data[int(dir_) * self._numPoints + i])

    def set(self, *args) -> None:
        """Write counterpart to __call__ -- see the NOTE(port) there.

        Call as `set(ind, value)` (flat-index form) or
        `set(dir, i, value)` (direction/point-index form).
        """
        if len(args) == 2:
            ind, value = args
            assert ind >= 0 and ind < self._numPoints * int(Direction.XY)
            self._data[ind] = value
            return
        dir_, i, value = args
        assert dir_ >= Direction.X and dir_ <= Direction.Y
        assert i >= 0 and i < self._numPoints
        self._data[int(dir_) * self._numPoints + i] = value

    def begin(self, dir: Optional[int] = None) -> int:
        """Returns an index that refers to the first element (overall, or
        the first element in direction dir if given)."""
        if dir is None:
            return 0
        assert dir >= Direction.X and dir <= Direction.Y
        return int(dir) * self._numPoints

    def end(self, dir: Optional[int] = None) -> int:
        """Returns an index one past the last element (overall, or one past
        the last element in direction dir if given)."""
        if dir is None:
            return self._numPoints * int(Direction.XY)
        assert dir >= Direction.X and dir <= Direction.Y
        return (int(dir) + 1) * self._numPoints

    def getIndex(self, dir: int, i: int) -> int:
        """Returns an index for the value in direction dir at point i."""
        assert dir >= Direction.X and dir <= Direction.Y
        assert i >= 0 and i < self._numPoints
        return int(dir) * self._numPoints + i

    def flatten(self) -> FloatArray:
        """Return the underlying data, expressed as a flat array.

        NOTE(port): C++ returns `double*`, a raw pointer to the underlying
        buffer (declared but its .cc implementation is not present in this
        codebase snapshot -- only the header declares it). The natural
        Python/numpy equivalent is to return the underlying ndarray
        directly (a view, not a copy), which callers can pass to any numpy
        or C-extension API expecting a flat buffer.
        """
        return self._data

    def dot(self, f: "BoundaryVector") -> float:
        """Return the dot product of *this and the argument.

        NOTE(port): the C++ declaration for `dot()` exists in the header
        but has a `// TODO: Implement this, and write tests` comment next
        to it in BoundaryVector.h, and no definition appears in
        BoundaryVector.cc. There is therefore no C++ behavior to port
        faithfully; this implements the natural vectorized meaning (numpy
        dot product over the flat data), consistent with the free function
        `InnerProduct` below (which *is* implemented in the header).
        """
        return float(np.dot(self._data, f._data))

    def assign(self, other: Union["BoundaryVector", float]) -> "BoundaryVector":
        """Copy assignment, from another BoundaryVector or from a scalar.

        NOTE(port): see Scalar.assign() -- Python cannot overload plain
        `=`, so both C++ `operator=(const BoundaryVector&)` and
        `operator=(double)` are exposed through this one method.
        """
        if isinstance(other, BoundaryVector):
            assert other._numPoints == self._numPoints
            self._data[...] = other._data
        else:
            self._data[...] = other
        return self

    # -- arithmetic operators -------------------------------------------
    #
    # NOTE(port): C++ operates on the whole `_data` array via Blitz++
    # expression templates (`_data += f._data;` etc); the numpy operations
    # below are the direct, already-vectorized equivalent.

    def __iadd__(self, other: Union["BoundaryVector", float]) -> "BoundaryVector":
        if isinstance(other, BoundaryVector):
            assert other._numPoints == self._numPoints
            self._data += other._data
        else:
            self._data += other
        return self

    def __isub__(self, other: Union["BoundaryVector", float]) -> "BoundaryVector":
        if isinstance(other, BoundaryVector):
            assert other._numPoints == self._numPoints
            self._data -= other._data
        else:
            self._data -= other
        return self

    def __imul__(self, a: float) -> "BoundaryVector":
        self._data *= a
        return self

    def __itruediv__(self, a: float) -> "BoundaryVector":
        self._data /= a
        return self

    def __add__(self, other: Union["BoundaryVector", float]) -> "BoundaryVector":
        g = BoundaryVector(self)
        g += other
        return g

    def __sub__(self, other: Union["BoundaryVector", float]) -> "BoundaryVector":
        g = BoundaryVector(self)
        g -= other
        return g

    def __mul__(self, a: float) -> "BoundaryVector":
        g = BoundaryVector(self)
        g *= a
        return g

    def __truediv__(self, a: float) -> "BoundaryVector":
        g = BoundaryVector(self)
        g /= a
        return g

    def __neg__(self) -> "BoundaryVector":
        g = BoundaryVector(self)
        g *= -1
        return g

    def __rmul__(self, a: float) -> "BoundaryVector":
        """a * f"""
        g = BoundaryVector(self)
        g *= a
        return g


def InnerProduct(x: "BoundaryVector", y: "BoundaryVector") -> float:
    """Return the inner product of BoundaryVectors x and y."""
    # NOTE(port): C++ implements this with an explicit loop over indices
    # "using only public interface" (per its own comment), with a
    # commented-out, faster Blitz-array version `sum(x._data * y._data)`.
    # We take the vectorized numpy equivalent of that commented-out
    # version (equivalent to `x._data @ y._data`), since it is
    # mathematically identical and the task requires vectorizing rather
    # than manually looping where numpy vectorization is equivalent.
    return float(np.dot(x._data, y._data))
