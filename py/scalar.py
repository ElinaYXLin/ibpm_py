# scalar.py
#
# Python port of src/Scalar.h / src/Scalar.cc
#
# Store a 2D array of scalar values, located at interior cell nodes.
#
# For a grid with nx cells in the x-direction and ny cells in the
# y-direction, there are (nx+1)*(ny+1) nodes.
#
# There are (nx-1)*(ny-1) inner nodes, and 2*(nx+ny) boundary nodes.
# Only the interior nodes are stored in a Scalar, and the boundary nodes are
# always zero.

from __future__ import annotations

from typing import Optional, Union

import numpy as np

from .bc import BC
from .field import Field
from .grid import Grid

# NOTE(port): kept as a module-level type alias so that the ndarray dtype
# used for all Scalar/Flux/BC data is defined in exactly one place. When
# porting to JAX later, swap this alias (and the corresponding array
# constructors) for `jax.numpy` equivalents; nothing else in this class
# should need to change, since indexing here uses only basic (start:stop:
# step) numpy slicing, which JAX arrays also support.
FloatArray = np.ndarray


class Scalar(Field):
    """A 2D array of scalar values, located at interior cell nodes."""

    def __init__(self, arg: Optional[Union[Grid, "Scalar"]] = None) -> None:
        # NOTE(port): C++ has three constructors:
        #   Scalar(const Grid& grid)   -- allocate memory for the given grid
        #   Scalar()                    -- default, do not allocate memory
        #   Scalar(const Scalar& f)     -- allocate new array, copy the data
        # Collapsed here into one constructor dispatching on the runtime
        # type of `arg`: None reproduces Scalar(), a Grid reproduces
        # Scalar(const Grid&), and a Scalar reproduces the copy
        # constructor Scalar(const Scalar&).
        if arg is None:
            super().__init__()
            self._data: Optional[FloatArray] = None
            return
        if isinstance(arg, Scalar):
            super().__init__(arg.getGrid())
            self.resize(arg.getGrid())
            self._data[...] = arg._data
            return
        super().__init__(arg)
        self.resize(arg)

    def resize(self, grid: Grid) -> None:
        """Reassign the grid parameters and allocate memory based on the
        new grid."""
        self.setGrid(grid)
        # Allocate arrays for interior points:
        #    lev in 0..Ngrid-1
        #    i   in 1..Nx-1   (stored at array index i-1)
        #    j   in 1..Ny-1   (stored at array index j-1)
        #
        # NOTE(port): C++ `_data.Allocate(...)` (Blitz-style) allocates
        # uninitialized memory (no zero-fill), matching `malloc`/`new[]`
        # semantics. We deliberately zero-fill here with `np.zeros` instead
        # of `np.empty`, trading strict bit-for-bit fidelity with C++'s
        # uninitialized memory for reproducible, deterministic behavior
        # (relying on uninitialized memory is undefined behavior in C++
        # too, so there is no "faithful" value to copy here).
        self._data = np.zeros((self.Ngrid(), self.Nx() - 1, self.Ny() - 1), dtype=np.float64)

    def print(self) -> None:
        """Print the whole field to standard output."""
        nx = self.Nx()
        ny = self.Ny()
        for lev in range(self.Ngrid()):
            for j in range(ny - 1, 0, -1):
                row = " ".join(str(self._data[lev, i - 1, j - 1]) for i in range(1, nx))
                print(row)
            print()

    def coarsify(self) -> None:
        """"Coarsify" the Scalar quantity.
        - Fine grid is left unchanged
        - Coarse values that correspond to points on the fine grid are
          replaced by averaged value of fine gridpoints
        """
        # NOTE(port): The C++ triple loop (lev, i, j) has a genuine
        # sequential dependency across `lev` (level `lev` reads level
        # `lev-1`), so the outer loop over levels cannot be vectorized and
        # is kept as a Python loop (this loop only runs Ngrid-1 times, a
        # small constant, not over grid points). The inner (i, j) double
        # loop, however, is a fixed-stencil average with no
        # cross-iteration dependency, so it is fully vectorized with numpy
        # fancy indexing below instead of a manual element-by-element loop.
        nxExt = self.NxExt()
        nyExt = self.NyExt()
        for lev in range(1, self.Ngrid()):
            # Loop over interior gridpoints that correspond to the finer grid
            i_vals = np.arange(nxExt + 1, self.Nx() // 2 + nxExt)
            j_vals = np.arange(nyExt + 1, self.Ny() // 2 + nyExt)
            if i_vals.size == 0 or j_vals.size == 0:
                continue
            # Corresponding points on the fine grid: Grid.c2f(i,j) -> (ii,jj)
            ii_vals = (i_vals - nxExt) * 2
            jj_vals = (j_vals - nyExt) * 2

            di, dj = np.meshgrid(i_vals - 1, j_vals - 1, indexing="ij")
            fii, fjj = np.meshgrid(ii_vals - 1, jj_vals - 1, indexing="ij")

            fine = self._data[lev - 1]
            self._data[lev, di, dj] = (
                0.25 * fine[fii, fjj]
                + 0.125 * (fine[fii + 1, fjj] + fine[fii, fjj + 1] + fine[fii - 1, fjj] + fine[fii, fjj - 1])
                + 0.0625 * (
                    fine[fii + 1, fjj + 1] + fine[fii + 1, fjj - 1]
                    + fine[fii - 1, fjj + 1] + fine[fii - 1, fjj - 1]
                )
            )

    def __call__(self, lev: int, i: int, j: int) -> float:
        """f(lev, i, j) refers to the value at index (i,j)."""
        # NOTE(port): C++ `operator()(int,int,int)` returns `double&`, so it
        # is usable both to read and, via assignment, to write
        # (`f(lev,i,j) = value`). Python has no equivalent of returning a
        # reference from a function call, so this `__call__` only supports
        # reads; use `set(lev, i, j, value)` for writes (see below).
        assert lev >= 0 and lev < self.Ngrid()
        assert i >= 1 and i < self.Nx()
        assert j >= 1 and j < self.Ny()
        return float(self._data[lev, i - 1, j - 1])

    def set(self, lev: int, i: int, j: int, value: float) -> None:
        """Set the value at index (lev, i, j). Substitute for the write
        side of the C++ `double& operator()(int,int,int)`, which Python
        cannot express (see __call__ above)."""
        assert lev >= 0 and lev < self.Ngrid()
        assert i >= 1 and i < self.Nx()
        assert j >= 1 and j < self.Ny()
        self._data[lev, i - 1, j - 1] = value

    def __getitem__(self, lev: int) -> FloatArray:
        """f[lev] returns a 2d array of grid level lev.

        Returned array has shape (Nx()-1, Ny()-1), and is indexed as
        arr[i-1, j-1] for interior point (i,j) -- i.e. index 0 corresponds
        to i=1 / j=1, matching the storage offset in the underlying C++
        Blitz array. It is a view, not a copy, so in-place mutation
        (`f[lev][:] = ...`) writes through to the Scalar, mirroring the
        mutability of the C++ `Array::Array2<double>` returned by
        `operator[]`.
        """
        return self._data[lev]

    def assign(self, other: Union["Scalar", float]) -> "Scalar":
        """Copy assignment, from another Scalar or from a scalar value.

        NOTE(port): C++ overloads `operator=` for both `const Scalar&` and
        `double`, so `f = g` and `f = a` both work directly. Python cannot
        overload plain `=` (it only rebinds the name), so both cases are
        exposed through this single `assign()` method instead, dispatching
        on the argument type exactly as the two C++ overloads did.
        """
        if isinstance(other, Scalar):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data[...] = other._data
        else:
            self._data[...] = other
        return self

    def getBC(self, lev: int, bc: BC) -> None:
        """Compute the boundary values at level lev from the next coarser
        grid.

        lev: the grid level for which the boundary values are desired;
             must be in the range 0..Ngrid-2
        bc:  BC object that receives the boundary values (mutated in place)
        """
        assert self.Nx() == bc.Nx()
        assert self.Ny() == bc.Ny()
        assert lev >= 0 and lev < self.Ngrid() - 1

        # NOTE(port): The C++ loops here use the idiom
        #   for (int i=0; i<=Nx(); ++i) { ...; if (++i <= Nx()) {...} }
        # which relies on an extra `++i` inside the loop body to advance
        # the loop variable by 2 per iteration: even i are copied directly
        # from the coincident coarse point, odd i are interpolated between
        # the two coarse points on either side. That idiom has no direct,
        # vectorizable Python translation (mutating the loop variable
        # inside a `for i in range(...)` body has no effect in Python), so
        # it is re-expressed below as two numpy strided-slice assignments
        # (even-index direct copy, odd-index averaged interpolation),
        # which is behaviorally identical to the C++ loop for the even
        # Nx()/Ny() enforced by the class invariant (nx, ny multiples of 4
        # whenever Ngrid() > 1, which is required for getBC to be called).
        nx = self.Nx()
        ny = self.Ny()
        nxExt = self.NxExt()
        nyExt = self.NyExt()

        # top and bottom boundaries
        # if grid is shifted completely up or down, then all points on the
        # shared boundary must take a value of 0, as required on the
        # boundary of the outermost grid
        k = np.arange(0, nx // 2 + 1)  # i = 2*k
        ii = k + nxExt  # from Grid.f2c(i, 0) with i = 2*k
        jj = nyExt  # constant: Grid.f2c(i, 0) -> jj = NyExt() for all even i

        bottom_even = np.array([self(lev + 1, int(iv), jj) for iv in ii])
        top_even = np.array([self(lev + 1, int(iv), ny // 2 + jj) for iv in ii])

        bottom = np.zeros(nx + 1, dtype=np.float64)
        top = np.zeros(nx + 1, dtype=np.float64)
        bottom[0::2] = bottom_even
        top[0::2] = top_even
        bottom[1::2] = 0.5 * (bottom_even[:-1] + bottom_even[1:])
        top[1::2] = 0.5 * (top_even[:-1] + top_even[1:])
        for i in range(nx + 1):
            bc.setBottom(i, float(bottom[i]))
            bc.setTop(i, float(top[i]))

        # left and right boundaries
        m = np.arange(0, ny // 2 + 1)  # j = 2*m
        jj2 = m + nyExt  # from Grid.f2c(0, j) with j = 2*m
        ii2 = nxExt  # constant: Grid.f2c(0, j) -> ii = NxExt() for all even j

        left_even = np.array([self(lev + 1, ii2, int(jv)) for jv in jj2])
        right_even = np.array([self(lev + 1, nx // 2 + ii2, int(jv)) for jv in jj2])

        left = np.zeros(ny + 1, dtype=np.float64)
        right = np.zeros(ny + 1, dtype=np.float64)
        left[0::2] = left_even
        right[0::2] = right_even
        left[1::2] = 0.5 * (left_even[:-1] + left_even[1:])
        right[1::2] = 0.5 * (right_even[:-1] + right_even[1:])
        for j in range(ny + 1):
            bc.setLeft(j, float(left[j]))
            bc.setRight(j, float(right[j]))

    # -- arithmetic operators -------------------------------------------
    #
    # NOTE(port): The C++ elementwise loops (`for i<_data.Size(): _data(i)
    # op= ...`) are replaced by direct numpy array operations on the whole
    # `_data` array, which is the vectorized equivalent with no manual
    # Python loop over elements.

    def __iadd__(self, other: Union["Scalar", float]) -> "Scalar":
        if isinstance(other, Scalar):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data += other._data
        else:
            self._data += other
        return self

    def __isub__(self, other: Union["Scalar", float]) -> "Scalar":
        if isinstance(other, Scalar):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data -= other._data
        else:
            self._data -= other
        return self

    def __imul__(self, other: Union["Scalar", float]) -> "Scalar":
        if isinstance(other, Scalar):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data *= other._data
        else:
            self._data *= other
        return self

    def __itruediv__(self, other: Union["Scalar", float]) -> "Scalar":
        if isinstance(other, Scalar):
            assert other.Ngrid() == self.Ngrid()
            assert other.Nx() == self.Nx()
            assert other.Ny() == self.Ny()
            self._data /= other._data
        else:
            self._data /= other
        return self

    def __add__(self, other: Union["Scalar", float]) -> "Scalar":
        g = Scalar(self)
        g += other
        return g

    def __sub__(self, other: Union["Scalar", float]) -> "Scalar":
        g = Scalar(self)
        g -= other
        return g

    def __mul__(self, other: Union["Scalar", float]) -> "Scalar":
        g = Scalar(self)
        g *= other
        return g

    def __truediv__(self, other: Union["Scalar", float]) -> "Scalar":
        g = Scalar(self)
        g /= other
        return g

    def __neg__(self) -> "Scalar":
        g = Scalar(self)
        g *= -1
        return g

    def __radd__(self, a: float) -> "Scalar":
        """a + f"""
        g = Scalar(self)
        g += a
        return g

    def __rsub__(self, a: float) -> "Scalar":
        """a - f"""
        g = -self
        g += a
        return g

    def __rmul__(self, a: float) -> "Scalar":
        """a * f"""
        g = Scalar(self)
        g *= a
        return g

    def __rtruediv__(self, a: float) -> "Scalar":
        """a / f"""
        g = Scalar(self.getGrid())
        g.assign(a)
        g /= self
        return g
