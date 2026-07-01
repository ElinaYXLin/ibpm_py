# state.py
#
# Python port of src/State.h / src/State.cc
#
# Structure for grouping state variables.

from __future__ import annotations

import sys
from typing import BinaryIO, Optional

import numpy as np

from .boundary_vector import BoundaryVector
from .direction import Direction
from .flux import Flux
from .grid import Grid
from .scalar import Scalar

# NOTE(port): C++ reads/writes raw `int`/`double` values with fread/fwrite,
# using the platform's native size and byte order (no explicit endianness or
# width is specified in State.cc). `np.intc`/`np.float64` with default
# (native, '=') byte order are the numpy equivalents of C `int`/`double` on
# the same machine; like the original code, this format is not portable
# across machines with different native int width or byte order.
_INT_DTYPE = np.dtype(np.intc)
_DOUBLE_DTYPE = np.dtype(np.float64)


def _fread(fp: BinaryIO, dtype: np.dtype, count: int) -> Optional[np.ndarray]:
    """Read `count` values of the given dtype from `fp`.

    Returns None (rather than a short array) if fewer than `count` values
    could be read, mirroring `fread()` returning less than the requested
    item count on a short read/EOF.
    """
    arr = np.fromfile(fp, dtype=dtype, count=count)
    if arr.size != count:
        return None
    return arr


class State:
    """Structure for grouping state variables."""

    def __init__(self, grid: Optional[Grid] = None, numPoints: Optional[int] = None, filename: Optional[str] = None) -> None:
        # NOTE(port): C++ has three constructors: State() (default, no
        # allocation), State(const Grid&, int numPoints) (allocate), and
        # State(string filename) (load from file). Collapsed here into one
        # constructor dispatching on which optional arguments are given,
        # matching the pattern used throughout this port. `filename` takes
        # priority if given alongside grid/numPoints (which the C++
        # overload set never allows simultaneously anyway).
        self.timestep: int = 0
        self.time: float = 0.0
        # NOTE(port): C++ `q`, `omega`, `f` are public data members,
        # default-constructed (Flux()/Scalar()/BoundaryVector(), i.e. no
        # grid/no data) whenever the constructor body does not call
        # resize()/load(). Mirrored here by always constructing empty
        # instances up front.
        self.q: Flux = Flux()
        self.omega: Scalar = Scalar()
        self.f: BoundaryVector = BoundaryVector()

        if filename is not None:
            self.load(filename)
            return

        if grid is not None:
            self.resize(grid, numPoints)

    def resize(self, grid: Grid, numPoints: int) -> None:
        """Allocate memory, with the specified Grid and number of boundary
        points."""
        self.q.resize(grid)
        self.omega.resize(grid)
        self.f.resize(numPoints)

    def computeNetForce(self) -> "tuple[float, float]":
        """Routine for computing X & Y forces.

        Note that f is actually a body force (force per unit mass), so the
        net force is the integral over the domain. By a property of the
        discrete delta function, this equals a sum of the BoundaryVector
        values times dx^2.
        """
        # NOTE(port): C++ signature is `computeNetForce(double& xforce,
        # double& yforce)` (out-params); Python returns a tuple instead,
        # matching the convention used throughout this port.
        #
        # NOTE(port): the C++ loop summing f(X,i)/f(Y,i) over
        # `i=0..numPoints-1` is vectorized here as a `.sum()` over the
        # corresponding X/Y slices of the underlying flat array, instead of
        # a manual per-point Python loop.
        xforce = self.f._data[self.f.begin(Direction.X):self.f.end(Direction.X)].sum()
        yforce = self.f._data[self.f.begin(Direction.Y):self.f.end(Direction.Y)].sum()
        dx2 = self.omega.Dx() * self.omega.Dx()
        xforce *= dx2
        yforce *= dx2
        return float(xforce), float(yforce)

    def load(self, filename: str) -> bool:
        """Load the state from a file (e.g. as a restart file). Return True
        if successful."""
        sys.stderr.write(f"Reading restart file {filename}...")
        sys.stderr.flush()
        try:
            fp = open(filename, "rb")
        except OSError:
            return False

        with fp:
            # read Grid info
            header = _fread(fp, _INT_DTYPE, 3)
            # NOTE(port) -- judgment call: C++ does not check the return
            # value of these header freads at all (only the later, per-value
            # reads are combined into `success`), so a truncated header
            # would silently leave `nx`/`ny`/`ngrid`/etc. as indeterminate
            # (uninitialized) C++ local variables, and then almost certainly
            # crash or behave unpredictably comparing them against
            # q.Nx()/q.Ny()/etc. There is no well-defined "faithful" value to
            # reproduce for reading past EOF in Python, so this is treated
            # as an unconditional failure (return False) instead -- a
            # strictly safer, but behaviorally different, outcome for this
            # corrupt/truncated-file edge case only.
            if header is None:
                return False
            nx, ny, ngrid = (int(v) for v in header)

            dbl_header = _fread(fp, _DOUBLE_DTYPE, 3)
            if dbl_header is None:
                return False
            dx, x0, y0 = (float(v) for v in dbl_header)

            numPoints_arr = _fread(fp, _INT_DTYPE, 1)
            if numPoints_arr is None:
                return False
            numPoints = int(numPoints_arr[0])

            # check that Grid and Geometry in file match those expected
            success = True
            mismatch = (
                nx != self.q.Nx()
                or ny != self.q.Ny()
                or ngrid != self.q.Ngrid()
                or dx != self.q.Dx()
                or x0 != self.q.getXEdge(0, 0)
                or y0 != self.q.getYEdge(0, 0)
                or numPoints != self.f.getNumPoints()
            )
            if mismatch:
                # If old grid was previously allocated, print a warning and
                # set the return value to false
                if self.q.Nx() > 0:
                    sys.stderr.write("Warning: grids do not match.  Resizing grid.\n")
                    success = False
                newgrid = Grid(nx, ny, ngrid, dx * nx, x0, y0)
                self.resize(newgrid, numPoints)

            # read Flux q
            #
            # NOTE(port): the C++ loop over `qind = q.begin() .. q.end()`
            # reads the entire flat per-level data array in order (begin()
            # with no args is 0, end() with no args is the full flux count),
            # which is exactly `self.q._data[lev, :]` in storage order --
            # read here in one `np.fromfile` call per level instead of a
            # manual per-element Python loop.
            numFluxes = self.q.end()
            for lev in range(self.q.Ngrid()):
                vals = _fread(fp, _DOUBLE_DTYPE, numFluxes)
                if vals is None:
                    success = False
                    break
                self.q._data[lev, :] = vals

            # read Scalar omega
            #
            # NOTE(port): the C++ (i,j) double loop over interior points
            # matches Scalar's own storage order exactly ([i-1,j-1], i outer,
            # j inner, both in increasing order) -- read here as one
            # `np.fromfile` call per level, reshaped to (nx-1, ny-1), instead
            # of a manual per-element Python loop.
            for lev in range(self.q.Ngrid()):
                vals = _fread(fp, _DOUBLE_DTYPE, (nx - 1) * (ny - 1))
                if vals is None:
                    success = False
                    break
                self.omega._data[lev, :, :] = vals.reshape(nx - 1, ny - 1)

            # read BoundaryVector f
            #
            # NOTE(port): the C++ loop reads (fx,fy) interleaved per point,
            # which does *not* match BoundaryVector's internal storage layout
            # (all X components, then all Y components) -- read here as one
            # `np.fromfile` call, reshaped to (numPoints, 2), and then
            # de-interleaved into the X/Y blocks with two slice assignments
            # instead of a manual per-point Python loop.
            vals = _fread(fp, _DOUBLE_DTYPE, 2 * numPoints)
            if vals is None:
                success = False
            else:
                interleaved = vals.reshape(numPoints, 2)
                self.f._data[self.f.begin(Direction.X):self.f.end(Direction.X)] = interleaved[:, 0]
                self.f._data[self.f.begin(Direction.Y):self.f.end(Direction.Y)] = interleaved[:, 1]

            # read timestep and time
            timestep_arr = _fread(fp, _INT_DTYPE, 1)
            if timestep_arr is None:
                success = False
            else:
                self.timestep = int(timestep_arr[0])

            time_arr = _fread(fp, _DOUBLE_DTYPE, 1)
            if time_arr is None:
                success = False
            else:
                self.time = float(time_arr[0])

        sys.stderr.write("done\n")
        return success

    def save(self, filename: str) -> bool:
        """Save the state to a file (e.g. as a restart file). Return True
        if successful.

        WARNING: At this point, the xshift and yshift parameters are not
        saved and are not checked for compatibility when loading. Caution
        should be taken when working with shifted grids. This approach was
        taken to preserve backwards compatibility with previously saved
        binary files. In the future perhaps using HDF5 would prevent such
        problems.
        """
        sys.stderr.write(f"Writing restart file {filename}...")
        sys.stderr.flush()
        try:
            fp = open(filename, "wb")
        except OSError:
            return False

        with fp:
            # write Grid info
            grid = self.q.getGrid()
            nx = grid.Nx()
            ny = grid.Ny()
            ngrid = grid.Ngrid()
            dx = grid.Dx()
            x0 = grid.getXEdge(0, 0)
            y0 = grid.getYEdge(0, 0)

            np.array([nx, ny, ngrid], dtype=_INT_DTYPE).tofile(fp)
            np.array([dx, x0, y0], dtype=_DOUBLE_DTYPE).tofile(fp)

            # write Geometry info
            numPoints = self.f.getNumPoints()
            np.array([numPoints], dtype=_INT_DTYPE).tofile(fp)

            # write Flux q
            for lev in range(self.q.Ngrid()):
                self.q._data[lev, :].astype(_DOUBLE_DTYPE).tofile(fp)

            # write Scalar omega
            for lev in range(self.q.Ngrid()):
                self.omega._data[lev, :, :].astype(_DOUBLE_DTYPE).tofile(fp)

            # write BoundaryVector f
            xvals = self.f._data[self.f.begin(Direction.X):self.f.end(Direction.X)]
            yvals = self.f._data[self.f.begin(Direction.Y):self.f.end(Direction.Y)]
            interleaved = np.empty((numPoints, 2), dtype=_DOUBLE_DTYPE)
            interleaved[:, 0] = xvals
            interleaved[:, 1] = yvals
            interleaved.tofile(fp)

            # write timestep and time
            np.array([self.timestep], dtype=_INT_DTYPE).tofile(fp)
            np.array([self.time], dtype=_DOUBLE_DTYPE).tofile(fp)

        sys.stderr.write("done\n")
        return True
