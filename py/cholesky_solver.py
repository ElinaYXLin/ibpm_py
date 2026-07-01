# cholesky_solver.py
#
# Python port of src/CholeskySolver.h / src/CholeskySolver.cc
#
# Subclass of ProjectionSolver in which the system M f = b is solved directly,
# using a Cholesky factorization M = L L^T (L lower triangular). Assumes M is
# symmetric. The factorization is computed at init() time (or loaded from a
# file), then reused for each solve.
#
# Original author: Clancy Rowley (28 Aug 2008)
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES (see also inline NOTE(port) comments):
#
#   1. Storage.  C++ uses Array::array2<double> _lower (_size x _size) and
#      Array::array1<double> _diagonal (_size). These become numpy arrays of
#      shape (_size, _size) and (_size,), both float64. Like BoundaryVector /
#      Scalar in this port, they are zero-filled at allocation (C++ leaves
#      them uninitialized); only entries the algorithm actually reads matter,
#      so this is behaviorally safe.
#
#   2. Reduction summation order.  The inner `for k` reduction loops in
#      computeFactorization() and Minv() are replaced by numpy dot products
#      (the task requires vectorizing loops with an exact numpy equivalent).
#      The C++ loops accumulate from high k down to low k, whereas np.dot
#      accumulates low-to-high; the two are mathematically identical but may
#      differ in the last bit of floating-point rounding. The outer i (and j)
#      loops are NOT vectorized -- they carry a genuine sequential dependency
#      (each triangular entry depends on previously computed ones), so
#      unrolling them would restructure the algorithm.
#
#   3. Upper-triangle "garbage" writes preserved.  In computeFactorization the
#      C++ inner loop runs j over ALL columns; for j < i it writes
#      _lower(j,i) (an upper-triangle entry) by dividing by _diagonal(i),
#      which has not been computed yet at that point. Those upper-triangle
#      values are never read again (Minv only touches strictly-lower entries).
#      We reproduce the exact loop bounds for faithfulness; the stray division
#      is by the zero-initialized diagonal, yielding inf/nan in unused cells
#      (C++ divides by uninitialized memory instead). np.errstate silences the
#      resulting numpy warnings; the discarded values are irrelevant either
#      way. Flagged because it is a deliberate faithful reproduction of an
#      apparently wasteful C++ loop.
#
#   4. Exact float comparison in load().  C++ compares `alphaBeta_in !=
#      _alphaBeta` with exact double equality; ported verbatim. This is
#      fragile (a file written by a slightly different build could fail to
#      load), but it is the original behavior.
#
#   5. File format / precision.  save() writes with C++ `setprecision(17)` and
#      the default (defaultfloat) float format, which prints up to 17
#      significant digits and strips trailing zeros -- reproduced here with
#      Python's `f"{x:.17g}"`. Stream extraction (`infile >> x`) skips
#      arbitrary whitespace; reproduced by tokenizing the whole file on
#      whitespace and consuming tokens in order.
#
#   6. cerr progress messages preserved, written to sys.stderr (the direct
#      analogue of C++ std::cerr).
#
#   7. JAX-readiness.  The heavy numerics (factorization, triangular solves)
#      are expressed as numpy array ops with only the unavoidable sequential
#      outer loops in Python. A later JAX port would replace those with
#      jax.numpy plus jax.lax control flow; nothing here relies on numpy-only
#      semantics (in-place mutation aside, which is called out at each site).
# ---------------------------------------------------------------------------

from __future__ import annotations

import sys

import numpy as np

from .boundary_vector import BoundaryVector
from .grid import Grid
from .navier_stokes_model import NavierStokesModel
from .projection_solver import ProjectionSolver


class CholeskySolver(ProjectionSolver):
    """Solve M f = b directly via a Cholesky factorization of M."""

    def __init__(
        self, grid: Grid, model: NavierStokesModel, beta: float
    ) -> None:
        # Allocate memory for the Cholesky factorization, but do not compute
        # it.
        super().__init__(grid, model, beta)
        self._numPoints: int = model.getNumPoints()   # points in the geometry
        self._size: int = 2 * self._numPoints          # vector size: 2*numPoints
        # local copy of alpha*beta, as a check when reading factorizations
        # from files
        self._alphaBeta: float = model.getAlpha() * beta
        self._lower: np.ndarray = np.zeros((self._size, self._size), dtype=np.float64)
        self._diagonal: np.ndarray = np.zeros(self._size, dtype=np.float64)
        self._hasBeenInitialized: bool = False

    def init(self) -> None:
        """Compute the Cholesky decomposition of M."""
        # Return if CholeskySolver has already been initialized
        if self._hasBeenInitialized:
            return

        # Build matrix M explicitly, one column at a time
        matrixM = np.zeros((self._size, self._size), dtype=np.float64)
        self.computeMatrixM(matrixM)

        # Compute Cholesky factorization
        self.computeFactorization(matrixM)
        self._hasBeenInitialized = True

    def computeMatrixM(self, matrixM: np.ndarray) -> None:
        """Compute the matrix M to be factored, one column at a time."""
        e = BoundaryVector(self._numPoints)   # basis vector
        x = BoundaryVector(self._numPoints)   # x = M(e)

        print("Computing the matrix for Cholesky factorization...", end="", file=sys.stderr, flush=True)
        for j in range(self._size):
            # Compute j-th column of M
            e.assign(0)         # e = 0
            e.set(j, 1)         # j-th basis vector
            self.M(e, x)        # Compute x = M(e)
            # Copy into matrix M
            # NOTE(port): C++ inner `for i` copy vectorized into a column
            # assignment.
            matrixM[:, j] = x.flatten()
        print("done", file=sys.stderr)

    def computeFactorization(self, matrixM: np.ndarray) -> None:
        """Compute the Cholesky factorization M = L L^T (L lower triangular).

        Preconditions:  matrixM is symmetric.
        Postconditions: _lower holds the strictly lower triangular part of L
                        (no diagonal); _diagonal holds L's diagonal.
        """
        print("Computing Cholesky factorization...", end="", file=sys.stderr, flush=True)
        self._lower = matrixM.copy()
        # NOTE(port): see port note 3 for the errstate; the j < i branch may
        # divide by the not-yet-computed (zero) diagonal, writing unused
        # upper-triangle cells.
        with np.errstate(divide="ignore", invalid="ignore"):
            for i in range(self._size):
                for j in range(self._size):
                    # sum = _lower(i,j) - sum_{k=i-1..0} _lower(i,k)*_lower(j,k)
                    # NOTE(port): inner k-reduction -> np.dot (port note 2).
                    s = self._lower[i, j] - float(
                        np.dot(self._lower[i, :i], self._lower[j, :i])
                    )
                    if i == j:
                        assert s > 0
                        self._diagonal[i] = np.sqrt(s)
                    else:
                        self._lower[j, i] = s / self._diagonal[i]
        print("done", file=sys.stderr)

    def load(self, basename: str) -> bool:
        """Load a Cholesky decomposition from <basename>.cholesky.
        Returns True if successful."""
        filename = basename + ".cholesky"
        print(f"Loading Cholesky factorization from file {filename}...", end="", file=sys.stderr, flush=True)
        try:
            with open(filename, "r") as infile:
                contents = infile.read()
        except OSError:
            print("(failed: could not open file)", file=sys.stderr)
            return False

        # NOTE(port): reproduce C++ stream extraction (`infile >> ...`), which
        # skips arbitrary whitespace, by tokenizing on whitespace and
        # consuming tokens in order (port note 5).
        tokens = contents.split()
        pos = 0

        # read dimension of matrix in file
        n = int(tokens[pos]); pos += 1
        if n != self._size:
            print("(failed: wrong file size)", file=sys.stderr)
            return False

        # read value of alphaBeta in file
        alphaBeta_in = float(tokens[pos]); pos += 1
        # NOTE(port): exact float comparison, ported verbatim (port note 4).
        if alphaBeta_in != self._alphaBeta:
            print("(failed: wrong timestep or Re)", file=sys.stderr)
            return False

        # read in diagonal part
        for i in range(self._size):
            self._diagonal[i] = float(tokens[pos]); pos += 1

        # check the marker, to make sure we did not get off track
        c = tokens[pos]; pos += 1
        if c != "#":
            print("(failed: corrupt file)", file=sys.stderr)
            return False

        # read the lower triangular portion
        for i in range(self._size):
            for j in range(i):
                self._lower[i, j] = float(tokens[pos]); pos += 1
        self._hasBeenInitialized = True
        print("done", file=sys.stderr)
        return True

    def save(self, basename: str) -> bool:
        """Save a Cholesky decomposition to <basename>.cholesky, overwriting
        if necessary. Returns True if successful.

        Saves only the strictly lower triangular portion of _lower, since that
        is all that is needed for back substitution.
        """
        assert self._hasBeenInitialized
        filename = basename + ".cholesky"
        print(f"Saving Cholesky factorization to file {filename}...", end="", file=sys.stderr, flush=True)
        try:
            outfile = open(filename, "w")
        except OSError:
            print("(failed: could not open file)", file=sys.stderr)
            return False
        with outfile:
            outfile.write(f"{self._size}\n")
            # NOTE(port): setprecision(17) + defaultfloat -> f"{x:.17g}"
            # (port note 5).
            outfile.write(f"{self._alphaBeta:.17g}\n")
            # write the diagonal part
            for i in range(self._size):
                outfile.write(f"{self._diagonal[i]:.17g}\n")

            # insert marker, as a crude verification that we did not skip a
            # line when reading back in
            outfile.write("#\n")

            # write the lower triangular part
            for i in range(self._size):
                for j in range(i):
                    outfile.write(f"{self._lower[i, j]:.17g}\n")
        print("done", file=sys.stderr)
        return True

    def Minv(self, b: BoundaryVector, x: BoundaryVector) -> None:
        """Solve M x = b using the Cholesky factorization M = L L^T.
        Assumes M is symmetric.

        NOTE(port): `x` is an out-parameter, written in place.
        """
        assert self._hasBeenInitialized
        bd = b.flatten()
        xd = x.flatten()

        # Solve L y = b for y
        # (Here, y and x use the same memory, for efficiency)
        for i in range(self._size):
            # sum = b(i) - sum_{k=i-1..0} _lower(i,k)*x(k)
            # NOTE(port): inner k-reduction -> np.dot (port note 2).
            s = bd[i] - float(np.dot(self._lower[i, :i], xd[:i]))
            xd[i] = s / self._diagonal[i]

        # Solve L^T x = y for x
        for i in range(self._size - 1, -1, -1):
            # sum = x(i) - sum_{k=i+1..size-1} _lower(k,i)*x(k)
            s = xd[i] - float(np.dot(self._lower[i + 1 : self._size, i], xd[i + 1 : self._size]))
            xd[i] = s / self._diagonal[i]
