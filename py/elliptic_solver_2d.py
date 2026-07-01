# elliptic_solver_2d.py
#
# Python port of src/EllipticSolver2d.h / src/EllipticSolver2d.cc
#
# Class for solving Poisson and Helmholtz equations on a uniform grid.
#
# Uses a sin transform, or a nested family of sin transforms. Note that
# EllipticSolver2d is an abstract base class and cannot be instantiated.
#
# Original author: Clancy Rowley (17 Sep 2008)
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES (see also inline NOTE(port) comments):
#
#   1. FFTW -> scipy.fft.  The C++ code uses FFTW's real-to-real transform
#      FFTW_RODFT00, which is the type-I discrete sine transform (DST-I).
#      numpy has no DST, so the exact equivalent scipy.fft.dstn(..., type=1)
#      is used instead. scipy's DST-I uses the same (unnormalized) convention
#      as FFTW_RODFT00, so the normalization factor below is unchanged. This
#      is the one dependency beyond numpy; it is isolated in sinTransform() so
#      it is easy to swap later.
#
#   2. FFTW plan machinery dropped.  The C++ constructor builds an
#      fftw_plan (with FFTW_EXHAUSTIVE) and keeps a scratch array `_fft`;
#      sinTransform copies in/out of that scratch buffer. scipy.fft plans
#      internally and is functional, so `_FFTWPlan` and `_fft` have no Python
#      counterpart and are omitted. The destructor (which only destroyed the
#      plan) is therefore also omitted.
#
#   3. Array2d indexing.  C++ uses Array::Array2<double> with offsets (1,1),
#      i.e. valid indices are i in 1..nx-1, j in 1..ny-1. Here an Array2d is a
#      plain numpy array of shape (nx-1, ny-1) with 0-based indexing, so the
#      C++ element (i,j) is stored at [i-1, j-1]. All index arithmetic below
#      is shifted accordingly.
#
#   4. Overload collapse.  C++ overloads solve() (2-arg zero-BC form and
#      3-arg BC form). Python has no overloading, so both are collapsed into
#      a single solve(f, arg, u=None) that dispatches on whether the third
#      argument is supplied, mirroring the collapse pattern used in the other
#      ported files (Grid, Scalar, Flux).
#
#   5. Output-parameter semantics preserved.  C++ solve()/getRHS() write
#      their result into a caller-supplied output array (pass-by-reference).
#      The Python versions write into `u`/`rhs` in place (u[...] = ...) to
#      keep the same interface, and also return the array for convenience.
#      Aliasing (e.g. f is u, or rhs is u) is handled exactly as in C++.
#
#   6. Abstract base class.  C++ makes EllipticSolver2d abstract via a pure
#      virtual destructor and pure virtual getRHS(). Here it derives from
#      abc.ABC and getRHS is an @abstractmethod, so EllipticSolver2d cannot be
#      instantiated directly and getRHS must be provided by a subclass.
#
#   7. JAX-readiness.  All array work uses basic numpy operations (arange,
#      cos, elementwise arithmetic, slicing) that have jax.numpy equivalents.
#      The one exception is scipy.fft.dstn in sinTransform(); jax.scipy.fft
#      does not currently expose a DST, so that single call is the place that
#      will need attention when porting to JAX. It is deliberately kept in one
#      small method to make that swap local.
# ---------------------------------------------------------------------------

from __future__ import annotations

import abc

import numpy as np
import scipy.fft

from .bc import BC

# NOTE(port): the C++ typedef `Array::Array2<double> Array2d` becomes a plain
# numpy ndarray here (see port note 3). Kept as a module-level alias so the
# array type is defined in one place, matching the FloatArray convention used
# in the other ported files; swap for jax.numpy when porting to JAX.
Array2d = np.ndarray


class EllipticSolver2d(abc.ABC):
    """Abstract base class for solving Poisson and Helmholtz equations on a
    uniform grid, using a discrete sine transform."""

    def __init__(self, nx: int, ny: int, dx: float) -> None:
        # Need only interior points, so eigenvalues are (nx-1) by (ny-1).
        # NOTE(port): the C++ Array2d has offsets (1,1); here it is a 0-based
        # numpy array of shape (nx-1, ny-1). See also port note 2 for the
        # dropped FFTW plan/scratch members.
        self._eigenvaluesOfInverse: Array2d = np.zeros((nx - 1, ny - 1), dtype=np.float64)
        self._nx: int = nx
        self._ny: int = ny
        self._dx: float = dx

    def getLaplacianEigenvalues(self) -> Array2d:
        # calculate eigenvalues of Laplacian
        pi = 4.0 * np.arctan(1.0)
        bydx2 = 1.0 / (self._dx * self._dx)
        # i in 1..nx-1 (column vector), j in 1..ny-1 (row vector)
        i = np.arange(1, self._nx, dtype=np.float64).reshape(-1, 1)
        j = np.arange(1, self._ny, dtype=np.float64).reshape(1, -1)
        eig = 2.0 * (np.cos((pi * i) / self._nx) + np.cos((pi * j) / self._ny) - 2.0) * bydx2
        return eig

    # Take discrete sin transform of u, leaving result in v
    def sinTransform(self, u: Array2d, v: Array2d) -> Array2d:
        assert u.size == v.size
        # NOTE(port): FFTW_RODFT00 (DST-I) -> scipy.fft.dstn(type=1); see
        # port notes 1 and 2. This replaces the copy-in / fftw_execute /
        # copy-out sequence of the C++ code.
        v[...] = scipy.fft.dstn(u, type=1)
        return v

    # Take inverse sin transform of u, leaving result in v.
    # Note that inverse is the same as the forward transform, except for a
    # normalization factor.
    def sinTransformInv(self, u: Array2d, v: Array2d) -> Array2d:
        assert u.size == v.size
        self.sinTransform(u, v)
        normalizationFactor = 1.0 / (2 * self._nx * 2 * self._ny)
        v *= normalizationFactor
        return v

    def solve(self, f: Array2d, arg: object = None, u: Array2d = None) -> Array2d:
        # NOTE(port): overload collapse (port note 4).
        #   solve(f, u)        -> zero-BC form: `arg` is the output array u
        #   solve(f, bc, u)    -> BC form: `arg` is the BC, `u` is the output
        if u is None:
            # Zero-BC form: solve L u = f, assuming zero boundary conditions
            return self._solveZeroBC(f, arg)
        else:
            # BC form: solve L u = f, with specified boundary conditions on u
            return self._solveWithBC(f, arg, u)

    # Solve L u = f, single domain, assuming zero boundary conditions on u
    def _solveZeroBC(self, f: Array2d, u: Array2d) -> Array2d:
        self.sinTransform(f, u)
        u *= self._eigenvaluesOfInverse
        self.sinTransformInv(u, u)  # normalize on inverse transform
        return u

    # Solve L u = f, with specified boundary conditions on u.
    # Note that u contains only the interior points of the domain.
    def _solveWithBC(self, f: Array2d, bc: BC, u: Array2d) -> Array2d:
        rhs = u  # use u as storage for rhs of Poisson equation
        self.getRHS(f, bc, rhs)
        self._solveZeroBC(rhs, u)
        return u

    @abc.abstractmethod
    def getRHS(self, f: Array2d, bc: BC, rhs: Array2d) -> None:
        ...


# -----------------------------------------------------------------------------
# Poisson solver
# -----------------------------------------------------------------------------

class PoissonSolver2d(EllipticSolver2d):
    """Solve a Poisson equation L u = f on a uniform grid, with specified
    boundary conditions on u, where L is the Laplacian."""

    def __init__(self, nx: int, ny: int, dx: float) -> None:
        super().__init__(nx, ny, dx)
        self.calculateEigenvalues()

    def calculateEigenvalues(self) -> None:
        eigL = self.getLaplacianEigenvalues()
        self._eigenvaluesOfInverse[...] = 1.0 / eigL

    # Set rhs = f - L * bc, in preparation for Poisson solve
    def getRHS(self, f: Array2d, bc: BC, rhs: Array2d) -> None:
        assert f.shape[0] == rhs.shape[0]  # f.Nx() == rhs.Nx()
        assert f.shape[1] == rhs.shape[1]  # f.Ny() == rhs.Ny()
        # if input and output arrays are not the same, copy data to rhs
        if rhs is not f:
            rhs[...] = f
        # subtract L(bc) from rhs
        byDx2 = 1.0 / (self._dx * self._dx)
        _subtractBoundary(rhs, bc, byDx2, self._nx, self._ny)


# -----------------------------------------------------------------------------
# Helmholtz solver
# -----------------------------------------------------------------------------

class HelmholtzSolver2d(EllipticSolver2d):
    """Solve a Helmholtz equation (1 + alpha * L) u = f on a uniform grid,
    with specified boundary conditions on u, where L is the Laplacian."""

    def __init__(self, nx: int, ny: int, dx: float, alpha: float) -> None:
        super().__init__(nx, ny, dx)
        self._alpha: float = alpha
        self.calculateEigenvalues()

    def calculateEigenvalues(self) -> None:
        eigL = self.getLaplacianEigenvalues()
        self._eigenvaluesOfInverse[...] = 1.0 / (1 + self._alpha * eigL)

    # Set rhs = f - alpha * L * bc, in preparation for Helmholtz solve
    def getRHS(self, f: Array2d, bc: BC, rhs: Array2d) -> None:
        assert f.shape[0] == rhs.shape[0]  # f.Nx() == rhs.Nx()
        assert f.shape[1] == rhs.shape[1]  # f.Ny() == rhs.Ny()
        # if input and output arrays are not the same, copy data to rhs
        if rhs is not f:
            rhs[...] = f
        # subtract alpha * L(bc) from rhs
        alphaByDx2 = self._alpha / (self._dx * self._dx)
        _subtractBoundary(rhs, bc, alphaByDx2, self._nx, self._ny)


# -----------------------------------------------------------------------------
# Shared helper
# -----------------------------------------------------------------------------

def _subtractBoundary(rhs: Array2d, bc: BC, scale: float, nx: int, ny: int) -> None:
    """Subtract the (scaled) boundary contribution of the discrete Laplacian
    from rhs, reproducing the two boundary loops shared by the Poisson and
    Helmholtz getRHS() implementations.

    The C++ code (for i in 1..nx-1, j in 1..ny-1) is:
        rhs(i,1)    -= bc.bottom(i) * scale
        rhs(i,ny-1) -= bc.top(i)    * scale
        rhs(1,j)    -= bc.left(j)   * scale
        rhs(nx-1,j) -= bc.right(j)  * scale
    which in 0-based storage maps element (i,j) -> [i-1, j-1] (port note 3).
    """
    # NOTE(port): BC exposes only scalar accessors (bottom/top/left/right),
    # so the per-boundary values are gathered into numpy arrays first; the
    # subtractions themselves are then vectorized (no manual loop over the
    # interior). The corner points (e.g. [0,0]) are touched by both a
    # horizontal and a vertical boundary, exactly as in the C++ double loop;
    # subtraction commutes, so the order of the two updates does not matter.
    bottom = np.array([bc.bottom(i) for i in range(1, nx)], dtype=np.float64)
    top = np.array([bc.top(i) for i in range(1, nx)], dtype=np.float64)
    left = np.array([bc.left(j) for j in range(1, ny)], dtype=np.float64)
    right = np.array([bc.right(j) for j in range(1, ny)], dtype=np.float64)

    rhs[:, 0] -= bottom * scale          # rhs(i,1),    i = 1..nx-1
    rhs[:, ny - 2] -= top * scale        # rhs(i,ny-1), i = 1..nx-1
    rhs[0, :] -= left * scale            # rhs(1,j),    j = 1..ny-1
    rhs[nx - 2, :] -= right * scale      # rhs(nx-1,j), j = 1..ny-1
