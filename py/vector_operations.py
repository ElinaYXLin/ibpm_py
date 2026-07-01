# vector_operations.py
#
# Partial Python port of src/VectorOperations.h / src/VectorOperations.cc
#
# NOTE(port) -- scope judgment call: VectorOperations.{h,cc} declares a large
# collection of free functions (Curl, InnerProduct, FineGridInnerProduct,
# VorticityInnerProduct, CrossProduct, FluxToXVelocity, FluxToYVelocity,
# XVelocityToFlux, YVelocityToFlux, VelocityToFlux, FluxToVelocity, Laplacian).
# Per the porting instructions ("only port other files when there are
# dependencies"), only `Curl` is ported here, since NavierStokesModel is the
# only file being ported that needs anything from this module, and it only
# calls the two `Curl` overloads (Curl(q, omega) in B(), and
# Curl(streamfunction, q) in computeFluxWithoutBaseFlow()). The remaining
# functions declared in VectorOperations.h are not ported; add them here
# (following the same vectorization approach) if/when a dependent file needs
# them.

from __future__ import annotations

from typing import Optional, Union

import numpy as np

from .bc import BC
from .direction import Direction
from .flux import Flux
from .scalar import Scalar

FloatArray = np.ndarray


def _bc_left(bc: BC, j: FloatArray) -> FloatArray:
    """Vectorized equivalent of BC.left(j) for an array of indices j.

    NOTE(port): reads `bc._data` directly (rather than looping over the
    public, scalar `BC.left()` accessor) to vectorize what is a per-element
    C++ loop body in Curl(Scalar,Flux); see the identical convention used
    elsewhere in this port (e.g. regularizer.py, geometry.py) for reaching
    into a sibling class's private array when vectorizing bulk access.
    """
    return bc._data[np.asarray(j)]


def _bc_right(bc: BC, j: FloatArray) -> FloatArray:
    """Vectorized equivalent of BC.right(j); see _bc_left."""
    j = np.asarray(j)
    return bc._data[2 * bc._ny + bc._nx - j]


def _bc_top(bc: BC, i: FloatArray) -> FloatArray:
    """Vectorized equivalent of BC.top(i); see _bc_left."""
    i = np.asarray(i)
    return bc._data[bc._ny + i]


def _bc_bottom(bc: BC, i: FloatArray) -> FloatArray:
    """Vectorized equivalent of BC.bottom(i); see _bc_left.

    Reproduces BC.bottom()'s special case at i == 0 (stored at `_data[0]`
    rather than at the general `2*(nx+ny) - i` formula, which would be out
    of bounds at i == 0) via `np.where`.
    """
    i = np.asarray(i)
    # NOTE(port): compute the index arithmetic first and index bc._data only
    # once with the selected (always in-bounds) index -- np.where would
    # otherwise evaluate `bc._data[2*(nx+ny) - i]` for every element
    # (including i == 0, where that formula is out of bounds) before
    # selecting between the two branches.
    idx = np.where(i == 0, 0, 2 * (bc._nx + bc._ny) - i)
    return bc._data[idx]


def _curl_flux_to_scalar(q: Flux, f: Scalar) -> None:
    """Compute the curl of Flux q, as a Scalar object f (in place)."""
    assert q.Nx() == f.Nx()
    assert q.Ny() == f.Ny()
    assert q.Ngrid() == f.Ngrid()
    nx = q.Nx()
    ny = q.Ny()

    # Curl (u,v) = v_x - u_y
    #
    # NOTE(port): the C++ (i,j) double loop over interior nodes, for each
    # grid level, is fully vectorized below with numpy slicing (no manual
    # per-node Python loop); the outer level loop is kept as a plain Python
    # loop since Ngrid is small and each level's dx differs.
    for lev in range(q.Ngrid()):
        dx = q.Dx(lev)
        bydx2 = 1.0 / (dx * dx)

        X_arr = q._data[lev, q.begin(Direction.X):q.end(Direction.X)].reshape(nx + 1, ny)
        Y_arr = q._data[lev, q.begin(Direction.Y):q.end(Direction.Y)].reshape(nx, ny + 1)

        # f(lev,i,j) for i=1..nx-1, j=1..ny-1 (Scalar's full interior range)
        f._data[lev] = (
            Y_arr[1:nx, 1:ny] - Y_arr[0:nx - 1, 1:ny] + X_arr[1:nx, 0:ny - 1] - X_arr[1:nx, 1:ny]
        ) * bydx2


def _curl_scalar_to_flux(f: Scalar, q: Flux) -> None:
    """Compute the curl of Scalar f, as a Flux object q (in place)."""
    assert f.Nx() == q.Nx()
    assert f.Ny() == q.Ny()
    assert f.Ngrid() == q.Ngrid()
    nx = f.Nx()
    ny = f.Ny()

    # boundary condition object for use in computing curl on finer grids
    bc = BC(nx, ny)

    # From coarsest grid to finest
    #
    # NOTE(port): the level loop has a genuine sequential dependency (each
    # finer level's boundary conditions come from the next coarser level's
    # already-computed Scalar values via f.getBC()), so it is kept as a
    # Python loop; everything inside a given level (all (i,j) loops in the
    # C++ body) is vectorized with numpy slicing instead of manual
    # per-element loops.
    for lev in range(f.Ngrid() - 1, -1, -1):
        # For outermost grid, all boundaries are zero
        if lev == f.Ngrid() - 1:
            bc.assign(0.0)
        # Otherwise, get bc from next coarser grid
        else:
            f.getBC(lev, bc)

        fdata = f._data[lev]  # shape (nx-1, ny-1), index [i-1, j-1]

        X_arr = q._data[lev, q.begin(Direction.X):q.end(Direction.X)].reshape(nx + 1, ny)
        Y_arr = q._data[lev, q.begin(Direction.Y):q.end(Direction.Y)].reshape(nx, ny + 1)

        # X direction: u = df/dy

        # Compute all points except boundaries
        X_arr[1:nx, 1:ny - 1] = fdata[:, 1:ny - 1] - fdata[:, 0:ny - 2]

        # top and bottom boundaries
        i_arr = np.arange(1, nx)
        X_arr[1:nx, 0] = fdata[:, 0] - _bc_bottom(bc, i_arr)
        X_arr[1:nx, ny - 1] = _bc_top(bc, i_arr) - fdata[:, -1]

        # left and right boundaries
        j_arr = np.arange(0, ny)
        X_arr[0, :] = _bc_left(bc, j_arr + 1) - _bc_left(bc, j_arr)
        X_arr[nx, :] = _bc_right(bc, j_arr + 1) - _bc_right(bc, j_arr)

        # Y direction: v = -df/dx

        # Compute all points except boundaries
        Y_arr[1:nx - 1, 1:ny] = fdata[0:nx - 2, :] - fdata[1:nx - 1, :]

        # left and right boundaries
        j_arr2 = np.arange(1, ny)
        Y_arr[0, 1:ny] = _bc_left(bc, j_arr2) - fdata[0, :]
        Y_arr[nx - 1, 1:ny] = fdata[-1, :] - _bc_right(bc, j_arr2)

        # top and bottom boundaries
        i_arr2 = np.arange(0, nx)
        Y_arr[:, 0] = _bc_bottom(bc, i_arr2) - _bc_bottom(bc, i_arr2 + 1)
        Y_arr[:, ny] = _bc_top(bc, i_arr2) - _bc_top(bc, i_arr2 + 1)

        # NOTE(port): X_arr/Y_arr are numpy views obtained via .reshape() of
        # a contiguous slice of q._data, so in-place mutation above already
        # writes through to q._data. The explicit write-back below is
        # nonetheless kept (matching the equivalent explicit-write pattern
        # used in flux.py's setFlow()) so correctness does not silently rely
        # on reshape() returning a view rather than a copy.
        q._data[lev, q.begin(Direction.X):q.end(Direction.X)] = X_arr.reshape(-1)
        q._data[lev, q.begin(Direction.Y):q.end(Direction.Y)] = Y_arr.reshape(-1)


def Curl(arg: Union[Flux, Scalar], out: Optional[Union[Scalar, Flux]] = None) -> Optional[Union[Scalar, Flux]]:
    """Return the curl of a Flux (as a Scalar) or of a Scalar (as a Flux).

    NOTE(port): collapses all four C++ overloads
        Scalar Curl(const Flux& q);       void Curl(const Flux& q, Scalar& omega);
        Flux   Curl(const Scalar& f);     void Curl(const Scalar& f, Flux& q);
    into one function dispatching on the runtime type of `arg` (Flux vs.
    Scalar) and on whether `out` is supplied, matching the overload-collapse
    convention used throughout this port (see e.g. EllipticSolver.solve()).
    When `out` is given, it is mutated in place and `None` is returned,
    matching the C++ `void` overloads; when omitted, a new object is
    allocated, filled, and returned, matching the C++ value-returning
    overloads.
    """
    if isinstance(arg, Flux):
        q = arg
        if out is None:
            omega = Scalar(q.getGrid())
            _curl_flux_to_scalar(q, omega)
            return omega
        _curl_flux_to_scalar(q, out)
        return None
    elif isinstance(arg, Scalar):
        f = arg
        if out is None:
            q = Flux(f.getGrid())
            _curl_scalar_to_flux(f, q)
            return q
        _curl_scalar_to_flux(f, out)
        return None
    else:
        raise TypeError(f"Curl: unsupported argument type {type(arg)!r}")
