# vector_operations.py
#
# Partial Python port of src/VectorOperations.h / src/VectorOperations.cc
#
# NOTE(port) -- scope: VectorOperations.{h,cc} declares a large collection of
# free functions. This module is ported incrementally, as dependent files
# need pieces of it. Ported so far:
#   * Curl (both directions)            -- needed by NavierStokesModel
#   * Laplacian (Scalar forms)          -- needed by IBSolver
#   * CrossProduct (both forms)         -- needed by IBSolver
#   * FluxToXVelocity / FluxToYVelocity -- needed by CrossProduct
#   * XVelocityToFlux / YVelocityToFlux -- needed by VelocityToFlux
#   * VelocityToFlux / FluxToVelocity   -- convenience wrappers
#   * InnerProduct(Flux, Flux)           -- needed by OutputEnergy
# Still not ported (no dependent yet): InnerProduct(Scalar, Scalar),
# FineGridInnerProduct(Scalar, Scalar), VorticityInnerProduct,
# FineGridVorticityInnerProduct, the Laplacian(Array2, dx, BC, Array2)
# low-level form (used only inside the already-ported EllipticSolver2d,
# which reimplements it directly).
#
# NOTE(port) -- vectorization judgment call for the velocity<->flux
# conversions: the C++ FluxToXVelocity/FluxToYVelocity/XVelocityToFlux/
# YVelocityToFlux routines compute one large rectangular INTERIOR sweep on the
# finest grid (a uniform two-point stencil), then patch a ring of BORDER /
# INTERFACE / CORNER points on each coarser grid using non-rectangular index
# ranges and cross-grid interpolation via Grid.c2f()/f2c(). Only the finest-
# grid interior sweep (and the trivially-rectangular coarsest-grid zero-BC
# boundary blocks) is vectorized with numpy slicing here; the multigrid
# border/interface/corner loops are kept as explicit Python loops that mirror
# the C++ index-by-index, because (a) their strides and coarse<->fine coupling
# make a slice rewrite error-prone and not obviously equivalent, and (b) that
# code path only executes for Ngrid > 1. For the common Ngrid == 1 case the
# entire computation is the vectorized interior/boundary blocks. Flagged per
# the "flag judgment calls" instruction rather than silently choosing.

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


# ---------------------------------------------------------------------------
# Laplacian
# ---------------------------------------------------------------------------

def Laplacian(f: Scalar, g: Optional[Scalar] = None) -> Optional[Scalar]:
    """Compute the Laplacian of Scalar f (assumes boundary values of f are
    zero), as `Laplacian = -Curl(Curl(f))`.

    NOTE(port): collapses the two Scalar C++ overloads
        void   Laplacian(const Scalar& f, Scalar& g);
        Scalar Laplacian(const Scalar& f);
    into one function dispatching on whether `g` is supplied (in-place form
    returns None; convenience form allocates, fills, and returns a Scalar) --
    the same overload-collapse convention used by Curl above. The third C++
    overload Laplacian(const Array2&, double, const BC&, Array2&) is not
    ported here (see the module scope note).
    """
    if g is None:
        out = Scalar(f.getGrid())
        Laplacian(f, out)
        return out

    assert f.Nx() == g.Nx()
    assert f.Ny() == g.Ny()
    assert f.Ngrid() == g.Ngrid()

    # Laplacian = - Curl( Curl( ) )
    q = Curl(f)         # Flux
    Curl(q, g)          # writes g in place
    g *= -1
    return None


# ---------------------------------------------------------------------------
# Velocity <-> flux conversions (see module vectorization note)
# ---------------------------------------------------------------------------

def _flux_x_views(q: Flux) -> list:
    """Per-level views of the X-fluxes of q, each shaped (Nx+1, Ny) so that
    the C++ accessor q(lev, X, i, j) maps to view[lev][i, j]."""
    nx, ny = q.Nx(), q.Ny()
    return [
        q._data[lev, q.begin(Direction.X):q.end(Direction.X)].reshape(nx + 1, ny)
        for lev in range(q.Ngrid())
    ]


def _flux_y_views(q: Flux) -> list:
    """Per-level views of the Y-fluxes of q, each shaped (Nx, Ny+1) so that
    q(lev, Y, i, j) maps to view[lev][i, j]."""
    nx, ny = q.Nx(), q.Ny()
    return [
        q._data[lev, q.begin(Direction.Y):q.end(Direction.Y)].reshape(nx, ny + 1)
        for lev in range(q.Ngrid())
    ]


def _scalar_views(u: Scalar) -> list:
    """Per-level views of Scalar u, each shaped (Nx-1, Ny-1) so that the C++
    accessor u(lev, i, j) maps to view[lev][i-1, j-1]."""
    return [u._data[lev] for lev in range(u.Ngrid())]


def FluxToXVelocity(q: Flux, u: Scalar) -> None:
    """Convert x-fluxes through edges to velocities at vertices (in place).

    NOTE(port): `u` is an out-parameter. See the module vectorization note:
    the finest-grid interior sweep (A) is vectorized; the coarse-grid
    border/interface/corner loops (B-F) are kept explicit.
    """
    assert q.Nx() == u.Nx()
    assert q.Ny() == u.Ny()
    assert q.Ngrid() == u.Ngrid()
    nx = q.Nx()
    ny = q.Ny()
    nx2 = q.NxExt()
    ny2 = q.NyExt()
    oneOver2Delta = 1.0 / (2 * q.Dx())
    X = _flux_x_views(q)
    U = _scalar_views(u)

    # Compute interior points (A) -- finest grid only
    # u(0,i,j) = ( q(0,X,i,j) + q(0,X,i,j-1) ) * oneOver2Delta, i,j in [1,nx-1]
    U[0][:, :] = (X[0][1:nx, 1:ny] + X[0][1:nx, 0:ny - 1]) * oneOver2Delta

    # Compute border points for each coarse grid (B-F)
    #
    # NOTE(port): the C++ per-element (i,j) loops labeled B-E are vectorized
    # here with numpy slicing (F, the four corners, stays scalar -- it is
    # O(1) work). Regions B, C, D all apply the identical two-point average
    # (X(i,j) + X(i,j-1)) * 0.5 * bydx, so a single full-level averaged array
    # `avgU` is formed once and copied into each region's slice; region E is
    # the cross-grid interface interpolation, expressed with fancy indexing.
    # This reproduces exactly the same writes (same cells, same values) as
    # the loop version, verified numerically against it for Ngrid > 1.
    for lev in range(1, q.Ngrid()):
        bydx = 1.0 / q.Dx(lev)
        Xl = X[lev]
        Ul = U[lev]
        # avgU[i-1, j-1] = (X(i,j) + X(i,j-1)) * 0.5 * bydx, for i,j in [1, n-1]
        avgU = (Xl[1:nx, 1:ny] + Xl[1:nx, 0:ny - 1]) * 0.5 * bydx

        # left and right borders (excluding interface) (B)
        Ul[0:nx2 - 1, :] = avgU[0:nx2 - 1, :]
        Ul[nx // 2 + nx2:nx - 1, :] = avgU[nx // 2 + nx2:nx - 1, :]
        # top and bottom borders (excluding interfaces) (C)
        Ul[nx2 - 1:nx // 2 + nx2, 0:ny2 - 1] = avgU[nx2 - 1:nx // 2 + nx2, 0:ny2 - 1]
        Ul[nx2 - 1:nx // 2 + nx2, ny // 2 + ny2:ny - 1] = avgU[nx2 - 1:nx // 2 + nx2, ny // 2 + ny2:ny - 1]
        # left and right interfaces, excluding corners (D)
        Ul[nx2 - 1, ny2:ny // 2 + ny2 - 1] = avgU[nx2 - 1, ny2:ny // 2 + ny2 - 1]
        Ul[nx // 2 + nx2 - 1, ny2:ny // 2 + ny2 - 1] = avgU[nx // 2 + nx2 - 1, ny2:ny // 2 + ny2 - 1]

        # top and bottom interfaces, excluding corners (E)
        Xf = X[lev - 1]
        i_arr = np.arange(nx2 + 1, nx // 2 + nx2)
        if i_arr.size:
            ii = (i_arr - nx2) * 2  # c2f(i, ny2) -> (ii, 0); c2f(i, ny/2+ny2) -> (ii, ny)
            # bottom interface (j = ny2, jj = 0)
            Ul[i_arr - 1, ny2 - 1] = (
                Xl[i_arr, ny2 - 1] * 2.0 / 3 + Xf[ii, 0] * 1.0 / 3
                + (Xf[ii - 1, 0] + Xf[ii + 1, 0]) * 1.0 / 6
            ) * bydx
            # top interface (j = ny/2+ny2, jj = ny)
            Ul[i_arr - 1, ny // 2 + ny2 - 1] = (
                Xl[i_arr, ny // 2 + ny2] * 2.0 / 3 + Xf[ii, ny - 1] * 1.0 / 3
                + (Xf[ii - 1, ny - 1] + Xf[ii + 1, ny - 1]) * 1.0 / 6
            ) * bydx
        # corners (F)
        # lower left
        i = nx2
        j = ny2
        U[lev][i - 1, j - 1] = (
            X[lev][i, j - 1] * 8.0 / 15 + X[lev][i, j] * 6.0 / 15 + X[lev - 1][1, 0] * 2.0 / 15
        ) * bydx
        # lower right
        i = nx // 2 + nx2
        U[lev][i - 1, j - 1] = (
            X[lev][i, j - 1] * 8.0 / 15 + X[lev][i, j] * 6.0 / 15 + X[lev - 1][nx - 1, 0] * 2.0 / 15
        ) * bydx
        # upper left
        i = nx2
        j = ny // 2 + ny2
        U[lev][i - 1, j - 1] = (
            X[lev][i, j] * 8.0 / 15 + X[lev][i, j - 1] * 6.0 / 15 + X[lev - 1][1, ny - 1] * 2.0 / 15
        ) * bydx
        # upper right
        i = nx // 2 + nx2
        U[lev][i - 1, j - 1] = (
            X[lev][i, j] * 8.0 / 15 + X[lev][i, j - 1] * 6.0 / 15
            + X[lev - 1][nx - 1, ny - 1] * 2.0 / 15
        ) * bydx


def FluxToYVelocity(q: Flux, v: Scalar) -> None:
    """Convert y-fluxes through edges to velocities at vertices (in place).

    NOTE(port): `v` is an out-parameter; see FluxToXVelocity / the module
    vectorization note.
    """
    assert q.Nx() == v.Nx()
    assert q.Ny() == v.Ny()
    assert q.Ngrid() == v.Ngrid()
    nx = q.Nx()
    ny = q.Ny()
    nx2 = q.NxExt()
    ny2 = q.NyExt()
    oneOver2Delta = 1.0 / (2 * q.Dx())
    Y = _flux_y_views(q)
    V = _scalar_views(v)

    # Compute interior points (A) -- finest grid only
    # v(0,i,j) = ( q(0,Y,i-1,j) + q(0,Y,i,j) ) * oneOver2Delta, i,j in [1,nx-1]
    V[0][:, :] = (Y[0][0:nx - 1, 1:ny] + Y[0][1:nx, 1:ny]) * oneOver2Delta

    # NOTE(port): C++ per-element loops B-E vectorized with numpy slicing
    # (F corners kept scalar); see the analogous note in FluxToXVelocity.
    # Regions B, C, D share the average (Y(i,j) + Y(i-1,j)) * 0.5 * bydx,
    # formed once as `avgV`; E is the cross-grid interface interpolation.
    for lev in range(1, q.Ngrid()):
        bydx = 1.0 / q.Dx(lev)
        Yl = Y[lev]
        Vl = V[lev]
        # avgV[i-1, j-1] = (Y(i,j) + Y(i-1,j)) * 0.5 * bydx, for i,j in [1, n-1]
        avgV = (Yl[1:nx, 1:ny] + Yl[0:nx - 1, 1:ny]) * 0.5 * bydx

        # top and bottom borders (excluding interface) (B)
        Vl[:, 0:ny2 - 1] = avgV[:, 0:ny2 - 1]
        Vl[:, ny // 2 + ny2:ny - 1] = avgV[:, ny // 2 + ny2:ny - 1]
        # left and right borders (excluding interfaces) (C)
        Vl[0:nx2 - 1, ny2 - 1:ny // 2 + ny2] = avgV[0:nx2 - 1, ny2 - 1:ny // 2 + ny2]
        Vl[nx // 2 + nx2:nx - 1, ny2 - 1:ny // 2 + ny2] = avgV[nx // 2 + nx2:nx - 1, ny2 - 1:ny // 2 + ny2]
        # top and bottom interfaces, excluding corners (D)
        Vl[nx2:nx // 2 + nx2 - 1, ny2 - 1] = avgV[nx2:nx // 2 + nx2 - 1, ny2 - 1]
        Vl[nx2:nx // 2 + nx2 - 1, ny // 2 + ny2 - 1] = avgV[nx2:nx // 2 + nx2 - 1, ny // 2 + ny2 - 1]

        # left and right interfaces, excluding corners (E)
        Yf = Y[lev - 1]
        j_arr = np.arange(ny2 + 1, ny // 2 + ny2)
        if j_arr.size:
            jj = (j_arr - ny2) * 2  # c2f(nx2, j) -> (0, jj); c2f(nx/2+nx2, j) -> (nx, jj)
            # left interface (i = nx2, ii = 0)
            Vl[nx2 - 1, j_arr - 1] = (
                Yl[nx2 - 1, j_arr] * 2.0 / 3 + Yf[0, jj] * 1.0 / 3
                + (Yf[0, jj - 1] + Yf[0, jj + 1]) * 1.0 / 6
            ) * bydx
            # right interface (i = nx/2+nx2, ii = nx)
            Vl[nx // 2 + nx2 - 1, j_arr - 1] = (
                Yl[nx // 2 + nx2, j_arr] * 2.0 / 3 + Yf[nx - 1, jj] * 1.0 / 3
                + (Yf[nx - 1, jj - 1] + Yf[nx - 1, jj + 1]) * 1.0 / 6
            ) * bydx
        # corners (F)
        j = ny2
        i = nx2
        V[lev][i - 1, j - 1] = (
            Y[lev][i - 1, j] * 8.0 / 15 + Y[lev][i, j] * 6.0 / 15 + Y[lev - 1][0, 1] * 2.0 / 15
        ) * bydx
        j = ny // 2 + ny2
        V[lev][i - 1, j - 1] = (
            Y[lev][i - 1, j] * 8.0 / 15 + Y[lev][i, j] * 6.0 / 15 + Y[lev - 1][0, ny - 1] * 2.0 / 15
        ) * bydx
        j = ny2
        i = nx // 2 + nx2
        V[lev][i - 1, j - 1] = (
            Y[lev][i, j] * 8.0 / 15 + Y[lev][i - 1, j] * 6.0 / 15 + Y[lev - 1][nx - 1, 1] * 2.0 / 15
        ) * bydx
        j = ny // 2 + ny2
        V[lev][i - 1, j - 1] = (
            Y[lev][i, j] * 8.0 / 15 + Y[lev][i - 1, j] * 6.0 / 15
            + Y[lev - 1][nx - 1, ny - 1] * 2.0 / 15
        ) * bydx


def XVelocityToFlux(u: Scalar, q: Flux) -> None:
    """Convert u-velocities at vertices to x-fluxes through edges (in place).
    Does not touch the y-component of q.

    NOTE(port): `q` is an out-parameter (X-component only). All regions
    (finest-grid interior D, coarse-grid B/D/G, and boundaries A/C) are
    vectorized with numpy slicing/fancy-indexing; the C++ per-element (i,j)
    loops are reproduced exactly (verified numerically for Ngrid > 1).
    """
    assert u.Nx() == q.Nx()
    assert u.Ny() == q.Ny()
    assert u.Ngrid() == q.Ngrid()
    nx = u.Nx()
    ny = u.Ny()
    nx2 = u.NxExt()
    ny2 = u.NyExt()
    g = q.getGrid()
    X = _flux_x_views(q)
    U = _scalar_views(u)

    for lev in range(u.Ngrid()):
        dx = g.Dx(lev)
        # Interior points
        if lev == 0:
            # interior points on finest grid, minus top and bottom rows (D)
            # q(0,X,i,j) = ( u(0,i,j) + u(0,i,j+1) ) * 0.5 * dx, i in [1,nx-1], j in [1,ny-2]
            X[0][1:nx, 1:ny - 1] = (U[0][0:nx - 1, 0:ny - 2] + U[0][0:nx - 1, 1:ny - 1]) * 0.5 * dx
        else:  # not the finest grid
            # NOTE(port): C++ per-element (i,j) loops B, D, G vectorized with
            # numpy slicing (B and D share the average (u(i,j)+u(i,j+1))*0.5*dx,
            # formed once as `avgX`; G copies from the finer grid). Verified
            # numerically against the loop version for Ngrid > 1.
            Xl = X[lev]
            Ul = U[lev]
            # avgX[i, j] = (u(i,j) + u(i,j+1)) * 0.5 * dx for i in [1,nx-1], j in [1,ny-2]
            # (indexing Ul[i-1, j-1] for u(i,j)); built lazily per region below.
            # top and bottom portions, excluding outer interface (B)
            Xl[1:nx, 1:ny2] = (Ul[0:nx - 1, 0:ny2 - 1] + Ul[0:nx - 1, 1:ny2]) * 0.5 * dx
            Xl[1:nx, ny // 2 + ny2:ny - 1] = (
                Ul[0:nx - 1, ny // 2 + ny2 - 1:ny - 2] + Ul[0:nx - 1, ny // 2 + ny2:ny - 1]
            ) * 0.5 * dx
            # left and right portions of coarse grid (D)
            Xl[1:nx2 + 1, ny2:ny // 2 + ny2] = (
                Ul[0:nx2, ny2 - 1:ny // 2 + ny2 - 1] + Ul[0:nx2, ny2:ny // 2 + ny2]
            ) * 0.5 * dx
            Xl[nx // 2 + nx2:nx, ny2:ny // 2 + ny2] = (
                Ul[nx // 2 + nx2 - 1:nx - 1, ny2 - 1:ny // 2 + ny2 - 1]
                + Ul[nx // 2 + nx2 - 1:nx - 1, ny2:ny // 2 + ny2]
            ) * 0.5 * dx
            # get interior portion of coarse grid from fine grid (G)
            Xf = X[lev - 1]
            i_arr = np.arange(nx2 + 1, nx // 2 + nx2)
            j_arr = np.arange(ny2, ny // 2 + ny2)
            if i_arr.size and j_arr.size:
                ii = (i_arr - nx2) * 2  # c2f(i,j) -> (ii, jj)
                jj = (j_arr - ny2) * 2
                Xl[np.ix_(i_arr, j_arr)] = Xf[np.ix_(ii, jj)] + Xf[np.ix_(ii, jj + 1)]
        # Boundary points
        # left and right boundaries of coarsest grid are zero (A)
        if lev == g.Ngrid() - 1:
            X[lev][0, 0:ny] = 0
            X[lev][nx, 0:ny] = 0
        # left and right boundaries of finer grids take values from coarser grid (A)
        else:
            # NOTE(port): C++ `for j=0; j<ny-1; j+=2` loop vectorized; f2c(0,j)
            # gives ii=nx2 (unused) and jj=j//2+ny2.
            Uf = U[lev + 1]
            j_arr = np.arange(0, ny - 1, 2)
            jj = j_arr // 2 + ny2
            a = Uf[nx2 - 1, jj - 1]
            b = Uf[nx2 - 1, jj]
            c = Uf[nx // 2 + nx2 - 1, jj - 1]
            d = Uf[nx // 2 + nx2 - 1, jj]
            X[lev][0, j_arr] = (0.75 * a + 0.25 * b) * dx
            X[lev][nx, j_arr] = (0.75 * c + 0.25 * d) * dx
            X[lev][0, j_arr + 1] = (0.25 * a + 0.75 * b) * dx
            X[lev][nx, j_arr + 1] = (0.25 * c + 0.75 * d) * dx
        # outer interface (top/bottom), get values from coarser grid (or zero, for coarsest) (C)
        if lev == u.Ngrid() - 1:
            # on coarsest grid: zero bcs
            # q(lev,X,i,0) = u(lev,i,1)*0.5*dx ; q(lev,X,i,ny-1) = u(lev,i,ny-1)*0.5*dx, i in [1,nx-1]
            X[lev][1:nx, 0] = U[lev][0:nx - 1, 0] * 0.5 * dx
            X[lev][1:nx, ny - 1] = U[lev][0:nx - 1, ny - 2] * 0.5 * dx
        else:
            # on intermediate grid: get bcs from coarser grid.
            # NOTE(port): the C++ nests these two i-loops inside `for i=1;i<nx`
            # but re-declares `i` in the inner loops, so the outer i is unused;
            # the inner loops run once each (per lev). Both are vectorized here.
            Uf = U[lev + 1]
            Ul = U[lev]
            # points that correspond to coarse points (even i)
            ie = np.arange(2, nx, 2)
            iie = ie // 2 + nx2  # f2c(i,0) -> ii = i//2 + nx2
            X[lev][ie, 0] = (Ul[ie - 1, 0] + Uf[iie - 1, ny2 - 1]) * 0.5 * dx
            X[lev][ie, ny - 1] = (Ul[ie - 1, ny - 2] + Uf[iie - 1, ny // 2 + ny2 - 1]) * 0.5 * dx
            # points that do not correspond to coarse points (odd i)
            io = np.arange(1, nx, 2)
            iio = io // 2 + nx2
            X[lev][io, 0] = (
                0.5 * Ul[io - 1, 0] + 0.25 * Uf[iio - 1, ny2 - 1] + 0.25 * Uf[iio, ny2 - 1]
            ) * dx
            X[lev][io, ny - 1] = (
                0.5 * Ul[io - 1, ny - 2]
                + 0.25 * Uf[iio - 1, ny // 2 + ny2 - 1] + 0.25 * Uf[iio, ny // 2 + ny2 - 1]
            ) * dx


def YVelocityToFlux(v: Scalar, q: Flux) -> None:
    """Convert v-velocities at vertices to y-fluxes through edges (in place).
    Does not touch the x-component of q.

    NOTE(port): `q` is an out-parameter (Y-component only). All regions are
    vectorized with numpy slicing/fancy-indexing (see XVelocityToFlux);
    verified numerically against the loop version for Ngrid > 1.
    """
    assert v.Nx() == q.Nx()
    assert v.Ny() == q.Ny()
    assert v.Ngrid() == q.Ngrid()
    nx = v.Nx()
    ny = v.Ny()
    nx2 = v.NxExt()
    ny2 = v.NyExt()
    g = v.getGrid()
    Y = _flux_y_views(q)
    V = _scalar_views(v)

    for lev in range(g.Ngrid()):
        dx = g.Dx(lev)
        # Interior points
        if lev == 0:
            # interior points on finest grid, minus left and right rows (D)
            # q(0,Y,i,j) = ( v(0,i,j) + v(0,i+1,j) ) * 0.5 * dx, j in [1,ny-1], i in [1,nx-2]
            Y[0][1:nx - 1, 1:ny] = (V[0][0:nx - 2, 0:ny - 1] + V[0][1:nx - 1, 0:ny - 1]) * 0.5 * dx
        else:  # not the finest grid
            # NOTE(port): C++ per-element (i,j) loops B, D, G vectorized with
            # numpy slicing (B and D share (v(i,j)+v(i+1,j))*0.5*dx; G copies
            # from the finer grid). Verified numerically for Ngrid > 1.
            Yl = Y[lev]
            Vl = V[lev]
            # left and right portions, excluding outer interface (B)
            Yl[1:nx2, 1:ny] = (Vl[0:nx2 - 1, 0:ny - 1] + Vl[1:nx2, 0:ny - 1]) * 0.5 * dx
            Yl[nx // 2 + nx2:nx - 1, 1:ny] = (
                Vl[nx // 2 + nx2 - 1:nx - 2, 0:ny - 1] + Vl[nx // 2 + nx2:nx - 1, 0:ny - 1]
            ) * 0.5 * dx
            # top and bottom portions of coarse grid (D)
            Yl[nx2:nx // 2 + nx2, 1:ny2 + 1] = (
                Vl[nx2 - 1:nx // 2 + nx2 - 1, 0:ny2] + Vl[nx2:nx // 2 + nx2, 0:ny2]
            ) * 0.5 * dx
            Yl[nx2:nx // 2 + nx2, ny // 2 + ny2:ny] = (
                Vl[nx2 - 1:nx // 2 + nx2 - 1, ny // 2 + ny2 - 1:ny - 1]
                + Vl[nx2:nx // 2 + nx2, ny // 2 + ny2 - 1:ny - 1]
            ) * 0.5 * dx
            # get interior portion of coarse grid from fine grid (G)
            Yf = Y[lev - 1]
            i_arr = np.arange(nx2, nx // 2 + nx2)
            j_arr = np.arange(ny2 + 1, ny // 2 + ny2)
            if i_arr.size and j_arr.size:
                ii = (i_arr - nx2) * 2  # c2f(i,j) -> (ii, jj)
                jj = (j_arr - ny2) * 2
                Yl[np.ix_(i_arr, j_arr)] = Yf[np.ix_(ii, jj)] + Yf[np.ix_(ii + 1, jj)]
        # Boundary points
        # top and bottom boundaries of coarsest grid are zero (A)
        if lev == g.Ngrid() - 1:
            Y[lev][0:nx, 0] = 0
            Y[lev][0:nx, ny] = 0
        # top and bottom boundaries of finer grids take values from coarser grid (A)
        else:
            # NOTE(port): C++ `for i=0; i<nx-1; i+=2` loop vectorized; f2c(i,0)
            # gives ii=i//2+nx2 and jj=ny2 (unused).
            Vf = V[lev + 1]
            i_arr = np.arange(0, nx - 1, 2)
            ii = i_arr // 2 + nx2
            a = Vf[ii - 1, ny2 - 1]
            b = Vf[ii, ny2 - 1]
            c = Vf[ii - 1, ny // 2 + ny2 - 1]
            d = Vf[ii, ny // 2 + ny2 - 1]
            Y[lev][i_arr, 0] = (0.75 * a + 0.25 * b) * dx
            Y[lev][i_arr, ny] = (0.75 * c + 0.25 * d) * dx
            Y[lev][i_arr + 1, 0] = (0.25 * a + 0.75 * b) * dx
            Y[lev][i_arr + 1, ny] = (0.25 * c + 0.75 * d) * dx
        # outer interface (left/right), get values from coarser grid (or zero, for coarsest) (C)
        if lev == g.Ngrid() - 1:
            # on coarsest grid: zero bcs
            Y[lev][0, 1:ny] = V[lev][0, 0:ny - 1] * 0.5 * dx
            Y[lev][nx - 1, 1:ny] = V[lev][nx - 2, 0:ny - 1] * 0.5 * dx
        else:
            # on intermediate grid: get bcs from coarser grid.
            # NOTE(port): as in XVelocityToFlux, the C++ nests these inside an
            # outer `for j=1;j<ny` with a re-declared inner `j`; the outer j is
            # unused, so the inner loops run once. Both are vectorized here.
            Vf = V[lev + 1]
            Vl = V[lev]
            # points that correspond to coarse points (even j)
            je = np.arange(2, ny, 2)
            jje = je // 2 + ny2  # f2c(0,j) -> jj = j//2 + ny2
            Y[lev][0, je] = (Vl[0, je - 1] + Vf[nx2 - 1, jje - 1]) * 0.5 * dx
            Y[lev][nx - 1, je] = (Vl[nx - 2, je - 1] + Vf[nx // 2 + nx2 - 1, jje - 1]) * 0.5 * dx
            # points that do not correspond to coarse points (odd j)
            jo = np.arange(1, ny, 2)
            jjo = jo // 2 + ny2
            Y[lev][0, jo] = (
                0.5 * Vl[0, jo - 1] + 0.25 * Vf[nx2 - 1, jjo - 1] + 0.25 * Vf[nx2 - 1, jjo]
            ) * dx
            Y[lev][nx - 1, jo] = (
                0.5 * Vl[nx - 2, jo - 1]
                + 0.25 * Vf[nx // 2 + nx2 - 1, jjo - 1] + 0.25 * Vf[nx // 2 + nx2 - 1, jjo]
            ) * dx


def VelocityToFlux(u: Scalar, v: Scalar, q: Flux) -> None:
    """Convert u- and v-velocities at vertices to fluxes through edges."""
    XVelocityToFlux(u, q)
    YVelocityToFlux(v, q)


def FluxToVelocity(q: Flux, u: Scalar, v: Scalar) -> None:
    """Convert fluxes through edges to u- and v-velocities at vertices."""
    FluxToXVelocity(q, u)
    FluxToYVelocity(q, v)


# ---------------------------------------------------------------------------
# Cross products
# ---------------------------------------------------------------------------

def CrossProduct(q: Flux, arg: Union[Flux, Scalar]) -> Union[Flux, Scalar]:
    """Cross product of a Flux with a Scalar (returns a Flux) or with another
    Flux (returns a Scalar).

    NOTE(port): collapses the two C++ overloads
        Flux   CrossProduct(const Flux& q,  const Scalar& f);   // q x f = (f v, -f u)
        Scalar CrossProduct(const Flux& q1, const Flux& q2);    // q1 x q2 = u1 v2 - u2 v1
    into one function dispatching on the type of the second argument.
    """
    if isinstance(arg, Scalar):
        f = arg
        assert q.Nx() == f.Nx()
        assert q.Ny() == f.Ny()
        assert q.Ngrid() == f.Ngrid()

        u = Scalar(f.getGrid())
        v = Scalar(f.getGrid())

        FluxToXVelocity(q, u)
        u *= f
        u *= -1

        FluxToYVelocity(q, v)
        v *= f

        cross = Flux(q.getGrid())
        VelocityToFlux(v, u, cross)  # cross = ( f v, -f u )
        return cross

    elif isinstance(arg, Flux):
        q1 = q
        q2 = arg
        assert q1.Nx() == q2.Nx()
        assert q1.Ny() == q2.Ny()
        assert q1.Ngrid() == q2.Ngrid()

        grid = q1.getGrid()
        u = Scalar(grid)
        v = Scalar(grid)

        FluxToXVelocity(q1, u)
        FluxToYVelocity(q2, v)

        f = u * v  # f is now u1 * v2

        FluxToXVelocity(q2, u)
        FluxToYVelocity(q1, v)

        f -= u * v  # f is now u1 * v2 - u2 * v1

        f.coarsify()  # fill in overlapping grid regions
        return f

    else:
        raise TypeError(f"CrossProduct: unsupported second argument type {type(arg)!r}")


# ---------------------------------------------------------------------------
# Inner products (Flux form only ported so far -- see module scope note)
# ---------------------------------------------------------------------------

def _fine_grid_inner_product_flux(p: Flux, q: Flux) -> float:
    """Inner product of Flux p and Flux q, restricted to the finest grid
    (interior points only). Note this is not multiplied by dx * dx, since
    the Fluxes are already multiplied by dx (i.e. the inner product is
    really over *velocities*).

    NOTE(port): the C++ nested (i,j) loops over the X- and Y-flux interior
    ranges are vectorized below with numpy slicing/sum, instead of a manual
    per-element Python loop.
    """
    nx = p.Nx()
    ny = p.Ny()
    Xp = _flux_x_views(p)[0]
    Xq = _flux_x_views(q)[0]
    Yp = _flux_y_views(p)[0]
    Yq = _flux_y_views(q)[0]

    ip = float(np.sum(Xp[1:nx, :] * Xq[1:nx, :]))
    ip += float(np.sum(Yp[:, 1:ny] * Yq[:, 1:ny]))
    return ip


def InnerProduct(p: Flux, q: Flux) -> float:
    """Inner product of Flux p and Flux q, over all grid levels.

    NOTE(port): only the Flux overload of C++'s
        double InnerProduct(const Scalar& f, const Scalar& g);
        double InnerProduct(const Flux& p, const Flux& q);
    is ported here (needed by OutputEnergy); see the module scope note at
    the top of this file. If a future dependent needs the Scalar overload,
    add it here following the same pattern (and re-introduce the type
    dispatch collapsing both, matching Curl/Laplacian/CrossProduct above).

    NOTE(port): the C++ nested (lev, i, j) loops over the coarser-grid
    interface/border/coarse-point ranges are vectorized below with numpy
    slicing/sum on each level's X/Y views; the outer level loop is kept as
    a plain Python loop (Ngrid is small, and this path only executes for
    Ngrid > 1), matching the same convention used in FluxToXVelocity/
    FluxToYVelocity/XVelocityToFlux/YVelocityToFlux above.
    """
    assert p.Ngrid() == q.Ngrid()
    assert p.Nx() == q.Nx()
    assert p.Ny() == q.Ny()
    nx = p.Nx()
    ny = p.Ny()
    nx2 = p.NxExt()
    ny2 = p.NyExt()

    ip = _fine_grid_inner_product_flux(p, q)

    Xp_all = _flux_x_views(p)
    Xq_all = _flux_x_views(q)
    Yp_all = _flux_y_views(p)
    Yq_all = _flux_y_views(q)

    for lev in range(1, p.Ngrid()):
        Xp, Xq = Xp_all[lev], Xq_all[lev]
        Yp, Yq = Yp_all[lev], Yq_all[lev]

        # X-fluxes, coarser grids
        # left and right interfaces (edges)
        ip += float(np.sum(Xp[nx2, ny2:ny // 2 + ny2] * Xq[nx2, ny2:ny // 2 + ny2])) * 0.75
        ip += float(
            np.sum(Xp[nx // 2 + nx2, ny2:ny // 2 + ny2] * Xq[nx // 2 + nx2, ny2:ny // 2 + ny2])
        ) * 0.75
        # left and right coarse points
        ip += float(np.sum(Xp[1:nx2, :] * Xq[1:nx2, :]))
        ip += float(np.sum(Xp[nx // 2 + nx2 + 1:nx, :] * Xq[nx // 2 + nx2 + 1:nx, :]))
        # top and bottom coarse points
        ip += float(np.sum(Xp[nx2:nx // 2 + nx2 + 1, 0:ny2] * Xq[nx2:nx // 2 + nx2 + 1, 0:ny2]))
        ip += float(
            np.sum(Xp[nx2:nx // 2 + nx2 + 1, ny // 2 + ny2:ny] * Xq[nx2:nx // 2 + nx2 + 1, ny // 2 + ny2:ny])
        )

        # Y-fluxes, coarser grids
        # left and right interfaces (edges)
        ip += float(np.sum(Yp[nx2:nx // 2 + nx2, ny2] * Yq[nx2:nx // 2 + nx2, ny2])) * 0.75
        ip += float(
            np.sum(Yp[nx2:nx // 2 + nx2, ny // 2 + ny2] * Yq[nx2:nx // 2 + nx2, ny // 2 + ny2])
        ) * 0.75
        # left and right coarse points
        ip += float(np.sum(Yp[:, 1:ny2] * Yq[:, 1:ny2]))
        ip += float(np.sum(Yp[:, ny // 2 + ny2 + 1:ny] * Yq[:, ny // 2 + ny2 + 1:ny]))
        # top and bottom coarse points
        ip += float(np.sum(Yp[0:nx2, ny2:ny // 2 + ny2 + 1] * Yq[0:nx2, ny2:ny // 2 + ny2 + 1]))
        ip += float(
            np.sum(Yp[nx // 2 + nx2:nx, ny2:ny // 2 + ny2 + 1] * Yq[nx // 2 + nx2:nx, ny2:ny // 2 + ny2 + 1])
        )

    return ip
