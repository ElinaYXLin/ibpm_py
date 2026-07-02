# regularizer.py
#
# Python port of src/Regularizer.h / src/Regularizer.cc
#
# Define regularization operations between grid and boundary data.
#
# Regularization (smearing) lifts values defined on the boundary to values
# defined on a grid, and interpolation takes data defined on a grid and
# interpolates it to the boundary.
#
# Specifically, if u(x) are values of u defined on a grid, and if u(xi)
# denotes values on a boundary, then regularization defines
#
#     u(x) = \int_\Omega u(xi) delta(x-xi) dxi
#
# where a discrete approximation of the delta function is used, and
# interpolation defines
#
#     u(xi) = \int_\Omega u(x) delta(x-xi) dx
#
# These integrals are discretized using a regularized version of the delta
# function, with finite support, as in (14) of Taira & Colonius (J Comput
# Phys, 2007).

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .boundary_vector import BoundaryVector
from .direction import Direction
from .flux import Flux

if TYPE_CHECKING:
    from .geometry import Geometry
    from .grid import Grid

FloatArray = np.ndarray

# number of cells over which delta function has support
deltaSupportRadius: float = 1.5


def deltaFunction(r: FloatArray) -> FloatArray:
    """Return the value of the regularized delta function phi(r), where
    delta(x) ~= phi(x/h) / h, and h is the grid spacing.

    From Roma, Peskin, and Berger, JCP 1999, eq (22). Note that the input r
    is given in numbers of cells (e.g. normalized by the cell width).

    NOTE(port): C++ `deltaFunction` takes a scalar `double r` and is called
    once per (boundary point, grid cell) pair inside the triple-nested loop
    in `update()`. Here it is vectorized to accept a numpy array of r values
    (broadcast over all pairs at once), using `np.where` in place of the
    scalar if/else branches -- this mirrors the same
    scalar-branch-to-`np.where` translation used for `Scalar`/`Flux`
    elementwise operators elsewhere in this port.
    """
    r = np.asarray(r, dtype=np.float64)
    near = (1 + np.sqrt(np.maximum(1 - 3 * r * r, 0.0))) / 3
    far = (5 - 3 * r - np.sqrt(np.maximum(1 - 3 * (1 - r) * (1 - r), 0.0))) / 6
    val = np.where(r <= 0.5, near, far)
    return np.where(r > deltaSupportRadius, 0.0, val)


class Regularizer:
    """Define regularization operations between grid and boundary data."""

    def __init__(self, grid: "Grid", geometry: "Geometry") -> None:
        self._grid = grid
        self._geometry = geometry
        # NOTE(port): C++ stores the associations between boundary points
        # and nearby Flux values as `vector<Association>`, a list of
        # {boundaryIndex, fluxIndex, weight} structs, built and consumed
        # with explicit loops. Since the task calls for numpy vectorization
        # in place of manual loops, the association list is instead stored
        # as three parallel numpy arrays (`_boundaryIndex`, `_fluxIndex`,
        # `_weight`), populated in `update()` via vectorized numpy
        # computation and consumed in `toFlux`/`toBoundary` via
        # `np.add.at` scatter-adds instead of a Python for-loop over
        # `Association` objects. This changes the *storage/iteration*
        # mechanism only, not the algorithm: the same set of
        # (boundaryIndex, fluxIndex, weight) triples is computed and
        # accumulated, in the same input order per direction/body/cell.
        self._boundaryIndex: FloatArray = np.empty(0, dtype=np.int64)
        self._fluxIndex: FloatArray = np.empty(0, dtype=np.int64)
        self._weight: FloatArray = np.empty(0, dtype=np.float64)

    def update(self) -> None:
        """Update list of relationships between boundary points and cells,
        and the corresponding weights. Checks only the finest grid level,
        level=0."""
        h = self._grid.Dx()  # mesh spacing

        # Get the coordinates of the body
        bodyCoords = self._geometry.getPoints()
        numPoints = bodyCoords.getNumPoints()

        # NOTE(port): access to `bodyCoords._data` (rather than only the
        # public `BoundaryVector.__call__` interface) mirrors the existing
        # convention used elsewhere in this port (e.g. geometry.py's
        # `getPoints`/`getVelocities`), where sibling classes reach into
        # `_data` directly to vectorize bulk copies/computations.
        bx = bodyCoords._data[bodyCoords.begin(Direction.X):bodyCoords.end(Direction.X)]
        by = bodyCoords._data[bodyCoords.begin(Direction.Y):bodyCoords.end(Direction.Y)]

        boundaryIndex_parts = []
        fluxIndex_parts = []
        weight_parts = []

        # For each direction (x and y) -- matches the `for (dir = X; dir <= Y; ++dir)`
        # loop in the C++ implementation.
        for dir_ in (Direction.X, Direction.Y):
            nx = self._grid.Nx()
            ny = self._grid.Ny()
            if dir_ == Direction.X:
                # X fluxes: i in 0..nx, j in 0..ny-1
                i_vals = np.arange(nx + 1)
                j_vals = np.arange(ny)
                x0 = self._grid.getXEdge(0, 0)
                y0 = self._grid.getYCenter(0, 0)
                x_cells = x0 + i_vals * h
                y_cells = y0 + j_vals * h
                ny_stride = ny  # getIndex(X,i,j) = i*(ny+0)+j
            else:
                # Y fluxes: i in 0..nx-1, j in 0..ny
                i_vals = np.arange(nx)
                j_vals = np.arange(ny + 1)
                x0 = self._grid.getXCenter(0, 0)
                y0 = self._grid.getYEdge(0, 0)
                x_cells = x0 + i_vals * h
                y_cells = y0 + j_vals * h
                ny_stride = ny + 1  # getIndex(Y,i,j) = numXFluxes + i*(ny+1)+j

            # Build the (i,j) grid of cell coordinates and the corresponding
            # flat Flux index for each cell, matching Flux.getIndex(dir,i,j).
            X_cell, Y_cell = np.meshgrid(x_cells, y_cells, indexing="ij")
            I_idx, J_idx = np.meshgrid(i_vals, j_vals, indexing="ij")
            flux_index_grid = I_idx * ny_stride + J_idx
            if dir_ == Direction.Y:
                # Flux.begin(Y) == numXFluxes == nx*ny + ny (see Flux.resize)
                flux_index_grid = flux_index_grid + (nx * ny + ny)

            x_cells_flat = X_cell.ravel()
            y_cells_flat = Y_cell.ravel()
            flux_index_flat = flux_index_grid.ravel()

            # Pairwise distances (in units of cells) between every boundary
            # point and every cell in this direction -- this is the
            # vectorized equivalent of the innermost two loops
            # ("for each point on the boundary" / "for each cell") in the
            # C++ `update()`.
            dx = np.abs(x_cells_flat[None, :] - bx[:, None]) / h
            dy = np.abs(y_cells_flat[None, :] - by[:, None]) / h

            mask = (dx < deltaSupportRadius) & (dy < deltaSupportRadius)
            weight = deltaFunction(dx) * deltaFunction(dy)

            point_idx, cell_idx = np.nonzero(mask)
            if point_idx.size == 0:
                continue

            boundaryIndex_parts.append(bodyCoords.getIndex(dir_, 0) + point_idx)
            fluxIndex_parts.append(flux_index_flat[cell_idx])
            weight_parts.append(weight[point_idx, cell_idx])

        if boundaryIndex_parts:
            self._boundaryIndex = np.concatenate(boundaryIndex_parts).astype(np.int64)
            self._fluxIndex = np.concatenate(fluxIndex_parts).astype(np.int64)
            self._weight = np.concatenate(weight_parts).astype(np.float64)
        else:
            self._boundaryIndex = np.empty(0, dtype=np.int64)
            self._fluxIndex = np.empty(0, dtype=np.int64)
            self._weight = np.empty(0, dtype=np.float64)

    def toFlux(self, u: BoundaryVector) -> Flux:
        """Smear boundary data to grid.

        In particular, if u1 denotes the vectors along the boundary, and u2
        denotes the velocity vectors in the 2d domain, computes a discrete
        approximation to
            u2(x,y) = \\int u1(xi,eta) delta(x-xi) delta(y-eta) dxi deta
                ~= sum u1(xi,eta) delta(x-xi) delta(y-eta) * dx^2
        The flux returned is the corresponding flux through cell edges, or
        u2 * dx.
        """
        u1 = u
        # Allocate a new Flux field, initialized to zero
        u2 = Flux(self._grid)
        u2.assign(0)

        # For each association between cells and boundary points, add the
        # weight factor times the boundary value to the flux.
        #
        # NOTE(port): C++ loops over `_neighbors`, doing
        # `u2(0,a->fluxIndex) += a->weight * u1(a->boundaryIndex);` one
        # association at a time. `np.add.at` is the numpy scatter-add
        # equivalent for the (possibly repeated) `fluxIndex` targets --
        # plain fancy-index assignment (`u2._data[0, self._fluxIndex] +=
        # ...`) would silently drop all but one contribution when the same
        # flux index appears more than once (buffered, not accumulated),
        # which does not match the sequential C++ `+=` loop.
        #
        # NOTE(port) -- judgment call: the summation order of `np.add.at`
        # is not guaranteed to match the sequential per-association order
        # of the C++ loop. Since floating-point addition is not strictly
        # associative, results may differ from the C++ implementation at
        # the level of machine epsilon (this cannot be avoided while also
        # vectorizing the accumulation).
        np.add.at(u2._data[0], self._fluxIndex, self._weight * u1._data[self._boundaryIndex])

        # Multiply by grid spacing for correct dimension (vector -> Flux)
        u2 *= self._grid.Dx()

        # Return the new flux field
        return u2

    def toBoundary(self, u: Flux) -> BoundaryVector:
        """Interpolate grid data to boundary.

        In particular, if q denotes fluxes in the 2D domain, compute
        corresponding velocities u2 = q / dx, and define velocities u1 at
        boundary points by interpolation:
            u1(x,y) = \\int u2(xi,eta) delta(x-xi) delta(y-eta) dxi deta
                ~= sum u2(xi,eta) delta(x-xi) delta(y-eta) * dx^2
        """
        u2 = u
        # Allocate a new BoundaryVector, initialized to zero
        u1 = BoundaryVector(self._geometry.getNumPoints())
        u1.assign(0)

        # See the NOTE(port) in toFlux() regarding np.add.at vs a manual
        # per-association loop, and the summation-order judgment call.
        np.add.at(u1._data, self._boundaryIndex, self._weight * u2._data[0, self._fluxIndex])

        # Divide by grid spacing for correct dimension (Flux -> vector)
        u1 /= self._grid.Dx()

        # Return the new BoundaryVector
        return u1
