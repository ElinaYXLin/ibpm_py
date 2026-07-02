#!/usr/bin/env python3
"""Python-side raw-value dump matching cpp/dump_adjoint.cc, for
cross-validating py.vector_operations.FluxToYVelocity / YVelocityToFlux
against the C++ reference (verifies YVelocityToFlux is the adjoint of
FluxToYVelocity, per basis vector, across a 3-level multigrid).

Run from the repo root:  python3 -m py.tests.cross_validation.python.dump_adjoint
"""

from py.grid import Grid
from py.flux import Flux
from py.scalar import Scalar
from py.direction import Direction
from py.vector_operations import FluxToYVelocity, YVelocityToFlux
from .inner_products import InnerProductFlux, InnerProductScalar


def AdjYVelToFlux(q):
    grid = q.getGrid()
    x = Scalar(grid)
    x.assign(0.0)
    for lev in range(x.Ngrid()):
        for i in range(1, x.Nx()):
            for j in range(1, x.Ny()):
                e = Scalar(grid)
                e.assign(0.0)
                e._data[lev][i - 1, j - 1] = 1.0
                p = Flux(grid)
                p.assign(0.0)
                YVelocityToFlux(e, p)
                a = InnerProductFlux(q, p)
                normsq = InnerProductScalar(e, e)
                x._data[lev][i - 1, j - 1] = a / normsq
    return x


def dumpFluxY(name, q):
    print(name)
    for lev in range(q.Ngrid()):
        for i in range(0, q.Nx()):
            for j in range(0, q.Ny() + 1):
                print(repr(q(lev, Direction.Y, i, j)))


def main():
    nx, ny, ngrid = 8, 8, 3
    length = 0.1
    grid = Grid(nx, ny, ngrid, length, 0.0, 0.0)

    mag = Flux(grid)
    err2 = Flux(grid)
    mag.assign(0.0)
    err2.assign(0.0)

    for lev in range(ngrid):
        for ind in range(err2.begin(Direction.Y), err2.end(Direction.Y)):
            e = Flux(grid)
            e.assign(0.0)
            e.set(lev, ind, value=1.0)
            x0 = Scalar(grid)
            x0.assign(0.0)
            FluxToYVelocity(e, x0)
            x = AdjYVelToFlux(e)
            mag.set(lev, ind, value=InnerProductScalar(x0, x0))

            diff = Scalar(grid)
            diff.assign(x)
            diff._data[...] = x._data - x0._data  # full multi-level difference, matching C++ x-x0
            err2.set(lev, ind, value=InnerProductScalar(diff, diff))

    dumpFluxY("DUMP_mag", mag)
    dumpFluxY("DUMP_err2", err2)


if __name__ == "__main__":
    main()
