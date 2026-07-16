#!/usr/bin/env python3
"""Python-side raw-value dump matching cpp/dump_fluxtox.cc, for
cross-validating py.vector_operations.FluxToXVelocity and XVelocityToFlux
against the C++ reference.

Run from the repo root:  python3 -m py.tests.cross_validation.python.dump_fluxtox
"""

from py.grid import Grid
from py.flux import Flux
from py.scalar import Scalar
from py.direction import Direction
from py.vector_operations import FluxToXVelocity, XVelocityToFlux


def dumpScalar(name, s):
    print(name)
    for lev in range(s.Ngrid()):
        for i in range(1, s.Nx()):
            for j in range(1, s.Ny()):
                print(repr(float(s._data[lev][i - 1, j - 1])))


def dumpFluxX(name, q):
    print(name)
    for lev in range(q.Ngrid()):
        for i in range(0, q.Nx() + 1):
            for j in range(0, q.Ny()):
                print(repr(q(lev, Direction.X, i, j)))


def computeDependencies(lev, i, j, q):
    q.assign(0.0)
    for l in range(q.Ngrid()):
        for ind in range(q.begin(Direction.X), q.end(Direction.X)):
            e = Flux(q.getGrid())
            e.assign(0.0)
            e.set(l, ind, value=1.0)
            x = Scalar(q.getGrid())
            x.assign(0.0)
            FluxToXVelocity(e, x)
            q.set(l, ind, value=x(lev, i, j))


def main():
    nx, ny, ngrid = 8, 8, 3
    length = 8.0
    grid = Grid(nx, ny, ngrid, length, 0.0, 0.0)

    u = Scalar(grid)
    u.assign(1.0)
    p = Flux(grid)
    p.assign(0.0)
    XVelocityToFlux(u, p)
    FluxToXVelocity(p, u)

    dumpFluxX("DUMP_ToFlux_u", p)
    dumpScalar("DUMP_roundtrip_u", u)

    pts = [(0, 1, 1), (0, 4, 4), (0, 7, 7), (1, 2, 3), (1, 4, 4), (2, 1, 1), (2, 4, 4)]
    for (lev, i, j) in pts:
        q = Flux(grid)
        computeDependencies(lev, i, j, q)
        print(f"DUMP_dep_{lev}_{i}_{j}")
        for l in range(q.Ngrid()):
            for ii in range(0, q.Nx() + 1):
                for jj in range(0, q.Ny()):
                    print(repr(q(l, Direction.X, ii, jj)))


if __name__ == "__main__":
    main()
