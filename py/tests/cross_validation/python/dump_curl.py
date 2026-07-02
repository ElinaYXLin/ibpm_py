#!/usr/bin/env python3
"""Python-side raw-value dump matching cpp/dump_curl.cc, for cross-validating
py.vector_operations.Curl against the C++ reference (see the .cc file's
header comment for build/run instructions and py/tests/cross_validation/
README.md for the overall workflow).

Run from the repo root:  python3 -m py.tests.cross_validation.python.dump_curl
"""

from py.grid import Grid
from py.flux import Flux
from py.scalar import Scalar
from py.direction import Direction
from py.vector_operations import Curl
from .inner_products import InnerProductFlux, InnerProductScalar


def AdjointOfScalarCurl(q_, f):
    for lev in range(q_.Ngrid()):
        for i in range(1, q_.Nx()):
            for j in range(1, q_.Ny()):
                e = Scalar(q_.getGrid())
                e.assign(0.0)
                e._data[lev][i - 1, j - 1] = 1.0
                curlE = Flux(q_.getGrid())
                Curl(e, curlE)
                f._data[lev][i - 1, j - 1] = InnerProductFlux(q_, curlE) / InnerProductScalar(e, e)


def dumpScalar(name, s):
    print(name)
    for lev in range(s.Ngrid()):
        for i in range(1, s.Nx()):
            for j in range(1, s.Ny()):
                print(repr(float(s._data[lev][i - 1, j - 1])))


def dumpFlux(name, q_):
    print(name)
    for lev in range(q_.Ngrid()):
        for i in range(0, q_.Nx() + 1):
            for j in range(0, q_.Ny()):
                print(repr(q_(lev, Direction.X, i, j)))
        for i in range(0, q_.Nx()):
            for j in range(0, q_.Ny() + 1):
                print(repr(q_(lev, Direction.Y, i, j)))


def main():
    nx, ny, ngrid = 8, 8, 2
    length = 8.0
    grid = Grid(nx, ny, ngrid, length, -length / 2, -length * ny / nx / 2)
    u = Scalar(grid)
    v = Scalar(grid)
    q = Flux(grid)
    u.assign(0.0)
    v.assign(0.0)
    q.assign(0.0)

    for lev in range(ngrid):
        for i in range(1, nx):
            for j in range(1, ny):
                u._data[lev][i - 1, j - 1] = grid.getXEdge(lev, i)
                v._data[lev][i - 1, j - 1] = grid.getYEdge(lev, j)

    for lev in range(ngrid):
        dx = grid.Dx(lev)
        for i in range(0, nx + 1):
            for j in range(0, ny):
                q.set(lev, Direction.X, i, j, value=-2 * q.y(lev, Direction.X, j) * q.y(lev, Direction.X, j) * dx)
        for i in range(0, nx):
            for j in range(0, ny + 1):
                q.set(lev, Direction.Y, i, j, value=q.x(lev, Direction.Y, i) * q.x(lev, Direction.Y, i) * dx)

    curlU = Curl(u)
    curlV = Curl(v)
    omega = Curl(q)

    omega2 = Scalar(grid)
    omega2.assign(0.0)
    AdjointOfScalarCurl(q, omega2)
    err = Scalar(grid)
    err.assign(omega)
    err._data[...] = omega._data - omega2._data

    dumpScalar("DUMP_u", u)
    dumpFlux("DUMP_curlU", curlU)
    dumpScalar("DUMP_v", v)
    dumpFlux("DUMP_curlV", curlV)
    dumpFlux("DUMP_q", q)
    dumpScalar("DUMP_omega", omega)
    dumpScalar("DUMP_omega2", omega2)
    dumpScalar("DUMP_err", err)


if __name__ == "__main__":
    main()
