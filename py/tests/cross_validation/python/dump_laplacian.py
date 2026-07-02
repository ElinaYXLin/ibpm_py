#!/usr/bin/env python3
"""Python-side raw-value dump matching cpp/dump_laplacian.cc, for
cross-validating py.vector_operations.Laplacian and Scalar.coarsify()
against the C++ reference.

Run from the repo root:  python3 -m py.tests.cross_validation.python.dump_laplacian
"""

from py.grid import Grid
from py.scalar import Scalar
from py.vector_operations import Laplacian


def equal(u, v):
    tol = 1e-12
    err = 0.0
    for lev in range(u.Ngrid()):
        for i in range(1, u.Nx()):
            for j in range(1, u.Ny()):
                err += abs(u._data[lev][i - 1, j - 1] - v._data[lev][i - 1, j - 1])
    return err < tol


def checkLCEqualsCLC(grid, lev, i, j):
    u = Scalar(grid)
    u.assign(0.0)
    u._data[lev][i - 1, j - 1] = 1.0
    u.coarsify()
    LCu = Scalar(grid)
    Laplacian(u, LCu)
    CLCu = Scalar(LCu)
    CLCu.coarsify()
    return equal(CLCu, LCu)


def main():
    nx, ny, ngrid = 8, 8, 2
    length = 8.0
    grid = Grid(nx, ny, ngrid, length, -length / 2, -length * ny / nx / 2)

    status = Scalar(grid)
    status.assign(2.0)
    for lev in range(ngrid):
        for i in range(1, nx):
            for j in range(1, ny):
                status._data[lev][i - 1, j - 1] = 0.0 if checkLCEqualsCLC(grid, lev, i, j) else 1.0

    print("DUMP_status")
    for lev in range(ngrid):
        for i in range(1, nx):
            for j in range(1, ny):
                print(repr(float(status._data[lev][i - 1, j - 1])))


if __name__ == "__main__":
    main()
