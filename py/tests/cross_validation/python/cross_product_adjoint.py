#!/usr/bin/env python3
"""Self-contained (no C++ build needed) cross-validation of
py.vector_operations.CrossProduct: checks the discrete adjoint identity

    < a, q1 x q2 > = < q1, q2 x a >

which is the property CrossProduct is explicitly designed to satisfy (see
the docstring in src/VectorOperations.h), at Ngrid = 1, 2, and 3 (i.e.
exercising the single-grid case and the multigrid border/interface/corner
code paths in FluxToXVelocity/FluxToYVelocity/VelocityToFlux).

This does not require building or running any C++ code -- it validates the
port purely against its own documented mathematical invariant. Run from the
repo root:  python3 -m py.tests.cross_validation.python.cross_product_adjoint
"""

import numpy as np

from py.grid import Grid
from py.flux import Flux
from py.scalar import Scalar
from py.direction import Direction
from py.vector_operations import CrossProduct


def make(nx, ny, ngrid, seed):
    grid = Grid(nx, ny, ngrid, 4.0, 0.0, 0.0)
    rng = np.random.default_rng(seed)

    def rand_flux():
        q = Flux(grid)
        q._data[...] = rng.standard_normal(q._data.shape)
        return q

    def rand_scalar():
        s = Scalar(grid)
        s._data[...] = rng.standard_normal(s._data.shape)
        return s

    return grid, rand_flux, rand_scalar


def ip_scalar(f, g):
    # full multigrid InnerProduct(Scalar,Scalar) per VectorOperations.cc
    nx, ny = f.Nx(), f.Ny()
    nx2, ny2 = f.NxExt(), f.NyExt()
    F = [f._data[l] for l in range(f.Ngrid())]
    G = [g._data[l] for l in range(g.Ngrid())]
    fv = lambda l, i, j: F[l][i - 1, j - 1]
    gv = lambda l, i, j: G[l][i - 1, j - 1]
    dx2 = f.Dx() * f.Dx()
    ip = 0.0
    for i in range(1, nx):
        for j in range(1, ny):
            ip += fv(0, i, j) * gv(0, i, j) * dx2
    for lev in range(1, f.Ngrid()):
        dx2 = f.Dx(lev) ** 2
        for (i, j) in [(nx2, ny2), (nx // 2 + nx2, ny2), (nx2, ny // 2 + ny2), (nx // 2 + nx2, ny // 2 + ny2)]:
            ip += fv(lev, i, j) * gv(lev, i, j) * dx2 * 15.0 / 16
        for j in range(ny2 + 1, ny // 2 + ny2):
            for i in (nx2, nx // 2 + nx2):
                ip += fv(lev, i, j) * gv(lev, i, j) * dx2 * 0.75
        for i in range(nx2 + 1, nx // 2 + nx2):
            for j in (ny2, ny // 2 + ny2):
                ip += fv(lev, i, j) * gv(lev, i, j) * dx2 * 0.75
        for i in range(1, nx2):
            for j in range(1, ny):
                ip += fv(lev, i, j) * gv(lev, i, j) * dx2
        for i in range(nx // 2 + nx2 + 1, nx):
            for j in range(1, ny):
                ip += fv(lev, i, j) * gv(lev, i, j) * dx2
        for i in range(nx2, nx // 2 + nx2 + 1):
            for j in range(1, ny2):
                ip += fv(lev, i, j) * gv(lev, i, j) * dx2
            for j in range(ny // 2 + ny2 + 1, ny):
                ip += fv(lev, i, j) * gv(lev, i, j) * dx2
    return ip


def ip_flux(p, q):
    nx, ny = p.Nx(), p.Ny()
    nx2, ny2 = p.NxExt(), p.NyExt()
    PX = [p._data[l, p.begin(Direction.X):p.end(Direction.X)].reshape(nx + 1, ny) for l in range(p.Ngrid())]
    PY = [p._data[l, p.begin(Direction.Y):p.end(Direction.Y)].reshape(nx, ny + 1) for l in range(p.Ngrid())]
    QX = [q._data[l, q.begin(Direction.X):q.end(Direction.X)].reshape(nx + 1, ny) for l in range(q.Ngrid())]
    QY = [q._data[l, q.begin(Direction.Y):q.end(Direction.Y)].reshape(nx, ny + 1) for l in range(q.Ngrid())]
    ip = 0.0
    for j in range(0, ny):
        for i in range(1, nx):
            ip += PX[0][i, j] * QX[0][i, j]
    for i in range(0, nx):
        for j in range(1, ny):
            ip += PY[0][i, j] * QY[0][i, j]
    for lev in range(1, p.Ngrid()):
        for j in range(ny2, ny // 2 + ny2):
            ip += PX[lev][nx2, j] * QX[lev][nx2, j] * 0.75
            ip += PX[lev][nx // 2 + nx2, j] * QX[lev][nx // 2 + nx2, j] * 0.75
        for j in range(0, ny):
            for i in range(1, nx2):
                ip += PX[lev][i, j] * QX[lev][i, j]
            for i in range(nx // 2 + nx2 + 1, nx):
                ip += PX[lev][i, j] * QX[lev][i, j]
        for i in range(nx2, nx // 2 + nx2 + 1):
            for j in range(0, ny2):
                ip += PX[lev][i, j] * QX[lev][i, j]
            for j in range(ny // 2 + ny2, ny):
                ip += PX[lev][i, j] * QX[lev][i, j]
    for lev in range(1, p.Ngrid()):
        for i in range(nx2, nx // 2 + nx2):
            ip += PY[lev][i, ny2] * QY[lev][i, ny2] * 0.75
            ip += PY[lev][i, ny // 2 + ny2] * QY[lev][i, ny // 2 + ny2] * 0.75
        for i in range(0, nx):
            for j in range(1, ny2):
                ip += PY[lev][i, j] * QY[lev][i, j]
            for j in range(ny // 2 + ny2 + 1, ny):
                ip += PY[lev][i, j] * QY[lev][i, j]
        for j in range(ny2, ny // 2 + ny2 + 1):
            for i in range(0, nx2):
                ip += PY[lev][i, j] * QY[lev][i, j]
            for i in range(nx // 2 + nx2, nx):
                ip += PY[lev][i, j] * QY[lev][i, j]
    return ip


def main():
    for ngrid in (1, 2, 3):
        grid, rand_flux, rand_scalar = make(8, 8, ngrid, seed=ngrid)
        a = rand_scalar()
        q1 = rand_flux()
        q2 = rand_flux()
        lhs = ip_scalar(a, CrossProduct(q1, q2))
        rhs = ip_flux(q1, CrossProduct(q2, a))
        rel = abs(lhs - rhs) / max(abs(lhs), 1e-30)
        print(f"Ngrid={ngrid}: <a,q1xq2>={lhs:.10g}  <q1,q2xa>={rhs:.10g}  reldiff={rel:.3g}")
        assert rel < 1e-10, f"adjoint identity failed at Ngrid={ngrid}"
    print("ALL PASS")


if __name__ == "__main__":
    main()
