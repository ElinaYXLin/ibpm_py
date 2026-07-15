"""Test-only helper implementations of VectorOperations.cc's
InnerProduct(Scalar,Scalar) and InnerProduct(Flux,Flux), transcribed
directly from src/VectorOperations.cc lines 249-310 (Scalar) and 345-411
(Flux).

These are NOT part of the production port (py/vector_operations.py does not
define InnerProduct -- see that module's header note on scope: it is ported
incrementally, only as dependents need pieces of it, and nothing in the
production port currently calls InnerProduct(Scalar,Scalar) or
InnerProduct(Flux,Flux) directly). They exist here solely so the
cross-validation dump scripts in this directory can reproduce the exact
weighted multigrid inner product that test/CheckCurl.cc, test/CheckAdjoint.cc
and their C++ dump-driver counterparts use, without modifying any production
port file.

If a production file ever needs InnerProduct(Scalar,Scalar) or
InnerProduct(Flux,Flux), port it into py/vector_operations.py directly rather
than importing from here.
"""

from py.direction import Direction


def InnerProductFlux(p, q):
    """Port of VectorOperations.cc InnerProduct(const Flux&, const Flux&)."""
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


def InnerProductScalar(f, g):
    """Port of VectorOperations.cc InnerProduct(const Scalar&, const Scalar&)."""
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
