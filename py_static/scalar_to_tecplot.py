# scalar_to_tecplot.py
#
# Python port of src/ScalarToTecplot.h / src/ScalarToTecplot.cc
#
# Write one or more Scalar fields to an ASCII Tecplot file.

from __future__ import annotations

import sys
from typing import List, Sequence, Union

import numpy as np

from .scalar import Scalar

ScalarOrScalars = Union[Scalar, Sequence[Scalar]]
StrOrStrs = Union[str, Sequence[str]]


class _VarList:
    """Private helper mirroring the C++ `ScalarToTecplot.cc::VarList` class:
    a parallel list of Scalar fields and their variable names."""

    def __init__(self) -> None:
        self._vars: List[Scalar] = []
        self._names: List[str] = []

    def addVariable(self, var: Scalar, varName: str) -> None:
        if self.getNumVars() > 0:
            assert self._vars[0].getGrid().isEqualTo(var.getGrid())
        self._vars.append(var)
        self._names.append(varName)

    def getNumVars(self) -> int:
        return len(self._vars)

    def getName(self, i: int) -> str:
        return self._names[i]

    def getVariable(self, i: int) -> Scalar:
        return self._vars[i]


def _writeTecplotFileASCII(filename: str, title: str, varlist: _VarList, lev: int) -> bool:
    numVars = varlist.getNumVars()
    assert numVars > 0

    # Get grid information
    grid = varlist.getVariable(0).getGrid()
    assert lev < grid.Ngrid()
    nx = grid.Nx()
    ny = grid.Ny()

    # Write the header for the Tecplot file
    sys.stderr.write(f"Writing Tecplot file {filename}\n")
    try:
        fp = open(filename, "w")
    except OSError:
        return False

    with fp:
        fp.write(f'TITLE = "{title}"\n')
        fp.write("VARIABLES = ")
        for i in range(numVars):
            fp.write(f'"{varlist.getName(i)}" ')
        fp.write("\n")
        fp.write('ZONE T="Rectangular zone"\n')
        fp.write(f"I={nx - 1}, J={ny - 1}, K=1, ZONETYPE=Ordered\n")
        fp.write("DATAPACKING=POINT\n")
        fp.write("DT=(")
        for i in range(numVars):
            fp.write("SINGLE ")
        fp.write(")\n")

        # Write the data
        #
        # NOTE(port): the C++ (j outer, i inner) double loop over interior
        # points, writing all `numVars` values at each (i,j) before moving
        # on, is vectorized below: the per-variable interior arrays are
        # stacked into one (nx-1, ny-1, numVars) array, transposed/reshaped
        # into the C++ (j,i) iteration order, and written in one
        # `np.savetxt` call, instead of a manual per-point/-variable Python
        # loop. `fmt` includes a trailing space per value (matching the
        # C++ `fprintf(fp, "%.5e ", ...)` per value, followed by a separate
        # "\n").
        data = np.stack(
            [varlist.getVariable(ind)[lev] for ind in range(numVars)], axis=-1
        )  # shape (nx-1, ny-1, numVars); data[i-1, j-1, ind]
        ordered = np.transpose(data, (1, 0, 2)).reshape(-1, numVars)
        # NOTE(port): C++ writes `fprintf(fp, "%.5e ", ...)` once per value
        # (value + a single trailing space) with no additional separator, so
        # the per-value tokens are concatenated with "" (not " ".join, which
        # would insert a second space between values).
        fmt = "".join(["%.5e "] * numVars)
        np.savetxt(fp, ordered, fmt=fmt)

    return True


def ScalarToTecplot(
    var: ScalarOrScalars,
    varName: StrOrStrs,
    filename: str,
    title: str,
    lev: int = 0,
) -> bool:
    """Write one or more Scalar fields to an ASCII Tecplot file.

    NOTE(port): collapses the four C++ overloads
        bool ScalarToTecplot(const Scalar* var, string varName, string filename, string title);
        bool ScalarToTecplot(vector<const Scalar*> varVec, vector<string> varNameVec, string filename, string title);
        bool ScalarToTecplot(const Scalar* var, string varName, string filename, string title, int lev);
        bool ScalarToTecplot(vector<const Scalar*> varVec, vector<string> varNameVec, string filename, string title, int lev);
    into one function: `var`/`varName` may each be a single Scalar/str
    (matching the pointer-taking overloads) or a sequence of them (matching
    the vector-taking overloads), and `lev` defaults to 0 (matching the two
    C++ overloads that omit the `lev` argument).
    """
    if isinstance(var, Scalar):
        varVec: List[Scalar] = [var]
        varNameVec: List[str] = [varName]  # type: ignore[list-item]
    else:
        varVec = list(var)
        varNameVec = list(varName)  # type: ignore[arg-type]

    assert len(varVec) > 0
    assert len(varVec) == len(varNameVec)

    # Get grid dimensions
    grid = varVec[0].getGrid()
    nx = grid.Nx()
    ny = grid.Ny()
    ngrid = grid.Ngrid()
    assert lev < ngrid

    # Calculate the variables for output
    # Calculate the grid
    #
    # NOTE(port): the C++ triple loop (over _lev, i, j) assigning
    # x(_lev,i,j) = grid.getXEdge(_lev,i) and y(_lev,i,j) = grid.getYEdge(_lev,j)
    # is vectorized below with numpy broadcasting per level, instead of a
    # manual per-point Python loop. The level loop itself is kept (Ngrid is
    # small, and getXEdge/getYEdge depend on the per-level offset/spacing).
    x = Scalar(grid)
    y = Scalar(grid)
    for _lev in range(ngrid):
        dx = grid.Dx(_lev)
        i_arr = np.arange(1, nx)
        j_arr = np.arange(1, ny)
        x_edge = grid.getXEdge(_lev, 0) + i_arr * dx
        y_edge = grid.getYEdge(_lev, 0) + j_arr * dx
        x._data[_lev] = x_edge[:, None]
        y._data[_lev] = y_edge[None, :]

    # Store pointers to variables and corresponding names in vectors
    varlist = _VarList()
    varlist.addVariable(x, "x")
    varlist.addVariable(y, "y")
    for i in range(len(varVec)):
        varlist.addVariable(varVec[i], varNameVec[i])

    # Write the Tecplot file
    status = _writeTecplotFileASCII(filename, title, varlist, lev)
    return status
