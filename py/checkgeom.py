# checkgeom.py
#
# Python port of src/checkgeom.cc
#
# checkgeom - read in geometry file and write corresponding Tecplot file.
#
# Use this utility to check the spacing of the boundary points. If the
# boundary points are too sparse, the boundary will look "leaky".
#
# NOTE: If the points are too finely spaced, the solution will not converge
# when used in a timestepper.
#
# NOTE(port): see the module-level judgment-call note in ibpm.py for the
# `main(argv=...)` entry-point convention and the `cout`/`cerr` -> `print`/
# `sys.stderr` mapping used throughout this port; the same choices apply
# here.

from __future__ import annotations

import sys
from typing import List, Optional

from .base_flow import BaseFlow
from .geometry import Geometry
from .grid import Grid
from .output_tecplot import OutputTecplot
from .parm_parser import ParmParser
from .regularizer import Regularizer
from .state import State


def main(argv: Optional[List[str]] = None) -> int:
    """Read in a geometry file and write the corresponding Tecplot file."""
    if argv is None:
        argv = sys.argv

    print("Check geometry\n")

    parser = ParmParser(len(argv), argv)
    helpFlag = parser.getFlag("h", "print this help message and exit")
    nx = parser.getInt("nx", "number of gridpoints in x-direction", 200)
    ny = parser.getInt("ny", "number of gridpoints in y-direction", 200)
    ngrid = parser.getInt("ngrid", "number of grid levels for multi-domain scheme", 1)
    length = parser.getDouble("length", "length of finest domain in x-dir", 4.0)
    xOffset = parser.getDouble("xoffset", "x-coordinate of left edge of finest domain", -2.0)
    yOffset = parser.getDouble("yoffset", "y-coordinate of bottom edge of finest domain", -2.0)
    geomFile = parser.getString("geom", "filename for reading geometry", "ibpm.geom")
    outFileName = parser.getString("o", "filename for writing Tecplot file", "")

    if (not parser.inputIsValid()) or helpFlag:
        parser.printUsage(sys.stderr)
        sys.exit(1)

    # Setup grid
    print(
        "Grid parameters:\n"
        f"  nx      {nx}\n"
        f"  ny      {ny}\n"
        f"  ngrid   {ngrid}\n"
        f"  length  {length}\n"
        f"  xoffset {xOffset}\n"
        f"  yoffset {yOffset}"
    )
    grid = Grid(nx, ny, ngrid, length, xOffset, yOffset)

    # Setup geometry
    geom = Geometry()
    print(f"Reading geometry from file {geomFile}")
    if geom.load(geomFile):
        print(f"  {geom.getNumPoints()} points on the boundary")
    else:
        print("  There were errors reading the geometry file.")
        return 1

    if outFileName == "":
        return 0

    regularizer = Regularizer(grid, geom)
    regularizer.update()
    x = State(grid, geom.getNumPoints())

    # set the boundary force to 1 (in x- and y- directions)
    x.f.assign(1.0)
    # NOTE(port): C++ `x.q = regularizer.toFlux(x.f);` invokes
    # `Flux::operator=`, copying into the existing `x.q` Flux in place;
    # `.assign(...)` is used here (rather than rebinding `x.q =`) to match
    # that copy-assignment convention, consistent with the rest of this
    # port (see e.g. state.py / scalar.py / flux.py `.assign()`).
    x.q.assign(regularizer.toFlux(x.f))

    tecplot = OutputTecplot(outFileName, "Check geometry", False)
    q_potential = BaseFlow(grid, 0.0, 0.0)
    tecplot.doOutput(q_potential, x)

    return 0


if __name__ == "__main__":
    sys.exit(main())
