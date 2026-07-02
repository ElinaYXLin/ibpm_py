# output_probes.py
#
# Python port of src/OutputProbes.h / src/OutputProbes.cc
#
# Write velocities, fluxes, and vorticity, at given probe locations, to
# files. Each probe has a corresponding output file. All probes are
# supposed to be located at the interior nodes at the finest grid level
# (level 0). Probes are labelled as Probe 1, 2, ... . Probe information
# (probe #, position) is stored in a separate file.

from __future__ import annotations

import sys
from typing import TYPE_CHECKING, List, Optional, TextIO

from .direction import Direction
from .grid import Grid
from .output import Output
from .scalar import Scalar
from .vector_operations import FluxToVelocity

if TYPE_CHECKING:
    from .base_flow import BaseFlow
    from .state import State


class _Probe:
    """Private helper mirroring the C++ `OutputProbes::Probe` nested class."""

    def __init__(self, ii: int, jj: int) -> None:
        self.i = ii
        self.j = jj
        self.fp: Optional[TextIO] = None


class OutputProbes(Output):
    """Write velocities, fluxes, and vorticity, at given probe locations,
    to files."""

    # all probes at finest grid level
    _lev = 0
    # two-dimensional domain for now
    _dimen = 2

    def __init__(self, filename: str, grid: Grid) -> None:
        """Constructor.

        filename: filename to which probe data will be written. For
        instance, if filename = "out/probe%02d.dat" then the following
        files will be created:
          out/probe00.dat   (description of probe locations)
          out/probe01.dat   (first probe)
          out/probe02.dat   (second probe)
          ...
        """
        self._filename: str = filename
        # NOTE(port): C++ stores `Grid _grid` by value, so `_grid(grid)`
        # in the member-initializer list copy-constructs a private copy.
        # This port instead aliases the same Grid object passed in, matching
        # the convention already established for storing a passed-in Grid
        # elsewhere in this port (see Field.__init__ in field.py, which does
        # the same for a plain `Grid` argument) rather than introducing a
        # new deep-copy helper Grid lacks.
        self._grid: Grid = grid
        self._hasBeenInitialized: bool = False
        self._probes: List[_Probe] = []

    def init(self) -> bool:
        """Write a file with description of probe locations (this file has
        index 0). Also, open a file for each probe. If a file with the same
        name is already present, it is overwritten. Returns true if
        successful."""
        if not self._writeSummaryFile():
            return False

        # open files for each probe to output data
        for n in range(len(self._probes)):
            name = self._filename % (n + 1)
            try:
                self._probes[n].fp = open(name, "w")
            except OSError:
                return False
        self._hasBeenInitialized = True
        return True

    def cleanup(self) -> bool:
        """Close all the files. Returns true if successful."""
        status = True
        for probe in self._probes:
            if probe.fp is not None:
                probe.fp.close()
        return status

    def doOutput(self, *args) -> bool:
        """Write velocities u, v, fluxes q.x, q.y and vorticity omega for
        each probe, to the corresponding file with name (filename +
        probe#).

        NOTE(port): collapses the two C++ overloads
            bool doOutput(const State& x);
            bool doOutput(const BaseFlow& q, const State& x);
        into one method dispatching on the number of positional arguments,
        matching the overload-collapse convention used throughout this
        port.
        """
        if len(args) == 1:
            (state,) = args
        elif len(args) == 2:
            _q, state = args
            # Currently no use for baseflow, but this method is defined for
            # future flexibility
        else:
            raise TypeError("OutputProbes.doOutput: unsupported arguments")

        # TODO: Unnecessary to transform velocity fields everywhere, when
        # only a few probe points will be used
        u = Scalar(self._grid)
        v = Scalar(self._grid)
        FluxToVelocity(state.q, u, v)

        # Write u, v, qx, qy, omega, all at gridpoint/edge (i, j), for each
        # probe, in separate files
        # TODO: Why store qx and u? Seems like this is redundant.
        for probe in self._probes:
            assert probe.fp is not None
            i = probe.i
            j = probe.j
            probe.fp.write(
                "%5d %.5e %.14e %.14e %.14e %.14e %.14e\n"
                % (
                    state.timestep,
                    state.time,
                    u(self._lev, i, j),
                    v(self._lev, i, j),
                    state.q(self._lev, Direction.X, i, j),
                    state.q(self._lev, Direction.Y, i, j),
                    state.omega(self._lev, i, j),
                )
            )
            probe.fp.flush()

        return True

    def addProbeByIndex(self, i: int, j: int) -> None:
        """Add a probe by specifying its gridpoint indices."""
        if self._hasBeenInitialized:
            print(
                "Error: Addition of probes is allowed only "
                "before initialization."
            )
            sys.exit(1)

        if i < 1 or j < 1 or i > self._grid.Nx() - 1 or j > self._grid.Ny() - 1:
            print(
                f"Warning: invalid probe position: ({i},{j})\n"
                "Probes should be located at the inner nodes "
                "at the finest grid level."
            )
            sys.exit(1)
        probe = _Probe(i, j)
        self._probes.append(probe)

    def addProbeByPosition(self, xcord: float, ycord: float) -> None:
        """Add a probe by specifying its absolute coordinates."""
        i = self._grid.getXGridIndex(xcord)
        j = self._grid.getYGridIndex(ycord)
        self.addProbeByIndex(i, j)

    def addProbe(self, *args) -> None:
        """Add a probe by specifying either its gridpoint indices (two
        ints) or its absolute coordinates (two floats).

        NOTE(port): collapses the two C++ overloads
            void addProbe(int i, int j);
            void addProbe(double xcord, double ycord);
        into one method dispatching on argument type, matching the
        overload-collapse convention used throughout this port.
        """
        a, b = args
        if isinstance(a, int) and isinstance(b, int):
            self.addProbeByIndex(a, b)
        else:
            self.addProbeByPosition(a, b)

    def print(self) -> None:
        """Print out probe locations (by grid indices), for debugging."""
        print("\n-- Probe locations: by grid indices --")
        for n, probe in enumerate(self._probes):
            print(f" {n + 1}-th probe: grid index ( {probe.i}, {probe.j} )")
        print("---------- end of probe list --------")

    def getNumProbes(self) -> int:
        """Return the number of probes."""
        return len(self._probes)

    def getProbeIndexX(self, index: int) -> int:
        """Return the gridpoint index i of the corresponding probe."""
        assert index <= len(self._probes) and index >= 1
        return self._probes[index - 1].i

    def getProbeIndexY(self, index: int) -> int:
        """Return the gridpoint index j of the corresponding probe."""
        assert index <= len(self._probes)
        assert index >= 1
        return self._probes[index - 1].j

    def getProbeCoordX(self, index: int) -> float:
        """Return the gridpoint x coordinate of the corresponding probe."""
        assert index <= len(self._probes) and index >= 1
        return self._grid.getXEdge(self._lev, self._probes[index - 1].i)

    def getProbeCoordY(self, index: int) -> float:
        """Return the gridpoint y coordinate of the corresponding probe."""
        assert index <= len(self._probes) and index >= 1
        return self._grid.getYEdge(self._lev, self._probes[index - 1].j)

    def _writeSummaryFile(self) -> bool:
        """Write summary of probe info, to file with index 0."""
        fname = self._filename % 0
        try:
            fp = open(fname, "w")
        except OSError:
            return False
        # write: probe#, probe grid indices i, j, probe coordinates x, y
        with fp:
            for n in range(len(self._probes)):
                i = self.getProbeIndexX(n + 1)
                j = self.getProbeIndexY(n + 1)
                x = self.getProbeCoordX(n + 1)
                y = self.getProbeCoordY(n + 1)
                fp.write("%2d %3d %d %.5e %.5e\n" % (n + 1, i, j, x, y))
        return True
