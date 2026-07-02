# output_tecplot.py
#
# Python port of src/OutputTecplot.h / src/OutputTecplot.cc
#
# Output routines for writing ASCII Tecplot files.

from __future__ import annotations

from typing import TYPE_CHECKING

from .output import Output
from .scalar import Scalar
from .scalar_to_tecplot import ScalarToTecplot
from .vector_operations import FluxToVelocity

if TYPE_CHECKING:
    from .base_flow import BaseFlow
    from .state import State


class OutputTecplot(Output):
    """Output routines for writing ASCII Tecplot files."""

    def __init__(self, filename: str, title: str, TecplotAllGrids: bool) -> None:
        """Constructor.

        filename: filename in the standard printf format (e.g.
        "file%06d.plt"), where timestep will be supplied.
        title: title in the standard printf format.
        """
        self._filename: str = filename
        self._title: str = title
        self._TecplotAllGrids: bool = TecplotAllGrids

    def doOutput(self, *args) -> bool:
        """Write the Tecplot file.

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
            raise TypeError("OutputTecplot.doOutput: unsupported arguments")

        # Add timestep to filename and title
        #
        # NOTE(port) -- judgment call: C++ unconditionally does
        # `sprintf(filename, _filename.c_str(), state.timestep)` here, even
        # when `_TecplotAllGrids` is true and `_filename` is expected to
        # contain a *second* `%d` (for the per-grid-level index, filled in
        # inside the loop below) -- calling sprintf with fewer arguments
        # than format specifiers is undefined behavior in C++ (typically it
        # reads a garbage/stray value rather than crashing), and the
        # resulting `filename` is provably unused in that branch anyway (it
        # is unconditionally overwritten inside the loop before any file is
        # written). Python's `%` operator has no such "UB but survives"
        # behavior -- it raises `TypeError` immediately for too few
        # arguments. Since this value is dead code exactly when the
        # mismatch would occur, the exception is caught and `filename` left
        # as `None` in that case; it is unconditionally reassigned before
        # use in the `_TecplotAllGrids` branch below, and used as-is
        # (successfully formatted) in the non-`_TecplotAllGrids` branch.
        try:
            filename = self._filename % state.timestep
        except TypeError:
            filename = None
        title = self._title % state.timestep
        status = False
        grid = state.omega.getGrid()

        # Calculate velocities
        u = Scalar(state.omega.getGrid())
        v = Scalar(state.omega.getGrid())
        FluxToVelocity(state.q, u, v)

        # Create list of Scalar fields
        varVec = [u, v, state.omega]
        varNameVec = ["u", "v", "Vorticity"]

        # Write the tecplot file
        if self._TecplotAllGrids:
            status = True
            for i in range(grid.Ngrid()):
                filename = self._filename % (state.timestep, i)
                print(filename)
                status = status and ScalarToTecplot(varVec, varNameVec, filename, title, i)
        else:
            status = ScalarToTecplot(varVec, varNameVec, filename, title)
        return status

    def setFilename(self, filename: str) -> None:
        """Change the filename for the output file."""
        self._filename = filename

    def setTitle(self, title: str) -> None:
        """Change the title for the output file."""
        self._title = title
