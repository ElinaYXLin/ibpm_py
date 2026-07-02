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


def _sprintf1(fmt: str, value: int) -> str:
    """Mimic C's `sprintf(buf, fmt, value)` for a single trailing integer
    argument, tolerating a mismatch between the number of `%`-conversion
    specifiers in `fmt` and the number of arguments supplied -- the way
    C's varargs `sprintf` does (extra arguments are silently ignored; too
    few are undefined behavior that in practice leaves the format string's
    extra specifiers un-substituted rather than crashing). Python's `%`
    operator instead requires an exact match, raising `TypeError` in both
    the "too many" case (e.g. a title string like "Check geometry" with no
    `%d` at all -- a legitimate caller-supplied title, not just a
    C++-vs-Python translation artifact) and the "too few" case. Both are
    caught here and `fmt` is returned unmodified: this is exactly correct
    for the "too many" case (matching sprintf's literal-copy behavior when
    there is nothing to substitute), and a safe, non-crashing placeholder
    for the "too few" case, which only arises in `doOutput` below in a
    branch where the resulting string is provably unused (see the
    docstring there).
    """
    try:
        return fmt % value
    except TypeError:
        return fmt


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
        # `sprintf(filename, _filename.c_str(), state.timestep)` and
        # `sprintf(title, _title.c_str(), state.timestep)` here, even
        # though: (a) when `_TecplotAllGrids` is true, `_filename` is
        # expected to contain a *second* `%d` (for the per-grid-level
        # index, filled in inside the loop below), so this first
        # `filename` sprintf is called with fewer arguments than format
        # specifiers -- undefined behavior in C++ that in practice leaves
        # the extra specifier un-substituted rather than crashing, and the
        # resulting value is provably unused in that branch anyway (it is
        # unconditionally overwritten inside the loop before any file is
        # written); and (b) `title` may legitimately contain *no* `%d` at
        # all (e.g. checkgeom.cc's literal `"Check geometry"`), in which
        # case sprintf just ignores the extra `state.timestep` argument and
        # copies the string unchanged. Python's `%` operator has no
        # equivalent tolerance for either case -- it raises `TypeError` for
        # both too few and too many arguments. Both are handled uniformly
        # by `_sprintf1` above (see its docstring for why returning `fmt`
        # unmodified is correct/safe in each case).
        filename = _sprintf1(self._filename, state.timestep)
        title = _sprintf1(self._title, state.timestep)
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
