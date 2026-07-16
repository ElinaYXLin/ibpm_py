# output_restart.py
#
# Python port of src/OutputRestart.h / src/OutputRestart.cc
#
# Output routine for writing a restart file.

from __future__ import annotations

from typing import TYPE_CHECKING

from .output import Output

if TYPE_CHECKING:
    from .base_flow import BaseFlow
    from .state import State


class OutputRestart(Output):
    """Output routine for writing a restart file."""

    def __init__(self, formatString: str) -> None:
        """Constructor.

        formatString: filename in the standard printf format
        (e.g. "file%06d.bin"), where timestep will be substituted for %d.
        """
        self._formatString: str = formatString

    def doOutput(self, *args) -> bool:
        """Write the restart file.

        NOTE(port): collapses the two C++ overloads
            bool doOutput(const State& x);
            bool doOutput(const BaseFlow& q, const State& x);
        into one method dispatching on the number of positional arguments,
        matching the overload-collapse convention used throughout this
        port.
        """
        if len(args) == 1:
            (x,) = args
        elif len(args) == 2:
            _q, x = args
            # Currently no use for baseflow, but this method is defined for
            # future flexibility
        else:
            raise TypeError("OutputRestart.doOutput: unsupported arguments")

        # NOTE(port): C++ `sprintf(filename, _formatString.c_str(), x.timestep)`
        # is reproduced with Python's C-style `%` string formatting, which
        # supports the same printf-style numeric format specifiers (e.g.
        # "%06d") used by `_formatString` throughout this codebase.
        filename = self._formatString % x.timestep
        status = x.save(filename)
        return status

    def setFilename(self, formatString: str) -> None:
        """Change the filename for the output file."""
        self._formatString = formatString
