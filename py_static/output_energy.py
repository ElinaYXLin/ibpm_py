# output_energy.py
#
# Python port of src/OutputEnergy.h / src/OutputEnergy.cc

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, TextIO

from .output import Output
from .vector_operations import InnerProduct

if TYPE_CHECKING:
    from .base_flow import BaseFlow
    from .state import State


class OutputEnergy(Output):

    def __init__(self, filename: str) -> None:
        """Constructor.

        filename: file to which force data will be written.
        """
        self._filename: str = filename
        # NOTE(port): see the corresponding note in OutputForce -- C++
        # `_fp` is uninitialized until `init()` is called; initialized to
        # None here instead of left undefined.
        self._fp: Optional[TextIO] = None

    def init(self) -> bool:
        """Open the file for writing. If a file with the same name is
        already present, it is overwritten. Returns true if successful."""
        try:
            self._fp = open(self._filename, "w")
        except OSError:
            return False
        return True

    def cleanup(self) -> bool:
        """Close the file. Returns true if successful."""
        status = True
        if self._fp is not None:
            self._fp.close()
        return status

    def doOutput(self, *args) -> bool:
        """Write data to the energy file.

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
            raise TypeError("OutputEnergy.doOutput: unsupported arguments")

        energy = 0.5 * InnerProduct(x.q, x.q)

        if self._fp is None:
            return False
        self._fp.write("%5d %.5e %.5e\n" % (x.timestep, x.time, energy))
        self._fp.flush()
        return True
