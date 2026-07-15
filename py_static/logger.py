# logger.py
#
# Python port of src/Logger.h / src/Logger.cc
#
# Maintain a list of output routines, and call them when specified.

from __future__ import annotations

from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from .base_flow import BaseFlow
    from .output import Output
    from .state import State


class _Entry:
    """Private helper mirroring the C++ `Logger::Entry` nested struct."""

    def __init__(self, output: "Output", numSkip: int) -> None:
        self.output = output
        self.numSkip = numSkip

    def shouldBeCalled(self, x: "State") -> bool:
        return x.timestep % self.numSkip == 0


class Logger:
    """Maintain a list of output routines, and call them when specified."""

    def __init__(self) -> None:
        self._outputs: List[_Entry] = []
        self._hasBeenInitialized: bool = False

    def addOutput(self, output: "Output", numSkip: int) -> None:
        """Add the specified Output to the list of output routines. The
        caller is responsible for allocating and deallocating memory for
        Output.

        numSkip: specifies how often to call the given routine (number of
        timesteps).
        """
        self._outputs.append(_Entry(output, numSkip))

    def doOutput(self, *args) -> bool:
        """Call all output routines needed at the current timestep,
        optionally making use of a baseflow.

        NOTE(port): collapses the two C++ overloads
            bool doOutput(const State& x);
            bool doOutput(const BaseFlow& q, const State& x);
        into one method dispatching on the number of positional arguments,
        matching the overload-collapse convention used throughout this
        port (see e.g. output_force.py).
        """
        if len(args) == 1:
            (x,) = args
            call_args = args
        elif len(args) == 2:
            _q, x = args
            call_args = args
        else:
            raise TypeError("Logger.doOutput: unsupported arguments")

        assert self._hasBeenInitialized
        successful = True
        for entry in self._outputs:
            if entry.shouldBeCalled(x):
                result = entry.output.doOutput(*call_args)
                successful = successful and result
        return successful

    def init(self) -> bool:
        """Initialize all of the output routines."""
        self._hasBeenInitialized = True
        successful = True
        for entry in self._outputs:
            result = entry.output.init()
            successful = successful and result
        return successful

    def cleanup(self) -> bool:
        """Clean up all of the output routines."""
        assert self._hasBeenInitialized
        successful = True
        for entry in self._outputs:
            result = entry.output.cleanup()
            successful = successful and result
        return successful
