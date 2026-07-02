# output.py
#
# Python port of src/Output.h
#
# Abstract base class for output routines. Subclasses provide a callback
# routine for writing output in various forms.

from __future__ import annotations

from abc import ABC, abstractmethod


class Output(ABC):
    """Abstract base class for output routines."""

    def init(self) -> bool:
        """Provide initialization, if needed (e.g. opening a file).
        Returns true if successful."""
        return True

    def cleanup(self) -> bool:
        """Clean up, if needed (e.g. close a file). Returns true if
        successful."""
        return False

    @abstractmethod
    def doOutput(self, *args) -> bool:
        """Callback for performing the actual output from a State object,
        or from a BaseFlow and a State object.

        NOTE(port): C++ declares two pure-virtual overloads,
            virtual bool doOutput(const State& x) = 0;
            virtual bool doOutput(const BaseFlow& q, const State& x) = 0;
        Python has no method overloading, so both are collapsed into this
        single abstract `doOutput`, dispatching on the number/type of
        `*args`; each concrete subclass below implements that dispatch
        (matching the overload-collapse convention used elsewhere in this
        port, e.g. vector_operations.Curl).
        """
        raise NotImplementedError
