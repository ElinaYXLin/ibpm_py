# motion.py
#
# Python port of src/Motion.h
#
# Specify a position and orientation as a function of time.
# Abstract base class: must be subclassed to be instantiated.

from __future__ import annotations

from abc import ABC, abstractmethod

from .tangent_se2 import TangentSE2


class Motion(ABC):
    """Abstract base class specifying a position/orientation as a function
    of time."""

    def isStationary(self) -> bool:
        """True if the body is moving (default is False)."""
        return False

    @abstractmethod
    def getTransformation(self, time: float) -> TangentSE2:
        """Return a Euclidean transformation and its velocity (an element
        of TSE(2)) at the specified time."""
        raise NotImplementedError

    @abstractmethod
    def clone(self) -> "Motion":
        """Return a copy of this Motion: subclasses must override."""
        raise NotImplementedError
