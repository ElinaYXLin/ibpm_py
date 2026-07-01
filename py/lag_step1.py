# lag_step1.py
#
# Python port of src/LagStep1.h
#
# Subclass of Motion, for a first-order lag filtered square pulse in
# dtheta (with area under pulse equal to AMP).

from __future__ import annotations

import math

from .motion import Motion
from .tangent_se2 import TangentSE2


class LagStep1(Motion):
    """First-order lag filtered square pulse in dtheta."""

    def __init__(self, AMP: float, PW: float, TAU: float, t0: float) -> None:
        self._AMP = AMP
        self._PW = PW
        self._TAU = TAU
        self._t0 = t0

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation for the lag-filtered pulse:
        (0, 0, theta(t), 0, 0, thetadot(t))."""
        pi = 4.0 * math.atan(1.0)
        AMPrad = self._AMP * (pi / 180)
        tdiff = time - self._t0
        if tdiff <= 0:
            theta = 0.0
            thetadot = 0.0
        elif 0 < tdiff <= self._PW:
            theta = tdiff + self._TAU * math.exp(-tdiff / self._TAU) - self._TAU
            thetadot = 1 - math.exp(-tdiff / self._TAU)
        elif tdiff > self._PW:
            theta = self._PW + self._TAU * math.exp(-tdiff / self._TAU) - self._TAU * math.exp(-(tdiff - self._PW) / self._TAU)
            thetadot = -1.0 * math.exp(-tdiff / self._TAU) + math.exp(-(tdiff - self._PW) / self._TAU)
        else:
            # NOTE(port): this branch is unreachable (the three conditions
            # above are exhaustive over the reals), matching the C++
            # `else { cerr << ...; exit(1); }` dead-code fallback. C++
            # terminates the whole process with exit(1); Python has no
            # direct equivalent that fits a library function, so this
            # raises RuntimeError instead, which is the idiomatic Python
            # way to halt execution with a clear error -- a judgment call
            # on translating `exit(1)` semantics.
            raise RuntimeError("LagStep1 (ERROR):  time is not compatible")
        return TangentSE2(0.0, 0.0, theta * AMPrad / self._PW, 0.0, 0.0, thetadot * AMPrad / self._PW)

    def clone(self) -> "LagStep1":
        return LagStep1(self._AMP, self._PW, self._TAU, self._t0)
