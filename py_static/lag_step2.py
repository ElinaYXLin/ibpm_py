# lag_step2.py
#
# Python port of src/LagStep2.h
#
# Subclass of Motion, for a second-order lag filtered square pulse in
# dtheta (with area under pulse equal to AMP).

from __future__ import annotations

import math

from .motion import Motion
from .tangent_se2 import TangentSE2


class LagStep2(Motion):
    """Second-order lag filtered square pulse in dtheta."""

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
        TAU = self._TAU
        PW = self._PW
        if tdiff <= 0:
            theta = 0.0
            thetadot = 0.0
        elif 0 < tdiff <= PW:
            theta = tdiff + (2 * TAU + tdiff) * math.exp(-tdiff / TAU) - 2 * TAU
            thetadot = 1.0 - (1.0 + tdiff / TAU) * math.exp(-tdiff / TAU)
        elif tdiff > PW:
            theta = (
                PW
                + (2 * TAU + tdiff) * math.exp(-tdiff / TAU)
                + (PW - 2 * TAU) * math.exp(-(tdiff - PW) / TAU)
                - tdiff * math.exp(-(tdiff - PW) / TAU)
            )
            thetadot = -(1.0 + tdiff / TAU) * math.exp(-tdiff / TAU) + (1.0 + (tdiff - PW) / TAU) * math.exp(-(tdiff - PW) / TAU)
        else:
            # NOTE(port): unreachable, matching C++'s dead `exit(1)`
            # fallback -- see the identical note in lag_step1.py.
            raise RuntimeError("LagStep2 (ERROR):  time is not compatible")
        return TangentSE2(0.0, 0.0, theta * AMPrad / PW, 0.0, 0.0, thetadot * AMPrad / PW)

    def clone(self) -> "LagStep2":
        return LagStep2(self._AMP, self._PW, self._TAU, self._t0)
