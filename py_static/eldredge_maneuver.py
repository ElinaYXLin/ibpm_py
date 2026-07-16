# eldredge_maneuver.py
#
# Python port of src/EldredgeManeuver.h
#
# Subclass of Motion, for the canonical maneuver described by Eldredge.

from __future__ import annotations

import math

from .motion import Motion
from .tangent_se2 import TangentSE2


class EldredgeManeuver(Motion):
    """Eldredge's canonical smoothed-ramp maneuver, applied to theta."""

    def __init__(self, AMP: float, a: float, t1: float, t2: float, t3: float, t4: float) -> None:
        self._AMP = AMP
        self._a = a
        self._t1 = t1
        self._t2 = t2
        self._t3 = t3
        self._t4 = t4
        tt = (t2 + t3) / 2.0
        self._maxG = math.log(
            (math.cosh(a * (tt - t1)) * math.cosh(a * (tt - t4)))
            / (math.cosh(a * (tt - t2)) * math.cosh(a * (tt - t3)))
        )

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation for the Eldredge maneuver, applied to
        theta: (0, 0, theta(t), 0, 0, thetadot(t))."""
        pi = 4.0 * math.atan(1.0)
        AMPrad = self._AMP * (pi / 180)
        a = self._a
        t1d = time - self._t1
        t2d = time - self._t2
        t3d = time - self._t3
        t4d = time - self._t4
        Garg = (math.cosh(a * t1d) * math.cosh(a * t4d)) / (math.cosh(a * t2d) * math.cosh(a * t3d))
        G = math.log(Garg)
        dGdt = a * (math.tanh(a * t1d) - math.tanh(a * t2d) - math.tanh(a * t3d) + math.tanh(a * t4d))

        return TangentSE2(0.0, 0.0, G * AMPrad / self._maxG, 0.0, 0.0, dGdt * AMPrad / self._maxG)

    def clone(self) -> "EldredgeManeuver":
        return EldredgeManeuver(self._AMP, self._a, self._t1, self._t2, self._t3, self._t4)
