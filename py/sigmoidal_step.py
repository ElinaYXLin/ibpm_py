# sigmoidal_step.py
#
# Python port of src/SigmoidalStep.h
#
# Subclass of Motion, for a sigmoidal step in angle-of-attack.

from __future__ import annotations

import math

from .motion import Motion
from .tangent_se2 import TangentSE2


class SigmoidalStep(Motion):
    """A sigmoidal pitch-up centered about the origin:
        alpha(t) = AMP * 1/2 * (1 + erf(12*t/DUR - 6))
    """

    def __init__(self, AMP: float, DUR: float, t0: float) -> None:
        self._AMP = AMP
        self._DUR = DUR
        self._t0 = t0

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation for the sigmoidal step:
        (0, 0, theta(t), 0, 0, thetadot(t))."""
        pi = 4.0 * math.atan(1.0)
        tdiff = time - self._t0
        arg = (12.0 * tdiff / self._DUR) - 6.0
        coeff = 12.0 / (self._DUR * math.sqrt(pi))
        nerf = math.erf(arg)
        sig = (1.0 / 2) * self._AMP * (pi / 180) * (1.0 + nerf)
        sigdot = self._AMP * (pi / 180) * coeff * math.exp(-arg * arg)
        theta = sig
        thetadot = sigdot
        return TangentSE2(0.0, 0.0, theta, 0.0, 0.0, thetadot)

    def clone(self) -> "SigmoidalStep":
        return SigmoidalStep(self._AMP, self._DUR, self._t0)
