# eldredge_combined2.py
#
# Python port of src/EldredgeCombined2.h
#
# Subclass of Motion, combining an Eldredge-2-style (running-integral)
# plunge maneuver with an Eldredge-style pitch maneuver.

from __future__ import annotations

import math

from .motion import Motion
from .tangent_se2 import TangentSE2


class EldredgeCombined2(Motion):
    """Combined plunge (Eldredge2-style, via running integral) and pitch
    (Eldredge-style) maneuver. The `a*` parameters describe the plunge
    profile, the `b*` parameters describe the pitch profile."""

    def __init__(
        self,
        AMPa: float,
        a: float,
        a1: float,
        a2: float,
        a3: float,
        a4: float,
        AMPb: float,
        b: float,
        b1: float,
        b2: float,
        b3: float,
        b4: float,
    ) -> None:
        # a variables are for plunging, b variables are for pitching
        self._AMPa = AMPa
        self._a = a
        self._a1 = a1
        self._a2 = a2
        self._a3 = a3
        self._a4 = a4
        self._AMPb = AMPb
        self._b = b
        self._b1 = b1
        self._b2 = b2
        self._b3 = b3
        self._b4 = b4
        tta = (a2 + a3) / 2.0
        ttb = (b2 + b3) / 2.0
        self._intG = 0.0
        self._oldtime = 0.0
        self._maxGa = math.log(
            (math.cosh(a * (tta - a1)) * math.cosh(a * (tta - a4)))
            / (math.cosh(a * (tta - a2)) * math.cosh(a * (tta - a3)))
        )
        # NOTE(port): the C++ source uses `_a` (the plunge rate parameter),
        # not `_b`, inside the second cosh factor of _maxGb -- i.e.
        # `cosh(_a*(ttb-_b4))` rather than `cosh(_b*(ttb-_b4))`. This looks
        # like a copy-paste bug in the original (mixing plunge rate `a`
        # into the pitch-profile normalization `_maxGb`), but per the
        # "faithful port, don't fix bugs" requirement it is reproduced
        # exactly rather than "corrected" to use `_b`.
        self._maxGb = math.log(
            (math.cosh(b * (ttb - b1)) * math.cosh(a * (ttb - b4)))
            / (math.cosh(b * (ttb - b2)) * math.cosh(b * (ttb - b3)))
        )

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation: (0, h(t), theta(t), 0, hdot(t),
        thetadot(t)), where h (plunge) is obtained via a running integral
        across successive calls -- see the NOTE(port) below.

        NOTE(port): like Eldredge2.getTransformation, the plunge component
        has a genuine sequential dependency on the history of calls (it
        accumulates `_intG` using the wall-clock time delta since the
        previous call); the pitch component (Gb, dGdtb) has no such
        dependency and is a pure function of `time`, so only the plunge
        integration is inherently non-vectorizable here.
        """
        pi = 4.0 * math.atan(1.0)
        AMPrad = self._AMPb * (pi / 180)
        a = self._a
        b = self._b
        a1d = time - self._a1
        a2d = time - self._a2
        a3d = time - self._a3
        a4d = time - self._a4
        b1d = time - self._b1
        b2d = time - self._b2
        b3d = time - self._b3
        b4d = time - self._b4
        Garga = (math.cosh(a * a1d) * math.cosh(a * a4d)) / (math.cosh(a * a2d) * math.cosh(a * a3d))
        Gargb = (math.cosh(b * b1d) * math.cosh(b * b4d)) / (math.cosh(b * b2d) * math.cosh(b * b3d))
        Ga = math.log(Garga)
        Gb = math.log(Gargb)
        dGdtb = b * (math.tanh(b * b1d) - math.tanh(b * b2d) - math.tanh(b * b3d) + math.tanh(b * b4d))

        # The following is a hack to compute integral of G so that I can
        # step vertical velocity...
        dtt = time - self._oldtime
        self._intG = self._intG + Ga * dtt
        self._oldtime = time

        return TangentSE2(
            0.0,
            self._intG * self._AMPa / self._maxGa,
            Gb * AMPrad / self._maxGb,
            0.0,
            Ga * self._AMPa / self._maxGa,
            dGdtb * AMPrad / self._maxGb,
        )

    def clone(self) -> "EldredgeCombined2":
        # NOTE(port): as in Eldredge2.clone(), the mutable integration
        # state (_intG, _oldtime) is intentionally not carried over.
        return EldredgeCombined2(
            self._AMPa,
            self._a,
            self._a1,
            self._a2,
            self._a3,
            self._a4,
            self._AMPb,
            self._b,
            self._b1,
            self._b2,
            self._b3,
            self._b4,
        )
