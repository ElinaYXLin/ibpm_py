# eldredge2.py
#
# Python port of src/Eldredge2.h
#
# Subclass of Motion, for the canonical maneuver described by Eldredge,
# applied to y (plunge) via a running integral of G, so that ydot(t) is the
# canonical Eldredge profile and y(t) is (approximately) its integral.

from __future__ import annotations

import math

from .motion import Motion
from .tangent_se2 import TangentSE2


class Eldredge2(Motion):
    """Eldredge's canonical smoothed-ramp maneuver, applied to plunge
    velocity, with plunge position obtained by running (rectangle-rule)
    integration across successive getTransformation() calls."""

    def __init__(self, AMP: float, a: float, t1: float, t2: float, t3: float, t4: float) -> None:
        self._AMP = AMP
        self._a = a
        self._t1 = t1
        self._t2 = t2
        self._t3 = t3
        self._t4 = t4
        tt = (t2 + t3) / 2.0
        # NOTE(port): the original constructor prints two debug lines
        # ("got her" / "not here") bracketing the mutable-state
        # initialization, and has a stray unused local `tt` used only in
        # _maxG below. Reproduced verbatim (including the debug prints)
        # per the "faithful port, no improving the algorithm" requirement,
        # even though they look like leftover debugging output.
        print(" got her")
        self._intG = 0.0
        self._oldtime = 0.0
        print(" not here")
        self._maxG = math.log(
            (math.cosh(a * (tt - t1)) * math.cosh(a * (tt - t4)))
            / (math.cosh(a * (tt - t2)) * math.cosh(a * (tt - t3)))
        )

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation: (0, intG(t)*AMP/maxG, 0, 0,
        G(t)*AMP/maxG, 0), where intG is a running integral of G updated
        across successive calls (see NOTE(port) below).

        NOTE(port): this method has a genuine sequential dependency on the
        history of calls (it accumulates `_intG` using the wall-clock time
        delta since the previous call, exactly like FixedVelocity.
        getTransformation), so there is nothing here to numpy-vectorize.
        """
        a = self._a
        t1d = time - self._t1
        t2d = time - self._t2
        t3d = time - self._t3
        t4d = time - self._t4
        Garg = (math.cosh(a * t1d) * math.cosh(a * t4d)) / (math.cosh(a * t2d) * math.cosh(a * t3d))
        G = math.log(Garg)
        dGdt = a * (math.tanh(a * t1d) - math.tanh(a * t2d) - math.tanh(a * t3d) + math.tanh(a * t4d))

        # The following is a hack to compute integral of G so that I can
        # step vertical velocity...
        dtt = time - self._oldtime
        self._intG = self._intG + G * dtt
        self._oldtime = time
        print(f"dtt = {dtt}")
        print(f"intG = {self._intG}")
        print(f"G = {G * self._AMP / self._maxG}")
        print(f"intG*stuff = {self._intG * self._AMP / self._maxG}")
        return TangentSE2(0.0, self._intG * self._AMP / self._maxG, 0.0, 0.0, G * self._AMP / self._maxG, 0.0)

    def clone(self) -> "Eldredge2":
        # NOTE(port): as in FixedVelocity.clone(), the mutable integration
        # state (_intG, _oldtime) is intentionally not carried over --
        # matches C++, which only forwards the constructor arguments.
        return Eldredge2(self._AMP, self._a, self._t1, self._t2, self._t3, self._t4)
