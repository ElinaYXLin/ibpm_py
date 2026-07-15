# fixed_velocity.py
#
# Python port of src/FixedVelocity.h
#
# Subclass of Motion, for a constant (linear + angular) velocity.

from __future__ import annotations

from .motion import Motion
from .tangent_se2 import TangentSE2


class FixedVelocity(Motion):
    """A Motion with constant translational/angular velocity, integrated
    forward in time to produce position."""

    def __init__(self, xdot: float, ydot: float, thetadot: float) -> None:
        self._xdot = xdot
        self._ydot = ydot
        self._thetadot = thetadot
        # NOTE(port): C++ marks these `mutable` so they can be updated
        # inside the `const` method getTransformation(). Python has no
        # const-correctness to route around; they are ordinary instance
        # attributes, mutated directly.
        self._intX = 0.0
        self._intY = 0.0
        self._intTheta = 0.0
        self._oldtime = 0.0

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation for constant velocity, integrated in time
        since the previous call: (x(t), y(t), theta(t), xdot, ydot,
        thetadot)."""
        # NOTE(port): this method has a genuine, unavoidable sequential
        # dependency on the *history of calls* (it integrates velocity
        # using the wall-clock time delta since the previous call), not
        # just on `time` -- there is nothing to numpy-vectorize here, since
        # each call depends on mutated state left by the previous call.
        dtt = time - self._oldtime
        if dtt > 1:
            dtt = 0.0
        self._intX = self._intX + self._xdot * dtt
        self._intY = self._intY + self._ydot * dtt
        self._intTheta = self._intTheta + self._thetadot * dtt
        self._oldtime = time
        return TangentSE2(self._intX, self._intY, self._intTheta, self._xdot, self._ydot, self._thetadot)

    def clone(self) -> "FixedVelocity":
        # NOTE(port): C++ `clone()` only copies the constructor
        # arguments (_xdot, _ydot, _thetadot), *not* the mutable
        # integration state (_intX, _intY, _intTheta, _oldtime) -- the
        # clone restarts integration from zero. Reproduced faithfully here
        # rather than "improving" it to also carry over the integration
        # state.
        return FixedVelocity(self._xdot, self._ydot, self._thetadot)
