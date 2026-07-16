# pitch_plunge.py
#
# Python port of src/PitchPlunge.h
#
# Subclass of Motion, for a pitching and plunging body.

from __future__ import annotations

import math

from .motion import Motion
from .tangent_se2 import TangentSE2


class PitchPlunge(Motion):
    """Sinusoidal pitching and plunging, centered about the origin:
        y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
        theta(t)  = A_theta  sin( 2*pi * f_theta t + phi_theta)
    """

    def __init__(
        self,
        pitchAmplitude: float,
        pitchFrequency: float,
        pitchPhase: float,
        plungeAmplitude: float,
        plungeFrequency: float,
        plungePhase: float,
    ) -> None:
        self._pitchAmp = pitchAmplitude
        self._pitchFreq = pitchFrequency
        self._pitchPhase = pitchPhase
        self._plungeAmp = plungeAmplitude
        self._plungeFreq = plungeFrequency
        self._plungePhase = plungePhase

        twopi = 8.0 * math.atan(1.0)
        self._pitchFreq *= twopi
        self._plungeFreq *= twopi

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation for sinusoidal pitch/plunge:
        (0, y(t), theta(t), 0, ydot(t), thetadot(t))."""
        y = self._plungeAmp * math.sin(self._plungeFreq * time + self._plungePhase)
        ydot = self._plungeAmp * self._plungeFreq * math.cos(self._plungeFreq * time + self._plungePhase)
        theta = self._pitchAmp * math.sin(self._pitchFreq * time + self._pitchPhase)
        thetadot = self._pitchAmp * self._pitchFreq * math.cos(self._pitchFreq * time + self._pitchPhase)
        return TangentSE2(0.0, y, theta, 0.0, ydot, thetadot)

    def clone(self) -> "PitchPlunge":
        # NOTE(port): C++ divides _pitchFreq/_plungeFreq by 2*pi before
        # passing them back into the constructor, undoing the *=twopi done
        # in __init__ above -- otherwise a clone-of-a-clone would multiply
        # by twopi again. Reproduced exactly (not "simplified" by storing
        # the un-scaled frequency separately), to keep clone() a faithful
        # inverse of __init__ as in the original.
        twopi = 8.0 * math.atan(1.0)
        return PitchPlunge(
            self._pitchAmp,
            self._pitchFreq / twopi,
            self._pitchPhase,
            self._plungeAmp,
            self._plungeFreq / twopi,
            self._plungePhase,
        )
