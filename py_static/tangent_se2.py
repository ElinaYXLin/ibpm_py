# tangent_se2.py
#
# Python port of src/TangentSE2.h
#
# An abstraction of the group TSE(2) of transformations in 2D. Keeps track
# of translations and rotations in 2D, and their velocities.

from __future__ import annotations

import math
from typing import Tuple


class TangentSE2:
    """Element of the tangent bundle of SE(2): a base point (x, y, theta)
    together with its velocity (xdot, ydot, thetadot)."""

    def __init__(
        self,
        x: float,
        y: float,
        theta: float,
        xdot: float,
        ydot: float,
        thetadot: float,
    ) -> None:
        self._x = x
        self._y = y
        self._theta = theta
        self._xdot = xdot
        self._ydot = ydot
        self._thetadot = thetadot

    def setPosition(self, x: float, y: float, theta: float) -> None:
        self._x = x
        self._y = y
        self._theta = theta

    def getPosition(self) -> Tuple[float, float, float]:
        # NOTE(port): C++ returns via reference out-params (double& x,
        # double& y, double& theta); Python has no out-params, so this
        # returns a tuple (x, y, theta) instead. Callers must be updated
        # accordingly.
        return self._x, self._y, self._theta

    def setVelocity(self, xdot: float, ydot: float, thetadot: float) -> None:
        self._xdot = xdot
        self._ydot = ydot
        self._thetadot = thetadot

    def getVelocity(self) -> Tuple[float, float, float]:
        # NOTE(port): see getPosition() above regarding out-params -> tuple.
        return self._xdot, self._ydot, self._thetadot

    def mapPosition(self, a: float, b: float) -> Tuple[float, float]:
        """Given the point (a,b), compute the mapped point (a_new, b_new)."""
        cost = math.cos(self._theta)
        sint = math.sin(self._theta)
        a_new = self._x + a * cost - b * sint
        b_new = self._y + a * sint + b * cost
        return a_new, b_new

    def mapVelocity(self, a: float, b: float) -> Tuple[float, float]:
        """Given the point (a,b) (with zero initial velocity), compute the
        mapped velocity (u_new, v_new)."""
        cost = math.cos(self._theta)
        sint = math.sin(self._theta)
        u_new = -a * sint * self._thetadot - b * cost * self._thetadot + self._xdot
        v_new = a * cost * self._thetadot - b * sint * self._thetadot + self._ydot
        return u_new, v_new
