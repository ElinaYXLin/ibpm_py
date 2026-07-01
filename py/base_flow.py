# base_flow.py
#
# Python port of src/BaseFlow.h / src/BaseFlow.cc
#
# Structure for a BaseFlow (flux) that moves with time.

from __future__ import annotations

import math
from typing import Optional, Union

from .flux import Flux
from .grid import Grid
from .motion import Motion
from .tangent_se2 import TangentSE2

FloatOrFlux = Union[Flux, float]


class BaseFlow:
    """Structure for a BaseFlow (flux) that moves with time."""

    def __init__(
        self,
        grid: Optional[Grid] = None,
        mag: Optional[float] = None,
        alpha: Optional[float] = None,
        motion: Optional[Motion] = None,
    ) -> None:
        # NOTE(port): C++ has four constructors:
        #   BaseFlow()                                    -- no allocation
        #   BaseFlow(const Grid&)                         -- no motion/baseflow
        #   BaseFlow(const Grid&, double mag, double alpha)      -- no motion
        #   BaseFlow(const Grid&, const Motion&)          -- no constant baseflow
        #   BaseFlow(const Grid&, double mag, double alpha, const Motion&)
        # Note the C++ header also declares a `BaseFlow(const Grid&, const
        # Motion&)` overload (no mag/alpha), but no such constructor is
        # defined in BaseFlow.cc -- only the four bodies above exist there.
        # Collapsed here into a single constructor dispatching on which
        # optional arguments are supplied; `grid=None` reproduces the
        # default constructor, and passing `motion` without `mag`/`alpha`
        # would hit the header-only, undefined-in-.cc overload (there is
        # nothing to faithfully port for that case, so it is simply not
        # reachable here -- callers must supply mag and alpha whenever
        # motion is supplied, matching what BaseFlow.cc actually defines).
        self._xCenter: float = 0.0
        self._yCenter: float = 0.0
        self._isStationary: bool = True
        self._motion: Optional[Motion] = None
        self._time: float = 0.0
        self._mag: float = 0.0
        self._magBF: float = 0.0
        self._alpha: float = 0.0
        self._alphaBF: float = 0.0
        self._gamma: float = 0.0
        # NOTE(port): C++ `_q` is a `Flux` data member, default-constructed
        # (no grid, no data) whenever a BaseFlow constructor body does not
        # explicitly assign to it -- mirrored here by always constructing a
        # default (empty) Flux() up front, then overwriting it below for
        # every constructor variant that calls resize()/UniformFlow.
        self._q: Flux = Flux()

        if grid is None:
            return

        if motion is not None:
            assert mag is not None and alpha is not None
            # NOTE(port) -- judgment call: this faithfully reproduces an
            # apparent bug in BaseFlow.cc's four-argument constructor.
            # There, `_motion = motion.clone();` is called, then
            # `_isStationary = motion.isStationary();`, and then
            # immediately `_motion = NULL;` -- discarding the clone that was
            # just made (a leaked allocation in C++; simply
            # garbage-collected here) and leaving `_motion` NULL exactly as
            # in the default constructor. `_isStationary` is still set
            # correctly from `motion` before the clone is discarded. This
            # is reproduced verbatim below rather than "fixed" (e.g. by
            # keeping the clone), per the instruction to do a faithful
            # port only.
            self._motion = motion.clone()
            self._isStationary = motion.isStationary()
            self._motion = None
            self._time = 0.0
            self._mag = mag
            self._magBF = mag
            self._alpha = alpha
            self._alphaBF = alpha
            self._gamma = -1.0 * self._alphaBF
            self.resize(grid)
            self._q = Flux.UniformFlow(grid, self._magBF, self._alphaBF)
            return

        if mag is not None and alpha is not None:
            self._time = 0.0
            self._mag = mag
            self._magBF = mag
            self._alpha = alpha
            self._alphaBF = alpha
            self._gamma = -1.0 * self._alphaBF
            self.resize(grid)
            self._q = Flux.UniformFlow(grid, self._magBF, self._alphaBF)
            return

        # BaseFlow(const Grid& grid)
        self._time = 0.0
        self.resize(grid)
        self._q = Flux.UniformFlow(grid, self._magBF, self._alphaBF)

    def resize(self, grid: Grid) -> None:
        """Allocate memory, with the specified Grid and number of boundary
        points."""
        self._q.resize(grid)

    def isStationary(self) -> bool:
        """Return true if the body is not moving in time."""
        return self._isStationary

    def setMotion(self, motion: Motion) -> None:
        """Set the evolution of the current BaseFlow (which may be
        stationary or not)."""
        # NOTE(port): C++ explicitly `delete`s the old `_motion` pointer
        # before overwriting it (manual memory management); Python relies
        # on garbage collection, so the old Motion is simply dropped by
        # reassignment.
        self._motion = motion.clone()
        self._isStationary = motion.isStationary()

    def setCenter(self, x: float, y: float) -> None:
        """Set the center of the domain, about which rotations are
        defined."""
        self._xCenter = x
        self._yCenter = y

    def getCenter(self) -> "tuple[float, float]":
        """Get the center of the domain, about which rotations are
        defined."""
        # NOTE(port): C++ signature is `getCenter(double& x, double& y)`
        # (out-params); Python returns a tuple instead, matching the same
        # convention used throughout this port (see grid.py, rigid_body.py).
        return self._xCenter, self._yCenter

    def getAlpha(self) -> float:
        """Get the angle of attack."""
        return self._alpha

    def getMag(self) -> float:
        """Get the magnitude of base flow."""
        return self._mag

    def setAlphaMag(self, time: float) -> None:
        """Determine the magnitude and angle of the base flow, including
        rigid body motion."""
        g = self._motion.getTransformation(time)
        x, y, theta = g.getPosition()
        xdot, ydot, thetadot = g.getVelocity()
        xdotBF = self._magBF * math.cos(self._alphaBF)
        ydotBF = self._magBF * math.sin(self._alphaBF)
        xdotT = xdot - xdotBF
        ydotT = ydot - ydotBF
        self._gamma = math.atan2(ydotT, -1.0 * xdotT)
        self._alpha = -1.0 * theta - self._gamma
        self._mag = math.sqrt(xdotT * xdotT + ydotT * ydotT)

    def moveFlow(self, time: float) -> None:
        """Update the BaseFlow, based on the Motion."""
        if self._motion is None:
            return
        g = self._motion.getTransformation(time)
        x, y, theta = g.getPosition()
        xdot, ydot, thetadot = g.getVelocity()
        # The flow is decomposed into a base flow of magnitude _mag at an
        # angle _alpha = -theta-_gamma, and a purely rotational component
        # -thetadot centered at the body center of rotation. The formulae
        # are:
        #   _mag = (_magBF*cos(_alphaBF) - xdot, -_magBF*sin(_alphaBF) + ydot)
        #   _alpha = -theta - gamma
        self.setAlphaMag(time)

        xdot = self._mag * math.cos(self._alpha)  # using the xdot,ydot of the baseflow
        ydot = self._mag * math.sin(self._alpha)
        gnew = TangentSE2(x, y, theta, xdot, ydot, -1.0 * thetadot)
        # Update the baseFlow based on this new motion
        self._q.setFlow(gnew, self._xCenter, self._yCenter)

    def setFlux(self, f: FloatOrFlux) -> None:
        """Set the value of the flux _q."""
        # NOTE(port): collapses C++'s `setFlux(Flux& f)` and
        # `setFlux(double f)` overloads (the latter setting every element
        # of _q to a constant, via Flux::operator=(double)) into a single
        # method dispatching on the argument type, matching the pattern
        # used for other overloaded setters/constructors in this port.
        # NOTE(port): C++ `_q = f` invokes `Flux::operator=`, a deep copy of
        # the underlying data (not a pointer/reference assignment). `assign`
        # is used for both branches here (rather than `self._q = f`) to
        # preserve that copy semantics instead of aliasing `f`.
        self._q.assign(f)

    def getFlux(self) -> Flux:
        """Return flux _q."""
        return self._q
