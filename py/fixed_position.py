# fixed_position.py
#
# Python port of src/FixedPosition.h
#
# Subclass of Motion, for a stationary body: specifies the location of the
# center of the body, and an angle of rotation about the center.
#
# NOTE(port): FixedPosition is the only Motion subclass ported here. It is
# ported because it is the minimal concrete Motion needed to exercise
# RigidBody's motion machinery (setMotion/getMotion/moveBody) and because
# RigidBody::load()'s "motion fixed x y theta" command is documented
# directly in RigidBody.h's own docstring. The C++ codebase has ~10 other
# Motion subclasses (FixedVelocity, PitchPlunge, SigmoidalStep, LagStep1/2,
# EldredgeManeuver/Combined2/1/2, MotionFile, MotionFilePeriodic) that are
# only reachable through RigidBody::load()'s "motion <type> ..." command
# parser, not through any other public entry point of RigidBody/Geometry.
# Per the porting scope ("only port other files when there are
# dependencies"), those are treated as optional plugins of the config-file
# parser rather than hard dependencies of RigidBody/Geometry, and are left
# unported; RigidBody.load() raises NotImplementedError with a clear
# message if one of those motion types is requested (see rigid_body.py).

from __future__ import annotations

from .motion import Motion
from .tangent_se2 import TangentSE2


class FixedPosition(Motion):
    """A Motion corresponding to a fixed position and rotation about the
    center."""

    def __init__(self, x: float, y: float, theta: float) -> None:
        self._x = x
        self._y = y
        self._theta = theta

    def isStationary(self) -> bool:
        return True

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns a fixed transformation for all time: (x,y,theta,0,0,0)."""
        return TangentSE2(self._x, self._y, self._theta, 0.0, 0.0, 0.0)

    def clone(self) -> "FixedPosition":
        return FixedPosition(self._x, self._y, self._theta)
