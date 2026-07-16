# rigid_body.py
#
# Python port of src/RigidBody.h / src/RigidBody.cc
#
# Specify coordinates and center of a rigid body, and its motion.
#
# The points on the body are stored with respect to a reference
# configuration in _refPoints. The current locations of the points on the
# body, defined by the associated Motion, are contained in _currentPoints,
# and are updated whenever moveBody() is called.

from __future__ import annotations

import math
from typing import IO, List, Optional, Union

import numpy as np

from .boundary_vector import BoundaryVector
from .direction import Direction
from .eldredge1 import Eldredge1
from .eldredge2 import Eldredge2
from .eldredge_combined2 import EldredgeCombined2
from .eldredge_maneuver import EldredgeManeuver
from .fixed_position import FixedPosition
from .fixed_velocity import FixedVelocity
from .lag_step1 import LagStep1
from .lag_step2 import LagStep2
from .motion import Motion
from .motion_file import MotionFile
from .motion_file_periodic import MotionFilePeriodic
from .pitch_plunge import PitchPlunge
from .sigmoidal_step import SigmoidalStep
from .utils import EatWhitespace, MakeLowercase


class Point:
    """A point in 2d. (C++: `struct Point` in RigidBody.h.)"""

    def __init__(self, x_in: float, y_in: float) -> None:
        self.x = x_in
        self.y = y_in


def _parse_error(buf: str) -> None:
    print("WARNING: could not parse the following line:")
    print(buf)


class RigidBody:
    """Coordinates and center of a rigid body, and its motion."""

    def __init__(self, body: Optional["RigidBody"] = None) -> None:
        # NOTE(port): collapses the C++ default constructor and copy
        # constructor RigidBody(const RigidBody&) into one, dispatching on
        # whether `body` is given (see grid.py/field.py for the same
        # pattern used elsewhere in this port).
        if body is None:
            self._name: str = ""
            self._xCenter: float = 0.0
            self._yCenter: float = 0.0
            self._isStationary: bool = True
            self._refPoints: List[Point] = []
            self._currentPoints: List[Point] = []
            self._currentVelocities: List[Point] = []
            self._motion: Optional[Motion] = None
            return
        self._xCenter = body._xCenter
        self._yCenter = body._yCenter
        self._isStationary = body._isStationary
        self._name = body._name
        self._refPoints = [Point(p.x, p.y) for p in body._refPoints]
        self._currentPoints = [Point(p.x, p.y) for p in body._currentPoints]
        self._currentVelocities = [Point(p.x, p.y) for p in body._currentVelocities]
        self._motion = body._motion.clone() if body._motion is not None else None

    def addPoint(self, x: float, y: float) -> None:
        """Add the specified point to the list of points on the body's
        boundary."""
        self._refPoints.append(Point(x, y))
        self._currentPoints.append(Point(x, y))
        self._currentVelocities.append(Point(0.0, 0.0))

    def addCircle(self, xc: float, yc: float, radius: float, dx: float) -> None:
        """Add a circle with center (xc, yc) and the given radius with the
        specified (approximate) distance between points."""
        dTheta = dx / radius
        twopi = 8.0 * math.atan(1.0)
        # To round a value x, take floor( x + 0.5 )
        numPoints = int(math.floor(twopi / dTheta + 1 + 0.5))
        self.addCircle_n(xc, yc, radius, numPoints)

    def addCircle_n(self, xc: float, yc: float, radius: float, numPoints: int) -> None:
        """Add a circle with center (xc, yc) and the given radius with the
        specified number of points."""
        twopi = 8.0 * math.atan(1.0)
        dTheta = twopi / numPoints
        # NOTE(port): the coordinate math (which has no cross-iteration
        # dependency) is vectorized with numpy; the loop that calls
        # addPoint() is kept because addPoint() has the side effect of
        # appending to three separate lists (_refPoints, _currentPoints,
        # _currentVelocities), matching the original algorithm exactly.
        theta = np.arange(numPoints) * dTheta
        x = xc + radius * np.cos(theta)
        y = yc + radius * np.sin(theta)
        for i in range(numPoints):
            self.addPoint(float(x[i]), float(y[i]))

    def addLine(self, x1: float, y1: float, x2: float, y2: float, dx: float) -> None:
        """Add a line connecting (x1,y1) and (x2,y2) with the specified
        (approximate) distance between points."""
        length = math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
        # To round a value x, take floor( x + 0.5 )
        numPoints = int(math.floor(length / dx + 1 + 0.5))
        self.addLine_n(x1, y1, x2, y2, numPoints)

    def addLine_n(self, x1: float, y1: float, x2: float, y2: float, numPoints: int) -> None:
        """Add a line connecting (x1,y1) and (x2,y2) with the specified
        number of points."""
        deltaX = (x2 - x1) / (numPoints - 1)
        deltaY = (y2 - y1) / (numPoints - 1)
        i = np.arange(numPoints)
        x = x1 + i * deltaX
        y = y1 + i * deltaY
        for k in range(numPoints):
            self.addPoint(float(x[k]), float(y[k]))

    def addLine_aoa(
        self,
        l: float,
        xC: float,
        yC: float,
        alpha: float,
        numPoints: int,
    ) -> None:
        """Add a line with length l, centered at (0,0) with AoA alpha (in
        degrees) and the specified number of points."""
        pi = 3.141592653589793238462643383279502884197169399375
        a = alpha / 180 * pi
        cosa = math.cos(a)
        sina = math.sin(a)
        delta = l / (numPoints - 1)
        i = np.arange(numPoints)
        x0 = i * delta
        y0 = np.zeros(numPoints)  # (x0,y0) point before rotation
        x0r = x0 - xC
        y0r = y0 - yC  # (x0r,y0r) referenced to center (xC,yC)
        x1r = x0r * cosa - y0r * sina
        y1r = x0r * sina + y0r * cosa  # after rotation, referenced to center (xC,yC)
        x1 = x1r + xC
        y1 = y1r + yC
        for k in range(numPoints):
            self.addPoint(float(x1[k]), float(y1[k]))

    def load(self, in_stream: IO[str]) -> bool:
        """Load a list of commands from the specified input stream.

        See the module-level docstring / RigidBody.h for the input format.
        Whitespace at the beginning of the line is ignored. Returns False
        if invalid input was encountered.
        """
        # NOTE(port): C++ tokenizes each line with `istringstream::
        # operator>>`, which splits on arbitrary whitespace and leaves the
        # stream in a fail state on a type mismatch (checked via
        # RB_CHECK_FOR_ERRORS). Python's `str.split()` plus explicit
        # `float()`/`int()` conversions with a try/except reproduces the
        # same "split on whitespace, flag malformed tokens" behavior.
        error_found = False
        for buf in in_stream:
            buf = buf.rstrip("\n")
            tokens = buf.split()
            if not tokens:
                continue
            cmd = MakeLowercase(tokens[0])
            if cmd[0] == "#":
                continue
            elif cmd == "center":
                try:
                    x, y = float(tokens[1]), float(tokens[2])
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
                self.setCenter(x, y)
            elif cmd == "circle":
                try:
                    xc, yc, radius, dx = (float(t) for t in tokens[1:5])
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
                self.addCircle(xc, yc, radius, dx)
            elif cmd == "circle_n":
                try:
                    xc, yc, radius = (float(t) for t in tokens[1:4])
                    numPoints = int(tokens[4])
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
                self.addCircle_n(xc, yc, radius, numPoints)
            elif cmd == "end":
                break
            elif cmd == "line":
                try:
                    x0, y0, x1, y1, dx = (float(t) for t in tokens[1:6])
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
                self.addLine(x0, y0, x1, y1, dx)
            elif cmd == "line_n":
                try:
                    x0, y0, x1, y1 = (float(t) for t in tokens[1:5])
                    numPoints = int(tokens[5])
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
                self.addLine_n(x0, y0, x1, y1, numPoints)
            elif cmd == "line_aoa":
                try:
                    l, xC, yC, aoa = (float(t) for t in tokens[1:5])
                    numPoints = int(tokens[5])
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
                self.addLine_aoa(l, xC, yC, aoa, numPoints)
                self.setCenter(xC, yC)
            elif cmd == "motion":
                if len(tokens) < 2:
                    _parse_error(buf)
                    error_found = True
                    continue
                motionType = MakeLowercase(tokens[1])
                args = tokens[2:]
                try:
                    if motionType == "fixed":
                        x, y, theta = (float(t) for t in args[0:3])
                        self.setMotion(FixedPosition(x, y, theta))
                    elif motionType == "fixedvel":
                        xdot, ydot, thetadot = (float(t) for t in args[0:3])
                        self.setMotion(FixedVelocity(xdot, ydot, thetadot))
                    elif motionType == "pitchplunge":
                        amp1, freq1, phase1, amp2, freq2, phase2 = (float(t) for t in args[0:6])
                        self.setMotion(PitchPlunge(amp1, freq1, phase1, amp2, freq2, phase2))
                    elif motionType == "sigmoidalstep":
                        AMP, DUR, startTime = (float(t) for t in args[0:3])
                        self.setMotion(SigmoidalStep(AMP, DUR, startTime))
                    elif motionType == "lagstep1":
                        AMP, PW, TAU, T0 = (float(t) for t in args[0:4])
                        self.setMotion(LagStep1(AMP, PW, TAU, T0))
                    elif motionType == "lagstep2":
                        AMP, PW, TAU, T0 = (float(t) for t in args[0:4])
                        self.setMotion(LagStep2(AMP, PW, TAU, T0))
                    elif motionType == "eldredge":
                        AMP, a, t1, t2, t3, t4 = (float(t) for t in args[0:6])
                        self.setMotion(EldredgeManeuver(AMP, a, t1, t2, t3, t4))
                    elif motionType == "eldredgecombined2":
                        AMPa, a, a1, a2, a3, a4, AMPb, b, b1, b2, b3, b4 = (float(t) for t in args[0:12])
                        self.setMotion(EldredgeCombined2(AMPa, a, a1, a2, a3, a4, AMPb, b, b1, b2, b3, b4))
                    elif motionType == "eldredge1":
                        AMP, a, t1, t2, t3, t4 = (float(t) for t in args[0:6])
                        self.setMotion(Eldredge1(AMP, a, t1, t2, t3, t4))
                    elif motionType == "eldredge2":
                        AMP, a, t1, t2, t3, t4 = (float(t) for t in args[0:6])
                        self.setMotion(Eldredge2(AMP, a, t1, t2, t3, t4))
                    elif motionType == "motionfile":
                        filename = args[0]
                        self.setMotion(MotionFile(filename))
                    elif motionType == "motionfileperiodic":
                        filename = args[0]
                        period = float(args[1])
                        self.setMotion(MotionFilePeriodic(filename, period))
                    # NOTE(port): C++ has no final `else` here either -- an
                    # unrecognized motionType is silently ignored, matching
                    # that (arguably buggy) original behavior faithfully.
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
            elif cmd == "name":
                idx = buf.find(tokens[0])
                name = EatWhitespace(buf[idx + len(tokens[0]):])
                self.setName(name)
            elif cmd == "point":
                try:
                    x, y = float(tokens[1]), float(tokens[2])
                except (IndexError, ValueError):
                    _parse_error(buf)
                    error_found = True
                    continue
                self.addPoint(x, y)
            elif cmd == "raw":
                try:
                    filename = tokens[1]
                except IndexError:
                    _parse_error(buf)
                    error_found = True
                    continue
                self.loadRaw(filename)
            else:
                _parse_error(buf)
        return not error_found

    def loadRaw(self, fname: str) -> bool:
        """Load a list of points, in ASCII format, with one point per line.
        Assumes the center is (0,0). Returns False if invalid input was
        encountered."""
        # NOTE(port): C++ uses `ifstream >> n` then `ifstream >> x >> y`
        # repeatedly, which tokenizes on arbitrary whitespace (not
        # necessarily one point per line, despite the docstring). Python
        # reproduces that by reading the whole file and splitting on
        # whitespace, rather than iterating line-by-line.
        try:
            with open(fname) as f:
                tokens = f.read().split()
        except OSError:
            return False
        if not tokens:
            return False
        try:
            n = int(tokens[0])
        except ValueError:
            return False
        if len(tokens) < 1 + 2 * n:
            return False
        try:
            coords = [float(t) for t in tokens[1:1 + 2 * n]]
        except ValueError:
            return False
        for i in range(n):
            self.addPoint(coords[2 * i], coords[2 * i + 1])
        return True

    def saveRaw(self, out: IO[str]) -> None:
        """Save a list of points to the specified output stream."""
        # NOTE(port): C++ uses `setw(10)` with the stream's default float
        # formatting (roughly "%g" with 6 significant digits). Python has
        # no exact equivalent of ostream's default float formatting;
        # "{:10.6g}" is used here as the closest standard match (6
        # significant digits, right-justified in a 10-character field).
        n = self.getNumPoints()
        out.write(str(n))
        for p in self._refPoints:
            out.write(f"\n{p.x:10.6g}{p.y:10.6g}")

    def getNumPoints(self) -> int:
        """Return the number of points on the body's boundary."""
        return len(self._refPoints)

    @staticmethod
    def toBoundaryVector(points: List[Point]) -> BoundaryVector:
        n = len(points)
        BVList = BoundaryVector(n)
        # NOTE(port): the C++ loop assigns element-by-element via
        # BVList(X,i)=...; BVList(Y,i)=.... Vectorized here with a single
        # numpy assignment per component instead of a manual Python loop.
        BVList._data[BVList.begin(Direction.X):BVList.end(Direction.X)] = [p.x for p in points]
        BVList._data[BVList.begin(Direction.Y):BVList.end(Direction.Y)] = [p.y for p in points]
        return BVList

    def getPoints(self) -> BoundaryVector:
        """Return the list of coordinates for each point on the body."""
        return RigidBody.toBoundaryVector(self._currentPoints)

    def getVelocities(self) -> BoundaryVector:
        """Return the list of velocities at each point on the body."""
        return RigidBody.toBoundaryVector(self._currentVelocities)

    def isStationary(self) -> bool:
        """Return true if the body is not moving in time."""
        return self._isStationary

    def setMotion(self, motion: Motion) -> None:
        """Set the evolution of the current body (which may be stationary
        or not)."""
        # make a local copy of the new motion
        self._motion = motion.clone()
        self._isStationary = motion.isStationary()

    def clearMotion(self) -> None:
        """Set the motion of the current body to NULL."""
        self._motion = None
        self._isStationary = True

    def getMotion(self) -> Optional[Motion]:
        """Get a clone of the motion of the current body."""
        return self._motion.clone() if self._motion is not None else None

    def setCenter(self, x: float, y: float) -> None:
        """Set the center of the body, about which rotations are defined."""
        self._xCenter = x
        self._yCenter = y

    def getCenter(self) -> "tuple[float, float]":
        """Get the center of the body, about which rotations are defined."""
        # NOTE(port): C++ returns via reference out-params (double& x,
        # double& y); Python returns a tuple (x, y) instead (see grid.py
        # for the same pattern).
        return self._xCenter, self._yCenter

    def moveBody(self, time: float) -> None:
        """Update the position of the body, based on the Motion."""
        if self._motion is None:
            return
        g = self._motion.getTransformation(time)

        # NOTE(port): C++ loops over _refPoints, calling
        # g.mapPosition/g.mapVelocity (each an O(1) affine formula) once
        # per point. Rather than replicate that with a Python loop calling
        # into TangentSE2 per point (a manual loop where numpy
        # vectorization is equivalent, since mapPosition/mapVelocity are
        # simple affine functions of (a,b)), the same formulas are applied
        # here directly to numpy arrays of all reference points at once.
        gx, gy, gtheta = g.getPosition()
        gxdot, gydot, gthetadot = g.getVelocity()
        cost = math.cos(gtheta)
        sint = math.sin(gtheta)
        xs = np.array([p.x for p in self._refPoints]) - self._xCenter
        ys = np.array([p.y for p in self._refPoints]) - self._yCenter

        xnew = gx + xs * cost - ys * sint
        ynew = gy + xs * sint + ys * cost

        u = -xs * sint * gthetadot - ys * cost * gthetadot + gxdot
        v = xs * cost * gthetadot - ys * sint * gthetadot + gydot

        self._currentPoints = [
            Point(float(xnew[k]) + self._xCenter, float(ynew[k]) + self._yCenter)
            for k in range(len(self._refPoints))
        ]
        self._currentVelocities = [Point(float(u[k]), float(v[k])) for k in range(len(self._refPoints))]

    def setName(self, name: str) -> None:
        """Set the name of the body."""
        self._name = name

    def getName(self) -> str:
        """Return the name of the body."""
        return self._name
