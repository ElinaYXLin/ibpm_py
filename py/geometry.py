# geometry.py
#
# Python port of src/Geometry.h / src/Geometry.cc
#
# Create, load, and save geometries composed of RigidBody objects.

from __future__ import annotations

from typing import IO, List, Optional, Union

from .boundary_vector import BoundaryVector
from .direction import Direction
from .motion import Motion
from .rigid_body import RigidBody, _parse_error
from .utils import EatWhitespace, MakeLowercase


class Geometry:
    """Create, load, and save geometries composed of RigidBody objects."""

    def __init__(self, filename: Optional[str] = None) -> None:
        # NOTE(port): collapses C++'s Geometry() and Geometry(string
        # filename) constructors into one, dispatching on whether
        # `filename` is given.
        self._bodies: List[RigidBody] = []
        self._numPoints: int = 0
        self._isStationary: bool = True
        if filename is not None:
            self.load(filename)

    def addBody(self, body: RigidBody) -> None:
        """Append the given RigidBody to the list of bodies in the current
        geometry. Makes a copy of it internally."""
        self._bodies.append(RigidBody(body))
        self._numPoints += body.getNumPoints()
        self._isStationary = self._isStationary and body.isStationary()

    def getNumPoints(self) -> int:
        """Return number of boundary points."""
        return self._numPoints

    def getNumBodies(self) -> int:
        """Return number of bodies."""
        return len(self._bodies)

    def transferMotion(self) -> Optional[Motion]:
        """Return motion from one of the bodies, and clear that motion.

        Useful for unsteady baseflow, since we want to initialize
        baseFlow's motion from the motion in a geometry. We need to get
        that motion, and then remove it from the rigid body object it came
        from. We go through the bodies, and take motion from the first one
        that is moving.
        """
        # NOTE(port): faithfully reproduces the C++ implementation,
        # including its apparent bug: it only ever inspects/clears the
        # *first* body (`body = _bodies.begin()`, used once, never
        # advanced before the clear), regardless of which body is actually
        # moving, and calls `tempmotion->clone()` on the return line even
        # when `tempmotion` was never set (i.e. is NULL) if the first body
        # was stationary -- in C++, that's technically UB (calling a
        # method through a NULL pointer), but for this class hierarchy
        # `clone()` doesn't dereference `this` state that would crash in
        # practice on the compilers this project targets. We reproduce the
        # given intent by returning `None` in that case, since dereferencing
        # None similarly cannot be "faithfully" carried over into Python
        # (there is no analogous UB to exploit) -- this is a judgment call.
        body = self._bodies[0]
        tempmotion: Optional[Motion] = None
        if not body.isStationary():
            tempmotion = body.getMotion()
            body.clearMotion()  # removes the motion, resets isStationary to True
        # recompute whether or not remaining RigidBody objects are stationary
        self._isStationary = True
        for b in self._bodies:
            if not b.isStationary():
                self._isStationary = False
        return tempmotion.clone() if tempmotion is not None else None

    def transferCenter(self) -> "tuple[float, float]":
        """Return (x,y) with the center of rotation of the first
        RigidBody."""
        # NOTE(port): C++ signature is `transferCenter(double &x, double
        # &y)` (out-params); Python returns a tuple instead (see grid.py
        # for the same pattern used throughout this port).
        return self._bodies[0].getCenter()

    def getPoints(self) -> BoundaryVector:
        """Return the boundary points in the geometry."""
        coords = BoundaryVector(self._numPoints)
        ind = 0
        for body in self._bodies:
            bodyCoords = body.getPoints()
            n = bodyCoords.getNumPoints()
            # NOTE(port): the C++ loop copies element-by-element
            # (coords(X,ind)=bodyCoords(X,bodyInd); ...; ++ind). Vectorized
            # here with numpy slice assignment per body instead of a
            # manual per-point Python loop.
            coords._data[Direction.X * self._numPoints + ind: Direction.X * self._numPoints + ind + n] = (
                bodyCoords._data[bodyCoords.begin(Direction.X):bodyCoords.end(Direction.X)]
            )
            coords._data[Direction.Y * self._numPoints + ind: Direction.Y * self._numPoints + ind + n] = (
                bodyCoords._data[bodyCoords.begin(Direction.Y):bodyCoords.end(Direction.Y)]
            )
            ind += n
        return coords

    def getVelocities(self) -> BoundaryVector:
        """Return the velocities of the boundary points."""
        velocities = BoundaryVector(self._numPoints)
        ind = 0
        for body in self._bodies:
            bodyVel = body.getVelocities()
            n = bodyVel.getNumPoints()
            velocities._data[Direction.X * self._numPoints + ind: Direction.X * self._numPoints + ind + n] = (
                bodyVel._data[bodyVel.begin(Direction.X):bodyVel.end(Direction.X)]
            )
            velocities._data[Direction.Y * self._numPoints + ind: Direction.Y * self._numPoints + ind + n] = (
                bodyVel._data[bodyVel.begin(Direction.Y):bodyVel.end(Direction.Y)]
            )
            ind += n
        return velocities

    def isStationary(self) -> bool:
        """Return true if the body is not moving; false otherwise."""
        return self._isStationary

    def moveBodies(self, time: float) -> None:
        """Move the boundary points and update their velocities."""
        for body in self._bodies:
            body.moveBody(time)

    def load(self, source: Union[str, IO[str]]) -> bool:
        """Load a geometry from the specified input stream or filename.

        Returns False if invalid input was encountered.
        """
        # NOTE(port): C++ overloads `load` for `istream&` and for
        # `string filename` (the latter opening the file and delegating to
        # the former). Collapsed into one method dispatching on whether
        # `source` is a str (filename) or an open text stream, matching the
        # pattern used for other overloaded constructors/methods elsewhere
        # in this port.
        if isinstance(source, str):
            try:
                with open(source) as f:
                    return self.load(f)
            except OSError:
                print(f"Error: could not open {source} for input.")
                return False

        in_stream = source
        error_found = False
        for buf in in_stream:
            buf = buf.rstrip("\n")
            tokens = buf.split()
            if not tokens:
                continue
            cmd = MakeLowercase(tokens[0])
            if cmd[0] == "#":
                continue
            elif cmd == "body":
                idx = buf.find(tokens[0])
                name = EatWhitespace(buf[idx + len(tokens[0]):])
                body = RigidBody()
                body.setName(name)
                if body.load(in_stream):
                    self.addBody(body)
                else:
                    error_found = True
            elif cmd == "end":
                break
            else:
                _parse_error(buf)
                error_found = True
        return not error_found
