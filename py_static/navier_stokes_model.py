# navier_stokes_model.py
#
# Python port of src/NavierStokesModel.h / src/NavierStokesModel.cc
#
# Define operators for the Navier-Stokes equations.
#
# Note: in C++ this is documented as "an abstract base class, and may not be
# instantiated directly" -- but NavierStokesModel.h declares no pure virtual
# methods, so nothing actually prevents C++ from instantiating it either;
# this docstring note is carried over as-is (see the class docstring below)
# without adding any Python-side enforcement (e.g. abc.ABC), since there is
# no abstract method to require overriding.

from __future__ import annotations

from typing import Optional

from .base_flow import BaseFlow
from .boundary_vector import BoundaryVector
from .elliptic_solver import PoissonSolver
from .flux import Flux
from .geometry import Geometry
from .grid import Grid
from .regularizer import Regularizer
from .scalar import Scalar
from .state import State
from .vector_operations import Curl


class NavierStokesModel:
    """Define operators for the Navier-Stokes equations.

    Note: This is an abstract base class, and may not be instantiated
    directly.
    """

    def __init__(
        self,
        grid: Grid,
        geometry: Geometry,
        Reynolds: float,
        q_potential: Optional[BaseFlow] = None,
    ) -> None:
        # NOTE(port): collapses C++'s two constructors --
        #   NavierStokesModel(grid, geometry, Reynolds, q_potential)
        #   NavierStokesModel(grid, geometry, Reynolds)
        # -- into one, dispatching on whether `q_potential` is given,
        # matching the pattern used throughout this port.
        self._grid = grid
        self._geometry = geometry
        self._regularizer = Regularizer(grid, geometry)
        self._ReynoldsNumber = Reynolds
        self._poisson = PoissonSolver(grid)
        self._hasBeenInitialized = False

        if q_potential is not None:
            # NOTE(port): C++ `_baseFlow( q_potential )` invokes BaseFlow's
            # implicit (compiler-generated) copy constructor -- see the
            # NOTE(port) in base_flow.py's constructor for the judgment call
            # this entails (a shallow/aliased copy of the `_motion` member).
            self._baseFlow = BaseFlow(q_potential)
        else:
            self._baseFlow = BaseFlow(grid)
            self._baseFlow.setFlux(0.0)

    def init(self) -> None:
        """Perform initial calculations needed to use model."""
        if self._hasBeenInitialized:
            return  # do only once
        # Update regularizer
        self._regularizer.update()
        self._hasBeenInitialized = True

    def isTimeDependent(self) -> bool:
        """Return true if the geometry has moving bodies."""
        flag = False
        if (not self._geometry.isStationary()) or (not self._baseFlow.isStationary()):
            flag = True
        return flag

    def geTimeDependent(self) -> bool:
        """geometry TD?"""
        flag = False
        if not self._geometry.isStationary():
            flag = True
        return flag

    def bfTimeDependent(self) -> bool:
        """baseflow TD?"""
        flag = False
        if not self._baseFlow.isStationary():
            flag = True
        return flag

    def getNumPoints(self) -> int:
        """Return the number of points in the geometry."""
        return self._geometry.getNumPoints()

    def getConstraints(self) -> BoundaryVector:
        """Return the right-hand side b of the constraint equations.

        Here, this is the velocity of the bodies minus the base flow
        velocity.
        """
        b = self._geometry.getVelocities()
        b0 = self.getBaseFlowBoundaryVelocities()
        b -= b0
        return b

    def updateOperators(self, time: float) -> None:
        """Update operators, for time-dependent models."""
        if self.bfTimeDependent():
            self._baseFlow.moveFlow(time)
        if self.geTimeDependent():
            self._geometry.moveBodies(time)
            self._regularizer.update()

    def getGrid(self) -> Grid:
        """Return a pointer to the associated Grid."""
        return self._grid

    def B(self, f: BoundaryVector, omega: Scalar) -> None:
        """Compute omega = B(f) as in (14)."""
        # NOTE(port): C++ signature is `B(const BoundaryVector& f, Scalar&
        # omega)` -- `omega` is an out-param, mutated in place (Curl's
        # out-param form writes directly into `omega._data`), matching the
        # C++ reference-mutation semantics.
        assert self._hasBeenInitialized
        q = self._regularizer.toFlux(f)
        Curl(q, omega)

    def C(self, omega: Scalar, f: BoundaryVector) -> None:
        """Compute f = C(omega) as in (14)."""
        # NOTE(port): C++ signature is `C(const Scalar& omega, BoundaryVector&
        # f)` -- `f` is an out-param. C++ `f = _regularizer.toBoundary(q)`
        # invokes BoundaryVector::operator=, copying into the caller's
        # object; `f.assign(...)` reproduces that copy-through-reference
        # (rather than `f = ...`, which in Python would only rebind the
        # local name and not mutate the caller's object).
        assert self._hasBeenInitialized
        q = Flux(self._grid)
        self.computeFluxWithoutBaseFlow(omega, q)
        f.assign(self._regularizer.toBoundary(q))

    def computeFluxWithoutBaseFlow(self, omega: Scalar, q: Flux) -> None:
        assert self._hasBeenInitialized
        streamfunction = self.vorticityToStreamfunction(omega)
        Curl(streamfunction, q)

    def getAlpha(self) -> float:
        """Return the constant alpha = 1/ReynoldsNumber."""
        return 1.0 / self._ReynoldsNumber

    def getAlphaBF(self) -> float:
        """Return the angle of attack of the baseflow."""
        return self._baseFlow.getAlpha()

    def computeFlux(self, omega: Scalar, q: Flux) -> None:
        """Compute flux q from vorticity omega, including base flow q0."""
        # NOTE(port): `q` is an out-param, mutated in place via
        # computeFluxWithoutBaseFlow (Curl's out-param form) and then `+=`
        # (Flux.__iadd__, also in place), matching C++ reference-mutation
        # semantics.
        assert self._hasBeenInitialized
        self.computeFluxWithoutBaseFlow(omega, q)
        q += self._baseFlow.getFlux()

    def refreshState(self, x: State) -> None:
        """Compute flux q from the vorticity omega, including base flow."""
        self.computeFlux(x.omega, x.q)

    def vorticityToStreamfunction(self, omega: Scalar) -> Scalar:
        """Given the vorticity omega, return the streamfunction psi.

        Assumes psi = 0 on the boundary, and does not add in potential flow
        solution.
        """
        assert self._hasBeenInitialized
        # Solve L psi = omega, with zero Dirichlet bc's
        psi = -1.0 * omega
        psi.coarsify()
        self._poisson.solve(psi, psi)
        return psi

    def getBaseFlow(self) -> BaseFlow:
        """Return the BaseFlow."""
        return self._baseFlow

    def getBaseFlowBoundaryVelocities(self) -> BoundaryVector:
        assert self._hasBeenInitialized
        velocity = self._regularizer.toBoundary(self._baseFlow.getFlux())
        return velocity
