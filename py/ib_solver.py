# ib_solver.py
#
# Python port of src/IBSolver.h / src/IBSolver.cc
#
# Immersed-boundary time-integration schemes. IBSolver is an abstract base
# class implementing the fractional-step time advance; concrete subclasses
# (NonlinearIBSolver, LinearizedIBSolver, AdjointIBSolver,
# LinearizedPeriodicIBSolver, SFDSolver) supply the nonlinear term N(x).
#
# Original author: Clancy Rowley
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES (see also inline NOTE(port) comments):
#
#   1. Abstract base class.  C++ makes IBSolver abstract via the pure virtual
#      `Scalar N(const State&) const = 0`. Here IBSolver derives from abc.ABC
#      and N is an @abstractmethod, so it cannot be instantiated directly.
#
#   2. Constructor overload collapse + tol argument order.  Every IBSolver
#      subclass has two C++ constructors, one without a tolerance (default
#      1e-7) and one with. In the *derived* classes the tol argument sits in
#      the MIDDLE of the C++ parameter list (before baseFlow / x0periodic),
#      which Python cannot reproduce positionally alongside a no-tol overload.
#      They are collapsed into a single constructor with `tol` as a trailing
#      keyword argument defaulting to 1e-7 (e.g.
#      LinearizedIBSolver(grid, model, dt, scheme, baseFlow, tol=1e-7)).
#      Flagged because it reorders the tol parameter relative to C++.
#
#   3. `vector<ProjectionSolver*> _solver`.  Becomes a Python list. The C++
#      constructor sizes it to nsteps null pointers, then createAllSolvers()
#      fills it. Reproduced with `[None]*nsteps` + createAllSolvers().
#
#   4. Destructor / deleteAllSolvers.  C++ ~IBSolver() calls deleteAllSolvers()
#      to `delete` the owned ProjectionSolver pointers. Python is garbage
#      collected, so no __del__ is defined; deleteAllSolvers() is retained
#      (part of the protected interface) but only clears the Python references.
#      Note the C++ setTol()/createAllSolvers() overwrite _solver[i] WITHOUT
#      first deleting the old solver (a leak in C++); the Python version simply
#      rebinds, letting the GC reclaim the previous object.
#
#   5. Out-parameter State.  advance()/advanceSubstep() mutate the passed-in
#      State `x` in place (C++ takes `State&`); Scalar/State copy-assignments
#      (`_Nprev = nonlinear`, `_xhat = x`) are done via `.assign(...)` to match
#      C++ value-copy semantics (see state.py / scalar.py).
#
#   6. sprintf formats.  IBSolver.load/save build "..._%02d" suffixes and
#      SFDSolver.saveFilteredState applies a caller-supplied printf format to a
#      timestep; ported with Python `f"{i+1:02d}"` and `%`-formatting
#      respectively.
#
#   7. cerr/cout progress messages preserved (cerr -> sys.stderr,
#      cout -> print).
#
#   8. JAX-readiness.  All numerics route through Scalar/Flux/BoundaryVector
#      and the VectorOperations / ProjectionSolver ports; this file only
#      orchestrates them, so nothing here blocks a later JAX port.
# ---------------------------------------------------------------------------

from __future__ import annotations

import abc
import sys
from typing import List, Optional

from .boundary_vector import BoundaryVector
from .cholesky_solver import CholeskySolver
from .conjugate_gradient_solver import ConjugateGradientSolver
from .grid import Grid
from .navier_stokes_model import NavierStokesModel
from .projection_solver import ProjectionSolver
from .scalar import Scalar
from .scheme import Scheme
from .state import State
from .vector_operations import CrossProduct, Curl, Laplacian


# =============== #
#   Base class    #
# =============== #

class IBSolver(abc.ABC):
    """Abstract base class: immersed-boundary fractional-step time integrator."""

    def __init__(
        self,
        grid: Grid,
        model: NavierStokesModel,
        dt: float,
        scheme: "Scheme.SchemeType",
        tol: float = 1e-7,
    ) -> None:
        # NOTE(port): the two C++ constructors (with/without tol) are collapsed
        # into one with tol defaulting to 1e-7 (port note 2).
        self._grid: Grid = grid
        self._scheme: Scheme = Scheme(scheme)
        self._name: str = self._scheme.name()
        self._dt: float = dt
        self._model: NavierStokesModel = model
        self._Nprev: Scalar = Scalar(grid)
        self._Ntemp: Scalar = Scalar(grid)
        self._oldSaved: bool = False
        self._solver: List[Optional[ProjectionSolver]] = [None] * self._scheme.nsteps()
        self._tol: float = tol
        self.createAllSolvers()

    def getName(self) -> str:
        return self._name

    def getTimestep(self) -> float:
        return self._dt

    def init(self) -> None:
        for i in range(self._scheme.nsteps()):
            self._solver[i].init()

    def reset(self) -> None:
        self._oldSaved = False

    def load(self, basename: str) -> bool:
        successInit = False
        successTemp = True
        for i in range(self._scheme.nsteps()):
            num = f"{i + 1:02d}"
            filename = basename + "_" + num
            successTemp = bool(self._solver[i].load(filename)) and successTemp
            if i == 0:
                successInit = True
        return successInit and successTemp

    def save(self, basename: str) -> bool:
        successInit = False
        successTemp = True
        for i in range(self._scheme.nsteps()):
            num = f"{i + 1:02d}"
            filename = basename + "_" + num
            successTemp = bool(self._solver[i].save(filename)) and successTemp
            if i == 0:
                successInit = True
        return successInit and successTemp

    def createAllSolvers(self) -> None:
        for i in range(self._scheme.nsteps()):
            self._solver[i] = self.createSolver(
                (self._scheme.an(i) + self._scheme.bn(i)) * self._dt
            )

    def deleteAllSolvers(self) -> None:
        # NOTE(port): C++ `delete`s each owned solver; Python GC handles that,
        # so we only clear the references (port note 4).
        for i in range(len(self._solver)):
            self._solver[i] = None

    def createSolver(self, beta: float) -> ProjectionSolver:
        # Check whether all bodies are stationary:
        #   If so, return a CholeskySolver; if not, a ConjugateGradientSolver.
        if self._model.geTimeDependent():
            print(
                "Using ConjugateGradient solver for projection step\n"
                f"  tolerance = {self._tol}",
                file=sys.stderr,
            )
            return ConjugateGradientSolver(self._grid, self._model, beta, self._tol)
        else:
            print("Using Cholesky solver for projection step", file=sys.stderr)
            return CholeskySolver(self._grid, self._model, beta)

    def setTol(self, tol: float) -> None:
        self._tol = tol
        self.createAllSolvers()

    def advance(self, x: State, Bu: Optional[Scalar] = None) -> None:
        # NOTE(port): overload collapse of advance(State&) and
        # advance(State&, const Scalar& Bu). When Bu is given, the nonlinear
        # term has Bu added (the forced form); otherwise it is N(x) alone.
        for i in range(self._scheme.nsteps()):
            if Bu is None:
                nonlinear = self.N(x)
            else:
                nonlinear = self.N(x) + Bu
            self.advanceSubstep(x, nonlinear, i)

        x.time += self._dt
        x.timestep += 1

    def advanceSubstep(self, x: State, nonlinear: Scalar, i: int) -> None:
        # If the body is moving, update the positions of the bodies
        if self._model.isTimeDependent():
            self._model.updateOperators(x.time + self._scheme.cn(i) * self._dt)

        # Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
        a = Laplacian(x.omega)
        a *= 0.5 * self._model.getAlpha() * (self._scheme.an(i) + self._scheme.bn(i))
        a += self._scheme.an(i) * nonlinear

        if self._scheme.bn(i) != 0:
            # for ab2
            if self._oldSaved is False:
                self._Nprev.assign(nonlinear)

            a += self._scheme.bn(i) * self._Nprev

        a *= self._dt
        a += x.omega

        # Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
        b = self._model.getConstraints()

        # Call the ProjectionSolver to determine the vorticity and forces
        self._solver[i].solve(a, b, x.omega, x.f)

        # Update the state, for instance to compute the corresponding flux
        self._model.refreshState(x)
        self._Nprev.assign(nonlinear)

        if self._oldSaved is False:
            self._oldSaved = True

    # -- protected pure-virtual --------------------------------------------

    @abc.abstractmethod
    def N(self, x: State) -> Scalar:
        """Compute the nonlinear term N(x)."""
        ...


# =============== #
# Derived classes #
# =============== #

class NonlinearIBSolver(IBSolver):
    def __init__(
        self,
        grid: Grid,
        model: NavierStokesModel,
        dt: float,
        scheme: "Scheme.SchemeType",
        tol: float = 1e-7,
    ) -> None:
        super().__init__(grid, model, dt, scheme, tol)

    def N(self, x: State) -> Scalar:
        v = CrossProduct(x.q, x.omega)
        g = Curl(v)
        return g


class LinearizedIBSolver(IBSolver):
    def __init__(
        self,
        grid: Grid,
        model: NavierStokesModel,
        dt: float,
        scheme: "Scheme.SchemeType",
        baseFlow: State,
        tol: float = 1e-7,
    ) -> None:
        # NOTE(port): tol reordered to a trailing default (port note 2).
        super().__init__(grid, model, dt, scheme, tol)
        self._x0: State = State(baseFlow)

    def N(self, x: State) -> Scalar:
        v = CrossProduct(self._x0.q, x.omega)
        v += CrossProduct(x.q, self._x0.omega)
        g = Curl(v)
        return g


class AdjointIBSolver(IBSolver):
    def __init__(
        self,
        grid: Grid,
        model: NavierStokesModel,
        dt: float,
        scheme: "Scheme.SchemeType",
        baseFlow: State,
        tol: float = 1e-7,
    ) -> None:
        # NOTE(port): tol reordered to a trailing default (port note 2).
        super().__init__(grid, model, dt, scheme, tol)
        self._x0: State = State(baseFlow)

    def N(self, x: State) -> Scalar:
        g = Laplacian(CrossProduct(self._x0.q, x.q))
        g -= Curl(CrossProduct(x.q, self._x0.omega))
        return g


class LinearizedPeriodicIBSolver(IBSolver):
    """Navier-Stokes equations linearized about a periodic orbit."""

    def __init__(
        self,
        grid: Grid,
        model: NavierStokesModel,
        dt: float,
        scheme: "Scheme.SchemeType",
        x0periodic: "List[State]",
        period: int,
        tol: float = 1e-7,
    ) -> None:
        # NOTE(port): tol reordered to a trailing default (port note 2).
        super().__init__(grid, model, dt, scheme, tol)
        # NOTE(port): C++ stores `const vector<State> _x0periodic` by value
        # (a copy of the passed vector). We keep a shallow list copy here;
        # the States within are used read-only in N(), matching the C++
        # `const` usage.
        self._x0periodic: List[State] = list(x0periodic)
        self._period: int = period
        assert self._period == len(x0periodic)

    def N(self, x: State) -> Scalar:
        k = x.timestep % self._period
        print(f"At time step {x.timestep}, phase k = {k}")
        v = CrossProduct(self._x0periodic[k].q, x.omega)
        v += CrossProduct(x.q, self._x0periodic[k].omega)
        g = Curl(v)
        return g


# =========== #
# SFD methods #
# =========== #

class SFDSolver(IBSolver):
    """Selective frequency damping solver."""

    def __init__(
        self,
        grid: Grid,
        model: NavierStokesModel,
        dt: float,
        scheme: "Scheme.SchemeType",
        Delta: float,
        chi: float,
    ) -> None:
        # NOTE(port): SFDSolver has only one C++ constructor (no tol overload),
        # so the base tol keeps its default 1e-7.
        super().__init__(grid, model, dt, scheme)
        self._Delta: float = Delta          # inverse of cutoff frequency
        self._chi: float = chi              # sfd gain
        self._xhat: State = State(self._grid, self._model.getNumPoints())
        self._omegaTemp: Scalar = Scalar(grid)
        self._rhsPrev: Scalar = Scalar(grid)
        self._xhatSaved: bool = False
        self._rhsSaved: bool = False

    def N(self, x: State) -> Scalar:
        v = CrossProduct(x.q, x.omega)
        g = Curl(v)
        temp = Scalar(x.omega)  # because x is const here...hmmm
        g -= self._chi * (temp - self._xhat.omega)
        return g

    def advanceSubstep(self, x: State, nonlinear: Scalar, i: int) -> None:
        assert x.time == self._xhat.time

        # Initialize _xhat if necessary, save current vorticity field
        if self._xhatSaved is False:
            self._xhat.assign(x)
            self._xhatSaved = True
        self._omegaTemp.assign(x.omega)

        # Advance state x
        IBSolver.advanceSubstep(self, x, nonlinear, i)

        # Advance state _xhat
        rhs = (self._omegaTemp - self._xhat.omega) / self._Delta
        a = self._scheme.an(i) * rhs

        if self._scheme.bn(i) != 0:
            if self._rhsSaved is False:
                self._rhsPrev.assign(rhs)

            a += self._scheme.bn(i) * self._rhsPrev

        a *= self._dt
        self._xhat.omega += a

        if i == self._scheme.nsteps() - 1:
            self._xhat.time += self._dt
            self._xhat.timestep += 1

        self._rhsPrev.assign(rhs)

        if self._rhsSaved is False:
            self._rhsSaved = True

    def saveFilteredState(self, outdir: str, name: str, numDigitInFileName: str) -> None:
        # NOTE(port): C++ builds a printf format string and sprintf's the
        # timestep into it; reproduced with Python %-formatting (port note 6).
        formatString = outdir + name + numDigitInFileName + ".bin" + "_xhat"
        filename = formatString % self._xhat.timestep
        self._xhat.save(filename)

    def loadFilteredState(self, icFile: str) -> None:
        xhatFile = icFile + "_xhat"
        self._xhat.omega.assign(0.0)
        self._xhat.f.assign(0.0)
        self._xhat.q.assign(0.0)
        if xhatFile != "_xhat":
            print(f"Loading initial condition from file: {xhatFile}")
            if not self._xhat.load(xhatFile):
                print("  (failed: setting xhat = x)")
        else:
            print("Setting xhat = x")

        self._xhatSaved = True
