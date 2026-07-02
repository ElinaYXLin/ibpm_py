# ibpm.py
#
# Python port of src/ibpm.cc (+ src/ibpm.h's IBPM_VERSION macro)
#
# Sample main routine for IBFS code: set up a timestepper and advance the
# flow in time.
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES (see also inline NOTE(port) comments):
#
#   1. Entry point signature. C++ is `int main(int argc, char* argv[])`.
#      Ported as `main(argv: Optional[List[str]] = None) -> int`, defaulting
#      to `sys.argv` when omitted (the idiomatic Python equivalent), rather
#      than requiring a redundant explicit `argc`. `ParmParser` itself keeps
#      the C++ `(argc, argv)` signature verbatim (see parm_parser.py) --
#      `argc = len(argv)` is computed once here and passed through.
#
#   2. `cout`/`cerr` -> `print`/`sys.stderr.write`. Formatted output (e.g.
#      `setw`) is approximated with Python string formatting rather than
#      reproduced byte-for-byte; this only affects console cosmetics, not
#      program behavior.
#
#   3. `x.omega = 0.` / `x.f = 0.` / `x.q = 0.` (C++ `operator=(double)`,
#      zeroing an existing Scalar/BoundaryVector/Flux in place) are ported as
#      `x.omega.assign(0.0)` etc., *not* `x.omega = 0.0` -- the latter would
#      just rebind the Python attribute to the float `0.0`, discarding the
#      underlying Scalar/Flux/BoundaryVector object entirely, which is not
#      what the C++ assignment operator does.
#
#   4. `x00 = x0[0]` (C++ `State::operator=`, a deep copy) is ported as
#      `x00.assign(x0[0])`, matching the `.assign()` convention used
#      throughout this port for C++ copy-assignment.
#
#   5. Preserved-as-is oddities in the original C++ source (not "fixed"
#      here, per the faithful-port instruction):
#        - The "Writing energy every ... step(s)" message prints `iForce`
#          instead of `iEnergy` (see the `iEnergy > 0` branch below).
#        - `mag` is computed but the resulting `BaseFlow` magnitude
#          parameter naming shadows nothing unusual; kept as `magnitude`.
#        - `x.timestep % iRestart` inside the SFD branch divides/mods by
#          `iRestart`, which is a genuine crash (`ZeroDivisionError` in
#          Python; undefined behavior / SIGFPE in C++) if the user passes
#          `-restart 0` together with `-model sfd`. Not guarded against,
#          matching the original.
#
#   6. `NavierStokesModel* model = NULL; IBSolver* solver = NULL;
#      SFDSolver* SFDsolver = NULL;` are ported as plain `None`-initialized
#      Python names (`model: Optional[NavierStokesModel] = None`, etc.)
#      rather than reproducing C++ pointer semantics; the later
#      `assert model is not None` / `assert solver is not None` checks are
#      kept, matching the original asserts.
#
#   7. `delete solver;` at the end of `main` has no Python equivalent
#      (garbage-collected); omitted.
#
#   8. JAX-readiness: this file only orchestrates already-ported
#      Grid/State/NavierStokesModel/IBSolver/Output objects and does no
#      numerics of its own, so nothing here blocks a later JAX port.
# ---------------------------------------------------------------------------

from __future__ import annotations

import math
import os
import sys
from typing import List, Optional

from .base_flow import BaseFlow
from .flux import Flux
from .geometry import Geometry
from .grid import Grid
from .ib_solver import (
    AdjointIBSolver,
    IBSolver,
    LinearizedIBSolver,
    LinearizedPeriodicIBSolver,
    NonlinearIBSolver,
    SFDSolver,
)
from .logger import Logger
from .motion import Motion
from .navier_stokes_model import NavierStokesModel
from .output_energy import OutputEnergy
from .output_force import OutputForce
from .output_restart import OutputRestart
from .output_tecplot import OutputTecplot
from .parm_parser import ParmParser
from .scheme import Scheme, SchemeType
from .state import State
from .utils import AddSlashToPath, MakeLowercase
from .vector_operations import InnerProduct

IBPM_VERSION = "1.0"


class ModelType:
    """Port of the C++ file-local `enum ModelType { LINEAR, NONLINEAR,
    ADJOINT, LINEARPERIODIC, SFD, INVALID };`.

    NOTE(port): kept as a plain class of int constants (rather than
    `enum.IntEnum`) purely as a stylistic choice matching this being a
    private, file-local C++ enum with no dependents outside ibpm.cc; either
    representation would be equally faithful.
    """

    LINEAR = 0
    NONLINEAR = 1
    ADJOINT = 2
    LINEARPERIODIC = 3
    SFD = 4
    INVALID = 5


def str2model(modelName: str) -> int:
    """Return the type of model specified in the string modelName."""
    modelName = MakeLowercase(modelName)

    if modelName == "nonlinear":
        modelType = ModelType.NONLINEAR
    elif modelName == "linear":
        modelType = ModelType.LINEAR
    elif modelName == "adjoint":
        modelType = ModelType.ADJOINT
    elif modelName == "linearperiodic":
        modelType = ModelType.LINEARPERIODIC
    elif modelName == "sfd":
        modelType = ModelType.SFD
    else:
        sys.stderr.write(f"Unrecognized model: {modelName}\n")
        modelType = ModelType.INVALID
    return modelType


def str2scheme(integratorType: str) -> SchemeType:
    """Return the integration scheme specified in the string
    integratorType."""
    schemeName = MakeLowercase(integratorType)
    if schemeName == "euler":
        schemeType = SchemeType.EULER
    elif schemeName == "ab2":
        schemeType = SchemeType.AB2
    elif schemeName == "rk3":
        schemeType = SchemeType.RK3
    elif schemeName == "rk3b":
        schemeType = SchemeType.RK3b
    else:
        sys.stderr.write(f"Unrecognized integration scheme: {schemeName}")
        sys.stderr.write("    Exiting program.\n")
        sys.exit(1)
    return schemeType


def main(argv: Optional[List[str]] = None) -> int:
    """Main routine for IBFS code: set up a timestepper and advance the
    flow in time.

    NOTE(port): see module-level judgment-call note 1 for the `argv`
    default-argument vs. C++'s explicit `(argc, argv)` signature.
    """
    if argv is None:
        argv = sys.argv

    print(f"Immersed Boundary Projection Method (IBPM), version {IBPM_VERSION}\n")

    # Get parameters
    parser = ParmParser(len(argv), argv)
    helpFlag = parser.getFlag("h", "print this help message and exit")

    # Output parameters
    name = parser.getString("name", "run name", "ibpm")
    outdir = parser.getString("outdir", "directory for saving output", ".")
    iTecplot = parser.getInt("tecplot", "if >0, write a Tecplot file every n timesteps", 100)
    TecplotAllGrids = parser.getBool("tecplotallgrids", "Tecplot output for all grids, or not", False)
    iRestart = parser.getInt("restart", "if >0, write a restart file every n timesteps", 100)
    iForce = parser.getInt("force", "if >0, write forces every n timesteps", 1)
    iEnergy = parser.getInt("energy", "if >0, write energy every n timesteps", 0)
    numDigitInFileName = parser.getString(
        "numdigfilename", "number of digits for time representation in filename", "%05d"
    )

    # Grid parameters
    nx = parser.getInt("nx", "number of gridpoints in x-direction", 200)
    ny = parser.getInt("ny", "number of gridpoints in y-direction", 200)
    ngrid = parser.getInt("ngrid", "number of grid levels for multi-domain scheme", 1)
    length = parser.getDouble("length", "length of finest domain in x-dir", 4.0)
    xOffset = parser.getDouble("xoffset", "x-coordinate of left edge of finest domain", -2.0)
    yOffset = parser.getDouble("yoffset", "y-coordinate of bottom edge of finest domain", -2.0)
    xShift = parser.getDouble("xshift", "percentage offset between grid levels in x-direction", 0.0)
    yShift = parser.getDouble("yshift", "percentage offset between grid levels in y-direction", 0.0)
    alpha = parser.getDouble("alpha", "angle of attack of base flow", 0.0)

    # Simulation parameters
    geomFile = parser.getString("geom", "filename for reading geometry", name + ".geom")
    ubf = parser.getBool("ubf", "Use unsteady base flow, or not", False)
    Reynolds = parser.getDouble("Re", "Reynolds number", 100.0)
    modelName = parser.getString(
        "model", "type of model (linear, nonlinear, adjoint, linearperiodic, sfd)", "nonlinear"
    )
    baseFlow = parser.getString("baseflow", "base flow for linear/adjoint model", "")

    # Initial condition
    icFile = parser.getString("ic", "initial condition filename", "")
    resetTime = parser.getBool("resettime", "Reset time when subtracting ic by baseflow (1/0(true/false))", False)
    subtractBaseflow = parser.getBool("subbaseflow", "Subtract ic by baseflow (1/0(true/false))", False)

    # Integration parameters
    dt = parser.getDouble("dt", "timestep", 0.02)
    numSteps = parser.getInt("nsteps", "number of timesteps to compute", 250)
    integratorType = parser.getString("scheme", "timestepping scheme (euler,ab2,rk3,rk3b)", "rk3")

    # Linear-periodic model
    period = parser.getInt("period", "period of periodic baseflow", 1)
    periodStart = parser.getInt("periodstart", "start time of periodic baseflow", 0)
    periodBaseFlowName = parser.getString(
        "pbaseflowname",
        "name of periodic baseflow, e.g. 'flow/ibpmperiodic%05d.bin', "
        "with '%05d' as time, decided by periodstart/period",
        "",
    )

    # SFD
    chi = parser.getDouble("chi", "sfd gain", 0.02)
    Delta = parser.getDouble("Delta", "sfd cutoff frequency", 15.0)

    # NOTE(port) -- NOT IN C++ ibpm.cc: opt-in Cholesky regularization. When
    # > 0, factor (M + cholreg*I) instead of M, which restores an SPD matrix
    # for over-resolved boundaries where the plain (C++-faithful) factorization
    # produces NaNs. Default 0.0 reproduces the C++ behavior exactly.
    cholReg = parser.getDouble(
        "cholreg",
        "Tikhonov diagonal regularization for the Cholesky projection solver "
        "(0 = faithful to C++; try 1e-8 if the factorization yields NaNs)",
        0.0,
    )

    modelType = str2model(modelName)
    schemeType = str2scheme(integratorType)

    if (not parser.inputIsValid()) or modelType == ModelType.INVALID or helpFlag:
        parser.printUsage(sys.stderr)
        sys.exit(1)

    # modify this long if statement?
    if (modelType != ModelType.NONLINEAR) and (modelType != ModelType.SFD):
        if modelType != ModelType.LINEARPERIODIC and baseFlow == "":
            print("ERROR: for linear or adjoint models, must specify a base flow")
            sys.exit(1)
        elif modelType != ModelType.LINEARPERIODIC and periodBaseFlowName != "":
            print("WARNING: for linear or adjoint models, a periodic base flow is not needed")
            sys.exit(1)
        elif modelType == ModelType.LINEARPERIODIC and periodBaseFlowName == "":
            print("ERROR: for linear periodic model, must specify a periodic base flow")
            sys.exit(1)
        elif modelType == ModelType.LINEARPERIODIC and baseFlow != "":
            print("WARNING: for linear periodic model, a single baseflow is not needed")
            sys.exit(1)

    # create output directory if not already present
    outdir = AddSlashToPath(outdir)
    # NOTE(port): C++ `mkdir(outdir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO)`
    # ignores its return value, so a pre-existing directory (or any other
    # mkdir failure) is silently ignored; `os.mkdir` raises `OSError` in
    # those cases, so the call is wrapped to match the "ignore failure"
    # behavior. `os.mkdir` (not `os.makedirs`) matches C's single-level
    # `mkdir` (it does not create missing parent directories).
    try:
        os.mkdir(outdir, 0o777)
    except OSError:
        pass

    # output command line arguments
    cmd = parser.getParameters()
    print(f"Command:\n{cmd}\n")
    parser.saveParameters(outdir + name + ".cmd")

    # Name of this run
    print(f"Run name: {name}\n")

    # Setup grid
    print(
        "Grid parameters:\n"
        f"    nx      {nx}\n"
        f"    ny      {ny}\n"
        f"    ngrid   {ngrid}\n"
        f"    length  {length}\n"
        f"    xoffset {xOffset}\n"
        f"    yoffset {yOffset}\n"
        f"    xshift  {xShift}\n"
        f"    yshift  {yShift}\n"
    )
    grid = Grid(nx, ny, ngrid, length, xOffset, yOffset, xShift, yShift)

    # Setup geometry
    geom = Geometry()
    print(f"Reading geometry from file {geomFile}")
    if geom.load(geomFile):
        print(f"    {geom.getNumPoints()} points on the boundary\n")
    else:
        sys.exit(-1)

    # Setup equations to solve
    print(f"Reynolds number = {Reynolds}\n")
    print("Setting up Immersed Boundary Solver...", end="", flush=True)
    magnitude = 1.0
    pi = 4.0 * math.atan(1.0)
    alpha = alpha * pi / 180.0
    xC = 0.0
    yC = 0.0
    q_potential = BaseFlow(grid, magnitude, alpha)
    # See if unsteady base flow can be used. Only implemented for a single
    # RigidBody in motion. In the future, have a function
    # geom.ubfEligible() that will make sure that the first RigidBody is
    # moving.
    if (not geom.isStationary()) and (geom.getNumBodies() == 1) and ubf:
        m: Optional[Motion] = geom.transferMotion()  # pull motion from first RigidBody object
        xC, yC = geom.transferCenter()  # pull center of motion from RigidBody object
        q_potential.setMotion(m)
        q_potential.setCenter(xC, yC)
    if ubf and (geom.getNumBodies() != 1):
        print("Unsteady base flow is only supported for a single moving body.  Exiting program.")
        sys.exit(1)

    model: Optional[NavierStokesModel] = None
    solver: Optional[IBSolver] = None
    SFDsolver: Optional[SFDSolver] = None
    x00 = State(grid, geom.getNumPoints())

    if modelType == ModelType.NONLINEAR:
        model = NavierStokesModel(grid, geom, Reynolds, q_potential)
        solver = NonlinearIBSolver(grid, model, dt, schemeType)
    elif modelType == ModelType.LINEAR:
        if not x00.load(baseFlow):
            print("baseflow failed to load.  Exiting program.")
            sys.exit(1)
        model = NavierStokesModel(grid, geom, Reynolds)
        solver = LinearizedIBSolver(grid, model, dt, schemeType, x00)
    elif modelType == ModelType.ADJOINT:
        if not x00.load(baseFlow):
            print("baseflow failed to load.  Exiting program.")
            sys.exit(1)
        model = NavierStokesModel(grid, geom, Reynolds)
        solver = AdjointIBSolver(grid, model, dt, schemeType, x00)
    elif modelType == ModelType.LINEARPERIODIC:
        # load periodic baseflow files
        x0 = [State(x00) for _ in range(period)]
        pbf = periodBaseFlowName
        for i in range(period):
            pbffilename = pbf % (i + periodStart)
            x0[i].load(pbffilename)
        x00.assign(x0[0])
        model = NavierStokesModel(grid, geom, Reynolds)
        solver = LinearizedPeriodicIBSolver(grid, model, dt, schemeType, x0, period)
    elif modelType == ModelType.SFD:
        print("SFD parameters:")
        print(f"    chi =   {chi}")
        print(f"    Delta = {Delta}\n")
        model = NavierStokesModel(grid, geom, Reynolds, q_potential)
        SFDsolver = SFDSolver(grid, model, dt, schemeType, Delta, chi)
        solver = SFDsolver
    elif modelType == ModelType.INVALID:
        print("ERROR: must specify a valid modelType")
        sys.exit(1)

    assert model is not None
    assert solver is not None
    if modelType == ModelType.SFD:
        assert chi != 0
        assert SFDsolver is not None
    # NOTE(port) -- NOT IN C++: enable opt-in Cholesky regularization if the
    # user passed -cholreg > 0. Rebuilds the projection solvers with the shift.
    if cholReg > 0.0:
        print(f"Using Cholesky regularization {cholReg}", file=sys.stderr)
        solver.setCholeskyRegularization(cholReg)
    # NOTE: still need to initialize model, but wait until after loading the
    #       initial condition, so we know what the initial time is, for
    #       moving the bodies

    # Load initial condition
    x = State(grid, geom.getNumPoints())
    x.omega.assign(0.0)
    x.f.assign(0.0)
    x.q.assign(0.0)
    if icFile != "":
        print(f"Loading initial condition from file: {icFile}")
        if not x.load(icFile):
            print("    (failed: using zero initial condition)")
        if subtractBaseflow:
            print("    Subtracting initial condition by baseflow to form a linear initial perturbation")
            if modelType != ModelType.NONLINEAR:
                assert x.q.Ngrid() == x00.q.Ngrid()
                assert x.omega.Ngrid() == x00.omega.Ngrid()
                x.q -= x00.q
                x.omega -= x00.omega
                x.f.assign(0.0)
            else:
                print("Flag subbaseflow should be true only for linear cases")
                sys.exit(1)

        if modelType == ModelType.SFD:
            SFDsolver.loadFilteredState(icFile)

    else:
        print("Using zero initial condition")

    if resetTime:
        x.timestep = 0
        x.time = 0.0

    # update the geometry to the current time
    geom.moveBodies(x.time)

    # Initialize model and timestepper
    model.init()
    print(f"using {solver.getName()} timestepper")
    print(f"    dt = {dt}\n")
    if not solver.load(outdir + name):
        # Set the tolerance for a ConjugateGradient solver below
        # Otherwise default is tol = 1e-7
        # solver.setTol( 1e-8 )
        solver.init()
        solver.save(outdir + name)

    # Calculate flux for state, in case only vorticity was saved
    if not q_potential.isStationary():
        q_potential.setAlphaMag(x.time)
        alpha = q_potential.getAlpha()
    model.updateOperators(x.time)
    model.refreshState(x)

    print(f"\nInitial timestep = {x.timestep}\n")

    # Setup output routines
    tecplot = OutputTecplot(
        outdir + name + numDigitInFileName + ".plt",
        "Test run, step" + numDigitInFileName,
        TecplotAllGrids,
    )
    if TecplotAllGrids:
        tecplot.setFilename(outdir + name + numDigitInFileName + "_g%01d.plt")
    restart = OutputRestart(outdir + name + numDigitInFileName + ".bin")
    force = OutputForce(outdir + name + ".force")
    energy = OutputEnergy(outdir + name + ".energy")

    logger = Logger()
    # Output Tecplot file every timestep
    if iTecplot > 0:
        print(f"Writing Tecplot file every {iTecplot} step(s)")
        logger.addOutput(tecplot, iTecplot)
    if iRestart > 0:
        print(f"Writing restart file every {iRestart} step(s)")
        logger.addOutput(restart, iRestart)
    if iForce > 0:
        print(f"Writing forces every {iForce} step(s)")
        logger.addOutput(force, iForce)
    if iEnergy > 0:
        # NOTE(port): preserved as-is -- the C++ message prints `iForce`
        # here, not `iEnergy` (see module-level judgment-call note 5).
        print(f"Writing energy every {iForce} step(s)")
        logger.addOutput(energy, iEnergy)
    print()
    logger.init()
    logger.doOutput(q_potential, x)

    print(f"Integrating for {numSteps} steps")
    for i in range(1, numSteps + 1):
        print(f"\nstep {i}")
        xtemp = State(x)  # For SFD norm calculation
        solver.advance(x)
        xF, yF = x.computeNetForce()
        # If there is an unsteady base flow, transform body frame normal and
        # parallel forces into lab frame lift and drag
        if not q_potential.isStationary():
            q_potential.setAlphaMag(x.time)
            alpha = q_potential.getAlpha()
        drag = xF * math.cos(alpha) + yF * math.sin(alpha)
        lift = xF * -1.0 * math.sin(alpha) + yF * math.cos(alpha)
        print(f"    x force: {drag * 2:16g}, y force: {lift * 2:16g}\n")
        logger.doOutput(q_potential, x)

        # For SFD
        if modelType == ModelType.SFD:
            # Inner product of fluxes is equal to inner product of
            # vorticity (with weighted inner product for latter)
            dq = xtemp.q - x.q
            q = math.sqrt(InnerProduct(x.q, x.q))
            twoNorm = math.sqrt(InnerProduct(dq, dq)) / (q * dt)

            # NOTE(port): see module-level judgment-call note 5 -- this
            # divides/mods by `iRestart`, which crashes if `-restart 0` is
            # combined with `-model sfd` (matching the original C++).
            if (x.timestep % iRestart == 0) and (chi != 0.0):
                SFDsolver.saveFilteredState(outdir, name, numDigitInFileName)

            print(f"    ||dx||/||x||/dt = {twoNorm:13g}")

    logger.cleanup()

    return 0


if __name__ == "__main__":
    sys.exit(main())
