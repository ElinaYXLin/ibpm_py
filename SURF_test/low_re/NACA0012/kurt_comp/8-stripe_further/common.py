"""
common.py

Shared helpers for kurt_comp/8-stripe_further -- the 5 follow-up tests
proposed at the end of ../7-stripe_investigation/README.md's "Proposals"
section. Directly imports (not redefines) ../7-stripe_investigation's
common.py, since this folder is a tight continuation of that one and
duplicating ~200 lines of the same grid/geometry/upstream-metric helpers
would just be a maintenance hazard -- any fix there should apply here
too. Only genuinely new helpers (state editing/saving for the regrow
test, eigendecomposition helpers) are added locally.
"""
import importlib.util
import os
import pathlib
import sys
import types

import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
os.chdir(REPO)

HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "8-stripe_further"
DATA = HERE / "data"
FIGS = HERE / "figures"
RUNS = HERE / "runs"
for d in (DATA, FIGS, RUNS):
    d.mkdir(parents=True, exist_ok=True)

KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp"
KURT1 = KURT / "1-paper_based"
KURT5 = KURT / "5-leading_edge"
KURT6 = KURT / "6-edges_further"
KURT7 = KURT / "7-stripe_investigation"

# Load ../7-stripe_investigation/common.py under a DISTINCT module name
# (not "common" -- this file is itself imported as "common" by test
# scripts here, and re-using that name would self-reference the
# still-initializing module instead of Kurt7's file).
_spec = importlib.util.spec_from_file_location("stripe7_common", KURT7 / "common.py")
_stripe7 = importlib.util.module_from_spec(_spec)
sys.modules["stripe7_common"] = _stripe7
_spec.loader.exec_module(_stripe7)

BASE_GEOM_DX002 = _stripe7.BASE_GEOM_DX002
CYLINDER_GEOM_DX002 = _stripe7.CYLINDER_GEOM_DX002
CPP_BIN = _stripe7.CPP_BIN
PY_RUNNER = _stripe7.PY_RUNNER
DOMAIN = _stripe7.DOMAIN
RE = _stripe7.RE
State = _stripe7.State
Geometry = _stripe7.Geometry
Grid = _stripe7.Grid
NavierStokesModel = _stripe7.NavierStokesModel
Scalar = _stripe7.Scalar
grid_xy = _stripe7.grid_xy
load_omega = _stripe7.load_omega
load_state = _stripe7.load_state
load_geom_points = _stripe7.load_geom_points
nearest_index = _stripe7.nearest_index
window = _stripe7.window
is_done = _stripe7.is_done
run_case = _stripe7.run_case
upstream_mask = _stripe7.upstream_mask
upstream_profile = _stripe7.upstream_profile
upstream_scalar_metrics = _stripe7.upstream_scalar_metrics
reach_L_up = _stripe7.reach_L_up
upstream_mask_2d_rotated = _stripe7.upstream_mask_2d_rotated
upstream_enstrophy_rotated = _stripe7.upstream_enstrophy_rotated

if "py_static" not in sys.modules:
    _pkg = types.ModuleType("py_static")
    _pkg.__path__ = [str(REPO / "py_static")]
    sys.modules["py_static"] = _pkg
from py_static.cholesky_solver import CholeskySolver  # noqa: E402
from py_static.boundary_vector import BoundaryVector  # noqa: E402


def build_projection_matrix(grid, geom, re, beta=0.5 * 0.01):
    """Build the dense projection matrix M (same construction as
    ../6-edges_further/testF_conditioning.py) for a given grid/geometry,
    and return (M, numPoints). `beta` matches Group F's convention
    (0.5*dt, dt=0.01's first-RK-stage value) -- only the matrix's
    *eigenvectors* are used by this folder's tests, and M's eigenvectors
    are the same regardless of the exact beta scaling (beta only scales
    eigenvalues uniformly), so this choice doesn't affect what's measured
    here."""
    model = NavierStokesModel(grid, geom, re)
    model.init()
    solver = CholeskySolver(grid, model, beta)
    numPoints = geom.getNumPoints()
    size = 2 * numPoints
    M = np.zeros((size, size))
    e = BoundaryVector(numPoints)
    x = BoundaryVector(numPoints)
    for j in range(size):
        e.assign(0)
        e.set(j, 1)
        solver.M(e, x)
        M[:, j] = x.flatten()
    return M, numPoints


def save_edited_state(state, out_path):
    """Save a (possibly edited) State object to a restart-format .bin
    file, for use as a `-ic` initial condition in a new run."""
    ok = state.save(str(out_path))
    return ok
