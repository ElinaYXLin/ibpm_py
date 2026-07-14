"""
compute_le_residual.py

Test 4 of the LE vorticity investigation (see ../README.md).

SURF_test/vortall/1-baseline/inner/compute_residuals.py already showed the
no-slip constraint residual ||C(omega) - b|| is at floating-point-roundoff
level for the vortall cylinder case, as a SINGLE scalar norm over all
boundary points. That can't say WHERE on the body the projection step is
straining hardest. This script reruns that same no-slip-residual
computation for the NACA0012 Re=500, alpha=0, dx=0.02 case (this suite's
baseline geometry/grid), but keeps the residual PER BOUNDARY POINT
(x,y components at each of the ~102 Lagrangian points around the airfoil)
instead of collapsing it to a norm, and records each point's (x,y) location
and arc-length s alongside it -- so the residual's spatial distribution can
be checked directly: if it peaks right at the leading edge, that's direct
evidence the projection step is genuinely straining hardest there (not
something silently wrong elsewhere in the domain).

Uses the internal py/ API directly (Grid/Geometry/NavierStokesModel/
NonlinearIBSolver), same pattern as compute_residuals.py.

Usage: python3 SURF_test/low_re/NACA0012/leading_edge_investigation/compute_le_residual.py
Output: SURF_test/low_re/NACA0012/leading_edge_investigation/data/le_residual_spatial.csv
"""
import pathlib
import sys
import types

import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg

from py.base_flow import BaseFlow  # noqa: E402
from py.geometry import Geometry  # noqa: E402
from py.grid import Grid  # noqa: E402
from py.ib_solver import NonlinearIBSolver  # noqa: E402
from py.navier_stokes_model import NavierStokesModel  # noqa: E402
from py.scheme import SchemeType  # noqa: E402
from py.state import State  # noqa: E402

DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"
GEOM_FILE = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"

NX, NY, NGRID = 300, 150, 1
LENGTH, XOFFSET, YOFFSET = 6.0, -2.0, -1.5
RE = 500.0
ALPHA_DEG = 0.0
DT = 0.01
# Re=500 alpha=0 is steady/laminar (see ../../README.md); 1500 steps (t=15)
# is well past the ~t=5-10 transient seen in this suite's force traces, so
# the residual is recorded at (and near) the converged steady state.
RECORD_STEPS = [0, 100, 500, 1500]


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    grid = Grid(NX, NY, NGRID, LENGTH, XOFFSET, YOFFSET, 0.0, 0.0)
    geom = Geometry()
    assert geom.load(str(GEOM_FILE)), f"failed to load {GEOM_FILE}"

    alpha_rad = np.deg2rad(ALPHA_DEG)
    q_potential = BaseFlow(grid, 1.0, alpha_rad)

    model = NavierStokesModel(grid, geom, RE, q_potential)
    solver = NonlinearIBSolver(grid, model, DT, SchemeType.RK3)
    model.init()
    solver.init()

    x = State(grid, geom.getNumPoints())
    x.omega.assign(0.0)
    x.f.assign(0.0)
    x.q.assign(0.0)
    geom.moveBodies(x.time)
    model.updateOperators(x.time)
    model.refreshState(x)

    # boundary-point coordinates + arc length (fixed for a stationary body)
    pts = geom.getPoints()
    n = geom.getNumPoints()
    px = pts._data[0 * n:1 * n].copy()
    py_ = pts._data[1 * n:2 * n].copy()
    seg = np.hypot(np.diff(np.append(px, px[0])), np.diff(np.append(py_, py_[0])))
    s = np.concatenate([[0.0], np.cumsum(seg)])[:-1]
    i_le = np.argmin(px)
    s_le = s[i_le]

    rows = []  # step, time, point_idx, x, y, s, d_to_le, resx, resy, res_mag

    def record(step, tsim):
        b = model.getConstraints()
        f_check = type(b)(b.getNumPoints())
        model.C(x.omega, f_check)
        resx = f_check._data[0 * n:1 * n] - b._data[0 * n:1 * n]
        resy = f_check._data[1 * n:2 * n] - b._data[1 * n:2 * n]
        res_mag = np.hypot(resx, resy)
        perim = s[-1] + seg[-1]
        d_le = np.abs(s - s_le)
        d_le = np.minimum(d_le, perim - d_le)
        for i in range(n):
            rows.append([step, tsim, i, px[i], py_[i], s[i], d_le[i],
                         resx[i], resy[i], res_mag[i]])
        print(f"  step {step} t={tsim:.2f}: max|res|={res_mag.max():.3e} "
              f"at point {int(np.argmax(res_mag))} (x={px[np.argmax(res_mag)]:.3f}, "
              f"d_to_LE={d_le[np.argmax(res_mag)]:.4f}), "
              f"mean|res|={res_mag.mean():.3e}", flush=True)

    record(0, x.time)
    steps_done = 0
    for target in RECORD_STEPS[1:]:
        while steps_done < target:
            solver.advance(x)
            steps_done += 1
        record(steps_done, x.time)

    header = "step,time,point_idx,x,y,s,d_to_le,resx,resy,res_mag"
    out_csv = DATA_DIR / "le_residual_spatial.csv"
    np.savetxt(out_csv, np.array(rows), delimiter=",", header=header, comments="")
    print(f"wrote {out_csv} ({len(rows)} rows, {n} points x {len(RECORD_STEPS)} snapshots)")


if __name__ == "__main__":
    main()
