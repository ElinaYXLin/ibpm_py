"""
common.py

Shared helpers for the kurt_comp/5-leading_edge investigation (why are there
striped/ringing artifacts at the leading and trailing edge of the vorticity
fields, e.g. kurt_comp/1-paper_based/figures/wake_steady_paperframe.png).

Uses py_static (never modified) via its internal API, the same pattern as
../../2-leading_edge_investigation/compute_le_residual.py used for the old
py/ port. All new runs and analysis here use py_static/cpp_static exclusively
(kurt_comp convention), never src/ or py/.
"""
import os
import pathlib
import sys
import types

import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
# .geom files store their "raw <path>" line relative to the repo root (see
# run_kurt_suite.py's subprocess calls, which always set cwd=REPO) -- the
# internal Geometry().load() used directly in this folder's scripts needs
# the same cwd, or it silently loads zero points instead of raising.
os.chdir(REPO)
HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "5-leading_edge"
DATA = HERE / "data"
FIGS = HERE / "figures"
GEOMDIR = HERE / "geom"
RUNS = HERE / "runs"
for d in (DATA, FIGS, GEOMDIR, RUNS):
    d.mkdir(parents=True, exist_ok=True)

KURT1 = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "1-paper_based"
BASE_GEOM_DX002 = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
NACA0012_DAT = REPO / "SURF_test" / "low_re" / "NACA0012" / "1-basics" / "naca0012.dat.txt"

CPP_BIN = REPO / "build_static" / "ibpm"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"

sys.path.insert(0, str(REPO))
_pkg = types.ModuleType("py_static")
_pkg.__path__ = [str(REPO / "py_static")]
sys.modules["py_static"] = _pkg
from py_static.state import State  # noqa: E402
from py_static.geometry import Geometry  # noqa: E402

# kurt_comp's standard dx=0.02 domain (6c x 3c)
DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.5)
RE = 1000.0
R_LE_0012 = 1.1019 * 0.12 ** 2  # =~ 0.01587, NACA0012 LE radius of curvature (chord=1)


def grid_xy(dx, domain=DOMAIN):
    """Cell-center-ish coordinate arrays matching State.omega's shape,
    same convention used throughout kurt_comp's gen_kurt_figs.py."""
    nx = int(round(domain["length"] / dx))
    ny = int(round(3.0 / dx))
    xs = domain["xoffset"] + np.arange(1, nx) * dx
    ys = domain["yoffset"] + np.arange(1, ny) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return X, Y, dx


def load_omega(run_dir, step):
    """omega field (2D array, x-major) from a flowNNNNN.bin restart file."""
    return State(filename=str(pathlib.Path(run_dir) / f"flow{step:05d}.bin")).omega._data[0].copy()


def load_state(run_dir, step):
    return State(filename=str(pathlib.Path(run_dir) / f"flow{step:05d}.bin"))


def load_geom_points(geom_path):
    """Boundary point (x,y) coords + arc length s (cumulative, point 0 at s=0),
    plus perimeter and the index/arc-length of the LE (min-x) and TE (max-x)
    points. Points come straight from the .geom file -- fixed in the solver's
    body frame regardless of alpha (this solver imposes AoA by rotating the
    free-stream, not the body; see 1-paper_based/README.md's "Wake vorticity
    fields" section)."""
    geom = Geometry(str(geom_path))
    n = geom.getNumPoints()
    pts = geom.getPoints()
    x = pts._data[0 * n:1 * n].copy()
    y = pts._data[1 * n:2 * n].copy()
    seg = np.hypot(np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0])))
    s = np.concatenate([[0.0], np.cumsum(seg)])[:-1]
    perimeter = s[-1] + seg[-1]
    i_le = int(np.argmin(x))
    i_te = int(np.argmax(x))
    return dict(x=x, y=y, s=s, n=n, perimeter=perimeter,
                i_le=i_le, i_te=i_te, s_le=s[i_le], s_te=s[i_te])


def nearest_index(arr, value):
    return int(np.argmin(np.abs(arr - value)))


def window(X, Y, field, xlim, ylim):
    """Return the sub-array of `field` (and the matching X,Y) within a
    physical (x,y) box -- avoids needing to hand-compute grid indices
    every time a test wants to zoom on the LE or TE."""
    xs = X[:, 0]
    ys = Y[0, :]
    ix = np.where((xs >= xlim[0]) & (xs <= xlim[1]))[0]
    iy = np.where((ys >= ylim[0]) & (ys <= ylim[1]))[0]
    sub = field[np.ix_(ix, iy)]
    return X[np.ix_(ix, iy)], Y[np.ix_(ix, iy)], sub, ix, iy


def dat_to_raw_dx(dat_pts, dx, out_txt, spacing_factor=1.0):
    """Resample a closed (x,y) polyline to uniform arc-length spacing
    dx*spacing_factor and write an IBPM 'raw' point file -- thin wrapper
    around SURF_test/make_airfoil_raw.py's resample_uniform/write_raw so
    this folder doesn't duplicate that logic."""
    sys.path.insert(0, str(REPO / "SURF_test"))
    from make_airfoil_raw import resample_uniform, write_raw  # noqa: E402
    new_pts, perimeter = resample_uniform(dat_pts, dx * spacing_factor)
    write_raw(new_pts, out_txt)
    return len(new_pts), perimeter


def write_geom(out_geom, raw_path, center=(0.25, 0.0)):
    pathlib.Path(out_geom).write_text(
        f"body body\n  raw {raw_path}\n  center {center[0]} {center[1]}\nend\n"
    )


def run_case(impl, geom, outdir, nx, ny, dt, nsteps, alpha=0.0, re=RE,
             domain=DOMAIN, restart=0, force_every=1):
    """Launch one py_static or cpp_static case -- same CLI pattern as
    ../1-paper_based/run_kurt_suite.py's run_one(), factored out here so
    every new run in this folder goes through one place."""
    import subprocess
    import time
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    common = [
        "-geom", str(geom), "-name", "flow", "-outdir", str(outdir),
        "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
        "-length", str(domain["length"]), "-xoffset", str(domain["xoffset"]),
        "-yoffset", str(domain["yoffset"]), "-alpha", str(alpha), "-Re", str(re),
        "-dt", str(dt), "-nsteps", str(nsteps),
        "-tecplot", "0", "-restart", str(restart), "-force", str(force_every),
    ]
    if impl == "cpp":
        cmd = [str(CPP_BIN)] + common
    else:
        cmd = [sys.executable, "-u", str(PY_RUNNER)] + common
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    ok = proc.returncode == 0 and (outdir / f"flow{nsteps:05d}.bin").exists()
    return ok, elapsed


def is_done(outdir, nsteps):
    return (pathlib.Path(outdir) / f"flow{nsteps:05d}.bin").exists()
