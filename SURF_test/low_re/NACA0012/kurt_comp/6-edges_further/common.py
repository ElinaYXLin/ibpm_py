"""
common.py

Shared helpers for kurt_comp/6-edges_further -- follow-ups to
5-leading_edge's two flagged open discrepancies (with the older Re=500
../../2-leading_edge_investigation/) plus Test 3b's non-monotonic
thickness-trend confounds (Groups A-F). Same conventions as
../5-leading_edge/common.py (which this imports rather than duplicates):
py_static (never modified) via its internal API for analysis, py_static +
cpp_static subprocess launches for every new run.
"""
import os
import pathlib
import subprocess
import sys
import time
import types

import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
os.chdir(REPO)  # .geom "raw <path>" lines are repo-root-relative

HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "6-edges_further"
DATA = HERE / "data"
FIGS = HERE / "figures"
GEOMDIR = HERE / "geom"
RUNS = HERE / "runs"
for d in (DATA, FIGS, GEOMDIR, RUNS):
    d.mkdir(parents=True, exist_ok=True)

KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp"
KURT5 = KURT / "5-leading_edge"
PRIOR = REPO / "SURF_test" / "low_re" / "NACA0012" / "2-leading_edge_investigation"
NACA0012_DAT = REPO / "SURF_test" / "low_re" / "NACA0012" / "1-basics" / "naca0012.dat.txt"
BASE_GEOM_DX002 = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"

CPP_BIN = REPO / "build_static" / "ibpm"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"

sys.path.insert(0, str(REPO))
_pkg = types.ModuleType("py_static")
_pkg.__path__ = [str(REPO / "py_static")]
sys.modules["py_static"] = _pkg
from py_static.state import State  # noqa: E402
from py_static.geometry import Geometry  # noqa: E402

IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}
R_LE_0012 = 1.1019 * 0.12 ** 2  # =~ 0.01587c


def grid_xy(dx, length=6.0, xoffset=-2.0, yoffset=-1.5, yheight=3.0):
    nx = int(round(length / dx))
    ny = int(round(yheight / dx))
    xs = xoffset + np.arange(1, nx) * dx
    ys = yoffset + np.arange(1, ny) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return X, Y


def load_omega(run_dir, step):
    return State(filename=str(pathlib.Path(run_dir) / f"flow{step:05d}.bin")).omega._data[0].copy()


def load_state(run_dir, step):
    return State(filename=str(pathlib.Path(run_dir) / f"flow{step:05d}.bin"))


def load_geom_points(geom_path):
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
    xs = X[:, 0]
    ys = Y[0, :]
    ix = np.where((xs >= xlim[0]) & (xs <= xlim[1]))[0]
    iy = np.where((ys >= ylim[0]) & (ys <= ylim[1]))[0]
    sub = field[np.ix_(ix, iy)]
    return X[np.ix_(ix, iy)], Y[np.ix_(ix, iy)], sub, ix, iy


def write_geom(out_geom, raw_path, center=(0.25, 0.0)):
    pathlib.Path(out_geom).write_text(
        f"body body\n  raw {raw_path}\n  center {center[0]} {center[1]}\nend\n"
    )


def dat_to_raw_dx(dat_pts, dx, out_txt, spacing_factor=1.0):
    sys.path.insert(0, str(REPO / "SURF_test"))
    from make_airfoil_raw import resample_uniform, write_raw  # noqa: E402
    new_pts, perimeter = resample_uniform(dat_pts, dx * spacing_factor)
    write_raw(new_pts, out_txt)
    return len(new_pts), perimeter


def is_done(outdir, nsteps):
    return (pathlib.Path(outdir) / f"flow{nsteps:05d}.bin").exists()


def run_case(impl, geom, outdir, nx, ny, dt, nsteps, alpha=0.0, re=1000.0,
             domain=None, ngrid=1, restart=0, force_every=1):
    domain = domain or dict(length=6.0, xoffset=-2.0, yoffset=-1.5)
    outdir = pathlib.Path(outdir)
    if is_done(outdir, nsteps):
        return True, 0.0, "skip"
    outdir.mkdir(parents=True, exist_ok=True)
    common = [
        "-geom", str(geom), "-name", "flow", "-outdir", str(outdir),
        "-nx", str(nx), "-ny", str(ny), "-ngrid", str(ngrid),
        "-length", str(domain["length"]), "-xoffset", str(domain["xoffset"]),
        "-yoffset", str(domain["yoffset"]), "-alpha", str(alpha), "-Re", str(re),
        "-dt", str(dt), "-nsteps", str(nsteps),
        "-tecplot", "0", "-restart", str(restart), "-force", str(force_every),
    ]
    cmd = [str(CPP_BIN)] + common if impl == "cpp" else [sys.executable, "-u", str(PY_RUNNER)] + common
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    ok = proc.returncode == 0 and is_done(outdir, nsteps)
    return ok, elapsed, ("ok" if ok else "FAILED")
