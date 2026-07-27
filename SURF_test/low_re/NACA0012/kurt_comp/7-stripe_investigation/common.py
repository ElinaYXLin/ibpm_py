"""
common.py

Shared helpers for kurt_comp/7-stripe_investigation -- the mentor's
follow-up on `../6-edges_further`'s Group G finding (every LE/TE peak
lands in the fluid, not inside the body): quantify the striping that
extends *upstream* of the leading edge specifically (its area and its
severity), using the region x < x_LE as a null region where the true
physical vorticity is ~0 (diffusion length scale nu/U = c/Re ~ 0.001c at
Re=1000, ~20x smaller than a single dx=0.02 cell), so anything measured
there is numerical artifact with no real-physics signal mixed in.

Same conventions as ../6-edges_further/common.py (independently
redefined here rather than imported, matching that folder's own relation
to ../5-leading_edge/common.py): py_static (never modified) via its
internal API for analysis, py_static + cpp_static compared everywhere
data for both exists.
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

HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "7-stripe_investigation"
DATA = HERE / "data"
FIGS = HERE / "figures"
RUNS = HERE / "runs"
for d in (DATA, FIGS, RUNS):
    d.mkdir(parents=True, exist_ok=True)

KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp"
KURT1 = KURT / "1-paper_based"
KURT5 = KURT / "5-leading_edge"
KURT6 = KURT / "6-edges_further"
BASE_GEOM_DX002 = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
CYLINDER_GEOM_DX002 = REPO / "SURF_test" / "vortall" / "3-grid_refine" / "geom" / "cylinder_dx0.0200.geom"

CPP_BIN = REPO / "build_static" / "ibpm"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"

sys.path.insert(0, str(REPO))
_pkg = types.ModuleType("py_static")
_pkg.__path__ = [str(REPO / "py_static")]
sys.modules["py_static"] = _pkg
from py_static.state import State  # noqa: E402
from py_static.geometry import Geometry  # noqa: E402
from py_static.grid import Grid  # noqa: E402
from py_static.navier_stokes_model import NavierStokesModel  # noqa: E402
from py_static.scalar import Scalar  # noqa: E402

IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}

# standard kurt_comp dx=0.02 domain (6c x 3c), used by every source run this
# folder reuses (1-paper_based, 5-leading_edge Groups 2 and 3)
DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.5, yheight=3.0)
RE = 1000.0


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


def is_done(outdir, nsteps):
    return (pathlib.Path(outdir) / f"flow{nsteps:05d}.bin").exists()


def run_case(impl, geom, outdir, nx, ny, dt, nsteps, alpha=0.0, re=RE,
             domain=DOMAIN, restart=0, force_every=1, ngrid=1):
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    common = [
        "-geom", str(geom), "-name", "flow", "-outdir", str(outdir),
        "-nx", str(nx), "-ny", str(ny), "-ngrid", str(ngrid),
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


# --------------------------------------------------------------------------
# Upstream-region metrics (x < x_le - buffer). LE tip is at x~0 for every
# geometry used here (raw airfoil/cylinder coords span x in [0,1] or
# [-0.5,0.5]-recentered; callers pass the actual x_le from load_geom_points).
# --------------------------------------------------------------------------

def upstream_mask(X, x_le, buffer_dx, dx, x_domain_min=None, margin_cells=2):
    """Boolean mask for the upstream null region: x in
    [x_domain_min + margin, x_le - buffer_dx*dx], all y. `buffer_dx` is the
    exclusion width in units of dx (keeps genuine LE-wraparound boundary
    layer out); `margin_cells` similarly excludes a few cells next to the
    domain's own upstream edge."""
    xs = X[:, 0]
    lo = xs.min() + margin_cells * dx if x_domain_min is None else x_domain_min + margin_cells * dx
    hi = x_le - buffer_dx * dx
    return (xs >= lo) & (xs <= hi)


def upstream_profile(X, Y, om, x_le, dx, buffer_dx=2.0, ylim=(-0.5, 0.5)):
    """Per-x-column stats over y in ylim, for x in the upstream region.
    Returns dict of 1-D arrays (xs, max_abs, rms, int_abs_dy, int_signed_dy)."""
    m = upstream_mask(X, x_le, buffer_dx, dx)
    ys = Y[0, :]
    iy = np.where((ys >= ylim[0]) & (ys <= ylim[1]))[0]
    xs = X[m, 0]
    sub = om[np.ix_(np.where(m)[0], iy)]  # (nx_up, ny_win)
    max_abs = np.abs(sub).max(axis=1)
    rms = np.sqrt((sub ** 2).mean(axis=1))
    int_abs = np.trapz(np.abs(sub), ys[iy], axis=1)
    int_signed = np.trapz(sub, ys[iy], axis=1)
    return dict(xs=xs, max_abs=max_abs, rms=rms, int_abs_dy=int_abs, int_signed_dy=int_signed)


def upstream_scalar_metrics(X, Y, om, x_le, dx, buffer_dx=2.0, ylim=(-0.5, 0.5), thresholds=(1.0, 2.0, 4.0, 8.0)):
    """Single-number summary over the whole upstream window: total
    enstrophy (0.5*int om^2 dA), total |om| integral, signed integral,
    peak |om|, and threshold-area A(tau) for each tau in `thresholds`."""
    m = upstream_mask(X, x_le, buffer_dx, dx)
    ys = Y[0, :]
    iy = np.where((ys >= ylim[0]) & (ys <= ylim[1]))[0]
    sub = om[np.ix_(np.where(m)[0], iy)]
    dA = dx * dx
    enstrophy = 0.5 * float(np.sum(sub ** 2) * dA)
    int_abs = float(np.sum(np.abs(sub)) * dA)
    int_signed = float(np.sum(sub) * dA)
    peak = float(np.abs(sub).max()) if sub.size else float("nan")
    areas = {tau: float(np.sum(np.abs(sub) > tau) * dA) for tau in thresholds}
    n_up_cols = int(m.sum())
    window_extent_chord = float(x_le - X[m, 0].min()) if n_up_cols else 0.0
    return dict(enstrophy=enstrophy, int_abs=int_abs, int_signed=int_signed,
                peak=peak, areas=areas, window_extent_chord=window_extent_chord,
                window_extent_cells=window_extent_chord / dx)


def reach_L_up(xs, max_abs, x_le, noise_floor):
    """Given a per-column max|omega| profile (from upstream_profile,
    ordered by x, values only defined for x < x_le), find L_up: the
    distance upstream of x_le out to the *last* column (walking outward
    from the LE) that still exceeds noise_floor, i.e. the outermost point
    of the contiguous exceedance run touching x_le. Returns
    (L_up_chord, L_up_cells) -- cells needs dx supplied by the caller
    (L_up_chord / dx), kept as a separate division so callers control
    rounding. Returns (0.0, 0.0) if the column nearest x_le doesn't even
    exceed the floor."""
    order = np.argsort(-xs)  # descending x: nearest LE first
    xs_o = xs[order]
    vals_o = max_abs[order]
    exceeds = vals_o > noise_floor
    if not exceeds.size or not exceeds[0]:
        return 0.0, xs_o, exceeds
    # first index (from the LE side) where it drops below floor
    below = np.where(~exceeds)[0]
    last_idx = below[0] - 1 if below.size else len(exceeds) - 1
    L_up_chord = float(x_le - xs_o[last_idx])
    return L_up_chord, xs_o, exceeds
