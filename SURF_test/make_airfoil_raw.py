"""Resample a Selig-format airfoil .dat file to uniform arc-length spacing
matched to a given grid dx, and write an IBPM 'raw' point file.

The raw UIUC coordinate files have highly non-uniform point spacing
(much finer than dx near the leading edge, coarser mid-chord); feeding
them directly into py/ibpm.py causes the projection matrix to be
singular ('over-resolved boundary' -- points spaced far closer together
than the grid, a documented degeneracy in this codebase, see
py/cholesky_solver.py's module docstring / SURF_test/built_in_tests/README.md).
Standard
practice for immersed-boundary methods is to keep the boundary point
spacing close to dx (roughly 0.8-1.2 dx); this module re-parametrizes
the airfoil boundary by arc length and resamples at that spacing.
"""
import numpy as np
import pathlib


def load_dat(path):
    lines = pathlib.Path(path).read_text().splitlines()
    pts = []
    for l in lines[1:]:
        l = l.strip()
        if not l:
            continue
        x, y = l.split()
        pts.append((float(x), float(y)))
    return np.array(pts)


def resample_uniform(pts, ds):
    """Resample a closed polyline `pts` (N,2) to uniform arc-length
    spacing `ds`, returning a new (M,2) array (closed: last==first)."""
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    seg = np.diff(pts, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    s = np.concatenate([[0.0], np.cumsum(seglen)])
    perimeter = s[-1]
    n = max(int(round(perimeter / ds)), 8)
    s_new = np.linspace(0.0, perimeter, n, endpoint=False)
    x_new = np.interp(s_new, s, pts[:, 0])
    y_new = np.interp(s_new, s, pts[:, 1])
    return np.column_stack([x_new, y_new]), perimeter


def write_raw(pts, out_path):
    with open(out_path, "w") as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x:.6f} {y:.6f}\n")


def make_raw_for_dx(dat_path, dx, out_path, spacing_factor=1.0):
    pts = load_dat(dat_path)
    new_pts, perimeter = resample_uniform(pts, dx * spacing_factor)
    write_raw(new_pts, out_path)
    return len(new_pts), perimeter


if __name__ == "__main__":
    import sys
    dat_path, dx, out_path = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    n, p = make_raw_for_dx(dat_path, dx, out_path)
    print(f"wrote {out_path}: {n} points, perimeter={p:.4f}, ds={p/n:.5f} (target dx={dx})")
