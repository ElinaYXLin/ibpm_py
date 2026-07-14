"""
make_le_densified_geom.py

Test 3 of the LE vorticity investigation (see ../README.md).

Builds a NACA0012 boundary-point file with LOCALLY finer arc-length
spacing near the leading edge, while every other design choice (in
particular the background Eulerian grid dx) stays exactly at the ../..
baseline (dx=0.02). This isolates Lagrangian boundary-point density from
Eulerian grid resolution: ../../run_gridconv.py's dx sweep changes BOTH
at once (finer dx always comes with make_airfoil_raw.py resampling the
boundary to match, at spacing_factor=1.0 i.e. ds=dx), so it can't tell
whether the LE vorticity speck is driven by the background grid being too
coarse to resolve the curvature, or by the regularized delta-function
support (tied to boundary-POINT spacing, not grid dx) being too coarse
relative to curvature -- related but distinguishable mechanisms.

Point-spacing rule: ds(s) = dx everywhere except within an arc-length
window of +-2*r_LE around the leading edge (s_LE = arc-length location of
the point with minimum x), where it ramps down (half-cosine transition)
to ds = dx/LE_REFINE_FACTOR at the LE itself. r_LE = 1.1019*(0.12)^2 =
0.01587 =~ 0.016 (chord=1) is NACA0012's leading-edge radius of curvature
(the length scale the earlier investigation flagged as the relevant one).

This uses variable-density resampling: define a "phase" function
P(s) = integral_0^s ds'/ds(s'), which counts how many boundary points
should have been placed by arc length s; points are then placed at every
integer phase value (so wherever ds(s) is small, phase advances slowly,
i.e., points are placed close together in s).

Usage: python3 SURF_test/low_re/NACA0012/leading_edge_investigation/make_le_densified_geom.py
Output: SURF_test/geom/naca0012_dx0.0200_LEdense.{geom,txt}
"""
import pathlib
import sys

import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO / "SURF_test"))
from make_airfoil_raw import load_dat, write_raw  # noqa: E402

DAT_PATH = REPO / "SURF_test" / "low_re" / "NACA0012" / "naca0012.dat.txt"
GEOM_DIR = REPO / "SURF_test" / "geom"
DX_BASE = 0.02
R_LE = 1.1019 * 0.12 ** 2  # =~ 0.01587, NACA0012 leading-edge radius of curvature (chord=1)
WINDOW = 2.0 * R_LE  # half-width of the refined arc-length window around the LE
LE_REFINE_FACTOR = 4.0  # ds at the LE itself = DX_BASE / LE_REFINE_FACTOR
N_FINE = 20000  # dense sampling of the original polyline for the phase integral


def variable_ds_resample(pts, dx_base, s_le, window, refine_factor, n_fine=N_FINE):
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    seg = np.diff(pts, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    s = np.concatenate([[0.0], np.cumsum(seglen)])
    perimeter = s[-1]

    s_fine = np.linspace(0.0, perimeter, n_fine)

    def ds_of_s(ss):
        d = np.abs(ss - s_le)
        d = np.minimum(d, perimeter - d)  # wrap-around distance on the closed loop
        ds_min = dx_base / refine_factor
        # half-cosine ramp from ds_min (at d=0) to dx_base (at d>=window)
        frac = np.clip(d / window, 0.0, 1.0)
        ramp = 0.5 * (1 - np.cos(np.pi * frac))  # 0 at d=0, 1 at d>=window
        return ds_min + ramp * (dx_base - ds_min)

    ds_fine = ds_of_s(s_fine)
    density = 1.0 / ds_fine
    # phase P(s) = integral_0^s density(s') ds'
    phase = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(s_fine))])
    n_points = int(round(phase[-1]))
    phase_targets = np.arange(n_points) * (phase[-1] / n_points)
    s_new = np.interp(phase_targets, phase, s_fine)

    x_new = np.interp(s_new, s, pts[:, 0])
    y_new = np.interp(s_new, s, pts[:, 1])
    return np.column_stack([x_new, y_new]), perimeter, s_new


def main():
    pts = load_dat(str(DAT_PATH))
    # find LE (min-x point) arc-length location on the dense original polyline
    if not np.allclose(pts[0], pts[-1]):
        pts_closed = np.vstack([pts, pts[0]])
    else:
        pts_closed = pts
    seg = np.diff(pts_closed, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    s_all = np.concatenate([[0.0], np.cumsum(seglen)])
    i_le = np.argmin(pts_closed[:, 0])
    s_le = s_all[i_le]

    new_pts, perimeter, s_new = variable_ds_resample(pts_closed, DX_BASE, s_le, WINDOW, LE_REFINE_FACTOR)

    out_txt = GEOM_DIR / "naca0012_dx0.0200_LEdense.txt"
    out_geom = GEOM_DIR / "naca0012_dx0.0200_LEdense.geom"
    write_raw(new_pts, str(out_txt))
    out_geom.write_text(
        f"body NACA0012\n  raw {out_txt.relative_to(REPO)}\n  center 0.25 0.0\nend\n"
    )

    d_near_le = np.abs(s_new - s_le)
    d_near_le = np.minimum(d_near_le, perimeter - d_near_le)
    n_in_window = int(np.sum(d_near_le < WINDOW))
    print(f"wrote {out_geom.relative_to(REPO)}: {len(new_pts)} points total "
          f"({n_in_window} within +-{WINDOW:.4f} of LE, vs "
          f"{int(round(2 * WINDOW / DX_BASE))} at uniform dx={DX_BASE} spacing)")
    print(f"LE at (x,y)=({pts_closed[i_le,0]:.4f},{pts_closed[i_le,1]:.4f}), s_le={s_le:.4f}, perimeter={perimeter:.4f}")


if __name__ == "__main__":
    main()
