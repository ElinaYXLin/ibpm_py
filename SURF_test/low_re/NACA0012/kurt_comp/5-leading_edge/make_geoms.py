"""
make_geoms.py

Builds every non-standard geometry variant needed for Group 3 (Test 3a:
boundary-point spacing; Test 3b: curvature/bluntness sweep) of the LE/TE
striping investigation. All variants keep the background Eulerian grid at
dx=0.02 (Test 3a) or are resampled at dx=0.02 (Test 3b) -- only the
boundary POINTS or the airfoil SHAPE change, one variable at a time.

Test 3a variants (extends ../../2-leading_edge_investigation's LE-only,
denser-only densification to both LE+TE, and adds a sparser counterpart):
  - naca0012_dx0.0200_LTEdense.geom  -- ds -> dx/4 at both LE and TE
  - naca0012_dx0.0200_LTEsparse.geom -- ds -> 4*dx at both LE and TE

Test 3b variants:
  - naca0006_dx0.0200.geom / naca0018_dx0.0200.geom -- sharper/blunter nose
    (NACA00xx thickness scales linearly with the 2-digit thickness, so
    these are naca0012.dat.txt's own y-column scaled by 0.06/0.12 and
    0.18/0.12 respectively -- same x-distribution/style, no new digitizing)
  - naca0012_dx0.0200_roundTE.geom -- baseline NACA0012 with the last 5%
    chord's thickness floored (half-cosine blended) to a blunt 2%c
    half-thickness at the trailing edge, isolating whether the TE blobs
    are a sharp/near-cusped-TE artifact
  - cylinder anchor: reuses ../../../vortall/3-grid_refine/geom's existing
    cylinder_dx0.0200.geom (diameter=1, procedurally generated via
    `circle_n`, no new file needed) -- the "very blunt, constant curvature"
    end of the sweep.

Usage: python3 make_geoms.py
"""
import sys

import numpy as np

import common as c

sys.path.insert(0, str(c.REPO / "SURF_test"))
from make_airfoil_raw import load_dat, resample_uniform, write_raw  # noqa: E402

DX = 0.02
R_LE = c.R_LE_0012  # =~0.01587, NACA0012 LE radius of curvature (chord=1)


# ---------------------------------------------------------------- Test 3a
def variable_ds_resample(pts, dx_base, s_anchors, window, ds_at_anchor, n_fine=20000):
    """Generalizes ../../2-leading_edge_investigation/make_le_densified_geom.py's
    variable-density resampler to an arbitrary set of anchor arc-lengths
    (here: LE and TE), each with its own local target spacing ds_at_anchor
    (smaller = denser, larger = sparser than dx_base), half-cosine blended
    back to dx_base over `window` arc length on each side."""
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    seg = np.diff(pts, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    s = np.concatenate([[0.0], np.cumsum(seglen)])
    perimeter = s[-1]
    s_fine = np.linspace(0.0, perimeter, n_fine)

    def ds_of_s(ss):
        # each point uses whichever anchor is nearest (LE/TE windows don't
        # overlap for a 2*R_LE-scale window), applying THAT anchor's own
        # half-cosine ramp -- unlike a min()/max() combination rule, this
        # works correctly whether ds_at_anchor is smaller (denser) or
        # larger (sparser) than dx_base
        best_ds = np.full_like(ss, dx_base)
        best_d = np.full(ss.shape, np.inf)
        for s_anchor in s_anchors:
            d = np.abs(ss - s_anchor)
            d = np.minimum(d, perimeter - d)
            frac = np.clip(d / window, 0.0, 1.0)
            ramp = 0.5 * (1 - np.cos(np.pi * frac))  # 0 at d=0 (anchor), 1 at d>=window
            ds_here = ds_at_anchor + ramp * (dx_base - ds_at_anchor)
            closer = d < best_d
            best_ds = np.where(closer, ds_here, best_ds)
            best_d = np.where(closer, d, best_d)
        return best_ds

    ds_fine = ds_of_s(s_fine)
    density = 1.0 / ds_fine
    phase = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(s_fine))])
    n_points = int(round(phase[-1]))
    phase_targets = np.arange(n_points) * (phase[-1] / n_points)
    s_new = np.interp(phase_targets, phase, s_fine)
    x_new = np.interp(s_new, s, pts[:, 0])
    y_new = np.interp(s_new, s, pts[:, 1])
    return np.column_stack([x_new, y_new]), perimeter, s_new


def find_s_le_te(pts):
    if not np.allclose(pts[0], pts[-1]):
        pts_c = np.vstack([pts, pts[0]])
    else:
        pts_c = pts
    seg = np.diff(pts_c, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    s = np.concatenate([[0.0], np.cumsum(seglen)])
    i_le = int(np.argmin(pts_c[:, 0]))
    i_te = int(np.argmax(pts_c[:, 0]))
    return pts_c, s, s[i_le], s[i_te]


def make_spacing_variant(name, ds_at_anchor):
    pts = load_dat(str(c.NACA0012_DAT))
    pts_c, s, s_le, s_te = find_s_le_te(pts)
    new_pts, perimeter, s_new = variable_ds_resample(pts_c, DX, [s_le, s_te], 2.0 * R_LE, ds_at_anchor)
    out_txt = c.GEOMDIR / f"naca0012_dx0.0200_{name}.txt"
    out_geom = c.GEOMDIR / f"naca0012_dx0.0200_{name}.geom"
    write_raw(new_pts, str(out_txt))
    c.write_geom(out_geom, out_txt)
    print(f"wrote {out_geom.name}: {len(new_pts)} points (ds_at_LE/TE={ds_at_anchor:.4f}, "
          f"background dx={DX}, ds_at_LE_ratio={ds_at_anchor / DX:.2f}x)")


# ---------------------------------------------------------------- Test 3b
def make_thickness_variant(name, thickness_ratio):
    """NACA00xx: the thickness distribution y_t(x) scales linearly with
    the 2-digit thickness parameter, so scaling naca0012's own y-column by
    thickness_ratio (0.5 for NACA0006, 1.5 for NACA0018) reproduces the
    exact same shape family without needing a new .dat source."""
    pts = load_dat(str(c.NACA0012_DAT))
    pts_scaled = np.column_stack([pts[:, 0], pts[:, 1] * thickness_ratio])
    new_pts, perimeter = resample_uniform(pts_scaled, DX)
    out_txt = c.GEOMDIR / f"{name}_dx0.0200.txt"
    out_geom = c.GEOMDIR / f"{name}_dx0.0200.geom"
    write_raw(new_pts, str(out_txt))
    c.write_geom(out_geom, out_txt)
    r_le = R_LE * thickness_ratio ** 2  # r_LE ~ t^2 for NACA00xx
    print(f"wrote {out_geom.name}: {len(new_pts)} points, r_LE~{r_le:.5f}c "
          f"(vs NACA0012's {R_LE:.5f}c)")


def make_round_te_variant(x_blend=0.95, min_te_half_thick=0.02):
    pts = load_dat(str(c.NACA0012_DAT))
    x, y = pts[:, 0], pts[:, 1]
    y2 = y.copy()
    mask = x >= x_blend
    frac = (x[mask] - x_blend) / (1 - x_blend)
    ramp = 0.5 * (1 - np.cos(np.pi * frac))  # 0 at x_blend -> 1 at TE
    floor = min_te_half_thick * ramp
    sign = np.sign(y[mask])
    sign[sign == 0] = 1.0
    y2[mask] = sign * np.maximum(np.abs(y[mask]), floor)
    new_pts, perimeter = resample_uniform(np.column_stack([x, y2]), DX)
    out_txt = c.GEOMDIR / "naca0012_dx0.0200_roundTE.txt"
    out_geom = c.GEOMDIR / "naca0012_dx0.0200_roundTE.geom"
    write_raw(new_pts, str(out_txt))
    c.write_geom(out_geom, out_txt)
    print(f"wrote {out_geom.name}: {len(new_pts)} points, TE half-thickness floored to "
          f"{min_te_half_thick}c (blend from x={x_blend}c), vs NACA0012's native ~0.0013c")


def main():
    make_spacing_variant("LTEdense", DX / 4.0)
    make_spacing_variant("LTEsparse", DX * 4.0)
    make_thickness_variant("naca0006", 0.5)
    make_thickness_variant("naca0018", 1.5)
    make_round_te_variant()
    print("cylinder anchor: reusing "
          f"{(c.REPO / 'SURF_test' / 'vortall' / '3-grid_refine' / 'geom' / 'cylinder_dx0.0200.geom').relative_to(c.REPO)} "
          "(diameter=1, no new file needed)")


if __name__ == "__main__":
    main()
