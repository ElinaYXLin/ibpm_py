"""
make_geoms_6.py

Builds every non-standard geometry needed across this folder's tests:
reconciliation thread #2 (recon2), and Groups C/D/E of the Test-3b
thickness-trend-confound investigation (Group A is pure reanalysis of
../5-leading_edge's existing fields; Group B needs no new geometry, just
shifted grid xoffset/yoffset at run time; Group F is a static diagnostic
on already-existing geometries).

Usage: python3 make_geoms_6.py
"""
import sys

import numpy as np

import common as c

sys.path.insert(0, str(c.REPO / "SURF_test"))
from make_airfoil_raw import load_dat, resample_uniform, write_raw  # noqa: E402

DX = 0.02
R_LE = c.R_LE_0012


# ---------------------------------------------------------- shared helpers
def variable_ds_resample(pts, dx_base, s_anchors, window, ds_at_anchor, n_fine=20000):
    """Same nearest-anchor variable-density resampler as
    ../5-leading_edge/make_geoms.py (fixed there to handle sparse AND
    dense correctly) -- duplicated here so this folder doesn't depend on
    importing a sibling test folder's script."""
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    seg = np.diff(pts, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    s = np.concatenate([[0.0], np.cumsum(seglen)])
    perimeter = s[-1]
    s_fine = np.linspace(0.0, perimeter, n_fine)

    def ds_of_s(ss):
        best_ds = np.full_like(ss, dx_base)
        best_d = np.full(ss.shape, np.inf)
        for s_anchor in s_anchors:
            d = np.abs(ss - s_anchor)
            d = np.minimum(d, perimeter - d)
            frac = np.clip(d / window, 0.0, 1.0)
            ramp = 0.5 * (1 - np.cos(np.pi * frac))
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


def write_variant(name, new_pts):
    out_txt = c.GEOMDIR / f"{name}.txt"
    out_geom = c.GEOMDIR / f"{name}.geom"
    write_raw(new_pts, str(out_txt))
    c.write_geom(out_geom, out_txt)
    print(f"wrote {out_geom.name}: {len(new_pts)} points")
    return out_geom


# ------------------------------------------------------------- recon2: LE-only dense
def make_le_only_dense():
    """Matches ../2-leading_edge_investigation/make_le_densified_geom.py's
    exact recipe (R_LE window, factor 4) but as a standalone .txt/.geom
    here, so recon2 can run it fresh at Re=1000 (the older investigation
    only ever ran this shape at Re=500)."""
    pts = load_dat(str(c.NACA0012_DAT))
    pts_c, s, s_le, s_te = find_s_le_te(pts)
    new_pts, perimeter, s_new = variable_ds_resample(pts_c, DX, [s_le], 2.0 * R_LE, DX / 4.0)
    write_variant("naca0012_dx0.0200_LEonly_dense", new_pts)


# ------------------------------------------------------- Group C1: per-shape refinement
def make_thickness_variant_at_dx(name, thickness_ratio, dx):
    pts = load_dat(str(c.NACA0012_DAT))
    pts_scaled = np.column_stack([pts[:, 0], pts[:, 1] * thickness_ratio])
    new_pts, perimeter = resample_uniform(pts_scaled, dx)
    write_variant(f"{name}_dx{dx:.4f}", new_pts)


# ------------------------------------------------- Group C2: thickness family @ dx=0.02
THICKNESS_FAMILY = {
    "naca0004": 4 / 12, "naca0008": 8 / 12, "naca0010": 10 / 12,
    "naca0014": 14 / 12, "naca0016": 16 / 12, "naca0020": 20 / 12,
}


# ------------------------------------------------------ Group D1: LE point density levels
def make_le_density_level(factor):
    """factor = dx/ds at the LE (>1 denser, <1 sparser); factor=4 matches
    recon2's LEonly_dense exactly (not regenerated twice)."""
    pts = load_dat(str(c.NACA0012_DAT))
    pts_c, s, s_le, s_te = find_s_le_te(pts)
    new_pts, perimeter, s_new = variable_ds_resample(pts_c, DX, [s_le], 2.0 * R_LE, DX / factor)
    write_variant(f"naca0012_dx0.0200_LEdensity{factor:g}x", new_pts)


# ------------------------------------------- Group E1: decouple nose vs tail thickness
def blended_thickness(x, y, x0, x1, ratio):
    """Scale |y| by `ratio` for x<=x0, ramp (half-cosine) back to 1.0 by
    x1, unchanged (ratio=1) for x>=x1."""
    y2 = y.copy()
    frac = np.clip((x - x0) / (x1 - x0), 0.0, 1.0)
    ramp = 0.5 * (1 - np.cos(np.pi * frac))  # 0 at x0 -> 1 at x1
    local_ratio = ratio + ramp * (1.0 - ratio)
    mask = x <= x1
    y2[mask] = y[mask] * local_ratio[mask]
    return y2


def make_front_only_variant(name, ratio):
    pts = load_dat(str(c.NACA0012_DAT))
    x, y = pts[:, 0], pts[:, 1]
    y2 = blended_thickness(x, y, x0=0.2, x1=0.4, ratio=ratio)
    new_pts, perimeter = resample_uniform(np.column_stack([x, y2]), DX)
    write_variant(name, new_pts)


def make_te_only_variant(name, ratio):
    """Ramp from unchanged (ratio=1) at x=0.6 to `ratio` at the TE (x=1),
    holding the front (x<0.6) at native NACA0012 thickness."""
    pts = load_dat(str(c.NACA0012_DAT))
    x, y = pts[:, 0], pts[:, 1]
    frac = np.clip((x - 0.6) / (1.0 - 0.6), 0.0, 1.0)
    ramp = 0.5 * (1 - np.cos(np.pi * frac))
    local_ratio = 1.0 + ramp * (ratio - 1.0)
    y2 = y * local_ratio
    new_pts, perimeter = resample_uniform(np.column_stack([x, y2]), DX)
    write_variant(name, new_pts)


# --------------------------------------------------- Group E2: common-TE resweep
def make_round_te(x, y, x_blend=0.95, min_te_half_thick=0.02):
    y2 = y.copy()
    mask = x >= x_blend
    frac = (x[mask] - x_blend) / (1 - x_blend)
    ramp = 0.5 * (1 - np.cos(np.pi * frac))
    floor = min_te_half_thick * ramp
    sign = np.sign(y[mask])
    sign[sign == 0] = 1.0
    y2[mask] = sign * np.maximum(np.abs(y[mask]), floor)
    return y2


def make_thickness_with_round_te(name, thickness_ratio):
    pts = load_dat(str(c.NACA0012_DAT))
    x, y = pts[:, 0], pts[:, 1] * thickness_ratio
    y2 = make_round_te(x, y)
    new_pts, perimeter = resample_uniform(np.column_stack([x, y2]), DX)
    write_variant(name, new_pts)


def main():
    print("--- recon2: LE-only dense (Re=1000 companion to the Re=500 investigation's shape) ---")
    make_le_only_dense()

    print("\n--- Group C1: naca0006/naca0018 at dx=0.01, dx=0.005 ---")
    for name, ratio in [("naca0006", 0.5), ("naca0018", 1.5)]:
        for dx in (0.01, 0.005):
            make_thickness_variant_at_dx(name, ratio, dx)

    print("\n--- Group C2: thickness family at dx=0.02 ---")
    for name, ratio in THICKNESS_FAMILY.items():
        make_thickness_variant_at_dx(name, ratio, DX)

    print("\n--- Group D1: LE point-density levels (factor=4 reuses recon2's LEonly_dense) ---")
    for factor in (0.5, 2, 8, 16):
        make_le_density_level(factor)

    print("\n--- Group E1: front-only / TE-only thickness variants ---")
    make_front_only_variant("naca0012_dx0.0200_frontsharp_TEnative", 0.5)
    make_front_only_variant("naca0012_dx0.0200_frontblunt_TEnative", 1.5)
    make_te_only_variant("naca0012_dx0.0200_TEsharp_frontnative", 0.5)

    print("\n--- Group E2: common-TE resweep (naca0006/0018 with the same rounded TE) ---")
    make_thickness_with_round_te("naca0006_dx0.0200_roundTE", 0.5)
    make_thickness_with_round_te("naca0018_dx0.0200_roundTE", 1.5)


if __name__ == "__main__":
    main()
