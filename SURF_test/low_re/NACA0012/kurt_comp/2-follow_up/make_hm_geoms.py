"""
make_hm_geoms.py

Builds the one non-standard geometry Tests H-M need: a pitching (f=4Hz,
alpha0=0) variant of ../5-leading_edge/geom/naca0012_dx0.0200_LTEdense.geom
(boundary points refined to ds=dx/4 at the LE+TE, background grid
unchanged at dx=0.02). ../5-leading_edge only ever ran this geometry
STEADY; Test J needs the dynamic (pitching) case, built the same way
../1-paper_based/run_kurt_suite.py's ensure_geom() builds pitching variants
from a steady raw point file.

Usage: python3 make_hm_geoms.py
Output: geom/naca0012_dx0.0200_LTEdense_f4hz.geom
"""
import pathlib

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "2-follow_up"
GEOMDIR = HERE / "geom"
GEOMDIR.mkdir(exist_ok=True)

LTEDENSE_RAW = (REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" /
                "5-leading_edge" / "geom" / "naca0012_dx0.0200_LTEdense.txt")

PITCH_AMP = 0.0174532925199433  # 1 deg in radians
PITCH_PHASE = 3.14159265358979  # pi -> effective AoA = alpha0 + A sin(2 pi f t)
FREQ_F4HZ = 2.73972602739726


def main():
    out_geom = GEOMDIR / "naca0012_dx0.0200_LTEdense_f4hz.geom"
    out_geom.write_text(
        "body body\n"
        f"  raw {LTEDENSE_RAW}\n"
        "  center 0.25 0.0\n"
        f"  motion pitchplunge {PITCH_AMP:.16g} {FREQ_F4HZ:.16g} {PITCH_PHASE:.16g} 0 0 0\n"
        "end\n"
    )
    print(f"wrote {out_geom}")


if __name__ == "__main__":
    main()
