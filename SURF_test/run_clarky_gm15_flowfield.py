"""Flowfield case for airfoils/ClarkY and airfoils/GM15 (Python), same
methodology/parameters as run_flowfield.py (SD7003/SD8000, now in
airfoils/): dx=0.02, dt=0.01, nsteps=3000 (t=30), restart every 250 steps
(7 snapshots at t=0,5,...,30).

Usage:  python3 SURF_test/run_airfoils_flowfield.py
Output: SURF_test/airfoils/{ClarkY,GM15}/_run_data/flowfield/
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from airfoil_driver import run_case

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

CASES = {
    "ClarkY": dict(geom=GEOMDIR / "clarky_dx0.0200.geom", Re=60700, alpha=4.16),
    "GM15": dict(geom=GEOMDIR / "gm15_dx0.0200.geom", Re=40600, alpha=4.61),
}

for name, cfg in CASES.items():
    outdir = REPO / "SURF_test" / "airfoils" / name / "_run_data" / "flowfield"
    fpath, elapsed = run_case(
        geom=cfg["geom"], name="flow", outdir=outdir, alpha=cfg["alpha"],
        nx=300, ny=150, Re=cfg["Re"], dt=0.01, nsteps=3000, restart=250, **DOMAIN,
    )
    print(name, "flowfield run done,", elapsed, "s")

print("done")
