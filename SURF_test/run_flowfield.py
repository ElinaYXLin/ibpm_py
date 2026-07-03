import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from airfoil_driver import run_case

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

CASES = {
    "SD7003": dict(geom=GEOMDIR / "sd7003_dx0.0200.geom", Re=61100, alpha=4.60),
    "SD8000": dict(geom=GEOMDIR / "sd8000_dx0.0200.geom", Re=60800, alpha=5.36),
}

for name, cfg in CASES.items():
    outdir = REPO / "results" / name / "_run_data" / "flowfield"
    fpath, elapsed = run_case(
        geom=cfg["geom"], name="flow", outdir=outdir, alpha=cfg["alpha"],
        nx=300, ny=150, Re=cfg["Re"], dt=0.01, nsteps=3000, restart=250, **DOMAIN,
    )
    print(name, "flowfield run done,", elapsed, "s")

print("done")
