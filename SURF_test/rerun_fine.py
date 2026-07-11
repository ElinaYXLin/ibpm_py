import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from airfoil_driver import run_case, time_avg_force
from make_airfoil_raw import make_raw_for_dx

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
AVG_FRAC = 0.6

# fine level diverged at dt=0.01; halve dt, double nsteps to reach the same t=30
DT = 0.005
NSTEPS = 6000

CASES = {
    "SD7003": dict(dat=REPO / "SURF_test" / "airfoils" / "SD7003" / "sd7003.dat.txt", Re=61100, conv_alpha=-0.09),
    "SD8000": dict(dat=REPO / "SURF_test" / "airfoils" / "SD8000" / "sd8000.dat.txt", Re=60800, conv_alpha=-0.81),
}
FINE = dict(tag="fine", dx=0.01, nx=600, ny=300)

results_path = REPO / "SURF_test" / "batch_results.json"
results = json.loads(results_path.read_text())

for name, cfg in CASES.items():
    raw_path = GEOMDIR / f"{name.lower()}_dx{FINE['dx']:.4f}.txt"
    n, perim = make_raw_for_dx(cfg["dat"], FINE["dx"], raw_path)
    geom_path = GEOMDIR / f"{name.lower()}_dx{FINE['dx']:.4f}.geom"
    geom_path.write_text(f"body {name}\n  raw {raw_path}\n  center 0.25 0.0\nend\n")

    outdir = REPO / "SURF_test" / "airfoils" / name / "_run_data" / "conv_fine_dt0005"
    fpath, elapsed = run_case(
        geom=geom_path, name="run", outdir=outdir, alpha=cfg["conv_alpha"],
        nx=FINE["nx"], ny=FINE["ny"], Re=cfg["Re"], dt=DT, nsteps=NSTEPS, **DOMAIN,
    )
    stats = time_avg_force(fpath, frac=AVG_FRAC)
    stats.update(alpha=cfg["conv_alpha"], elapsed=elapsed, nx=FINE["nx"], ny=FINE["ny"],
                 dx=FINE["dx"], tag=FINE["tag"], npts=n, dt=DT)
    print(name, "fine (dt=0.005):", stats)

    # replace the NaN fine entry in-place
    conv_list = results["convergence"][name]
    conv_list[:] = [c for c in conv_list if c["tag"] != "fine"] + [stats]
    results_path.write_text(json.dumps(results, indent=2))

print("done")
