import sys, json, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from airfoil_driver import run_case, time_avg_force
from make_airfoil_raw import make_raw_for_dx

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
GEOMDIR.mkdir(exist_ok=True)

DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
DT = 0.01
NSTEPS = 3000
AVG_FRAC = 0.6  # average over the last 60% of the run

CASES = {
    "SD7003": dict(
        dat=REPO / "SURF_test" / "airfoils" / "LSAT-SD7003" / "sd7003.dat.txt",
        Re=61100,
        polar_alphas=[-2.92, -0.09, 1.66, 4.60, 7.72],
        conv_alpha=-0.09,
    ),
    "SD8000": dict(
        dat=REPO / "SURF_test" / "airfoils" / "LSAT-SD8000" / "sd8000.dat.txt",
        Re=60800,
        polar_alphas=[-3.88, -0.81, 2.29, 5.36, 8.36],
        conv_alpha=-0.81,
    ),
}

CONV_LEVELS = [
    dict(tag="coarse", dx=0.04, nx=150, ny=75),
    dict(tag="medium", dx=0.02, nx=300, ny=150),
    dict(tag="fine",   dx=0.01, nx=600, ny=300),
]
POLAR_LEVEL = dict(dx=0.02, nx=300, ny=150)  # production resolution for the polar sweep


def geom_for_dx(name, dat, dx):
    raw_path = GEOMDIR / f"{name.lower()}_dx{dx:.4f}.txt"
    n, perim = make_raw_for_dx(dat, dx, raw_path)
    geom_path = GEOMDIR / f"{name.lower()}_dx{dx:.4f}.geom"
    geom_path.write_text(
        f"body {name}\n  raw {raw_path}\n  center 0.25 0.0\nend\n"
    )
    return geom_path, n, perim


def main():
    results = {"polar": {}, "convergence": {}}
    results_path = REPO / "SURF_test" / "batch_results.json"

    for name, cfg in CASES.items():
        print(f"=== {name}: polar sweep at dx={POLAR_LEVEL['dx']} ===", flush=True)
        geom_path, npts, perim = geom_for_dx(name, cfg["dat"], POLAR_LEVEL["dx"])
        print(f"  geometry: {npts} points, perimeter={perim:.4f}", flush=True)
        results["polar"][name] = []
        for alpha in cfg["polar_alphas"]:
            t0 = time.time()
            outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data" / f"polar_a{alpha:+.2f}"
            fpath, elapsed = run_case(
                geom=geom_path, name="run", outdir=outdir, alpha=alpha,
                nx=POLAR_LEVEL["nx"], ny=POLAR_LEVEL["ny"], Re=cfg["Re"],
                dt=DT, nsteps=NSTEPS, **DOMAIN,
            )
            stats = time_avg_force(fpath, frac=AVG_FRAC)
            stats.update(alpha=alpha, elapsed=elapsed, nx=POLAR_LEVEL["nx"], ny=POLAR_LEVEL["ny"])
            results["polar"][name].append(stats)
            print(f"  alpha={alpha:+.2f}  Cl={stats['cl_mean']:+.4f}±{stats['cl_std']:.4f}  "
                  f"Cd={stats['cd_mean']:+.4f}±{stats['cd_std']:.4f}  ({elapsed:.1f}s)", flush=True)
            results_path.write_text(json.dumps(results, indent=2))

        print(f"=== {name}: grid convergence at alpha={cfg['conv_alpha']} ===", flush=True)
        results["convergence"][name] = []
        for lvl in CONV_LEVELS:
            geom_path, npts, perim = geom_for_dx(name, cfg["dat"], lvl["dx"])
            outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data" / f"conv_{lvl['tag']}"
            fpath, elapsed = run_case(
                geom=geom_path, name="run", outdir=outdir, alpha=cfg["conv_alpha"],
                nx=lvl["nx"], ny=lvl["ny"], Re=cfg["Re"],
                dt=DT, nsteps=NSTEPS, **DOMAIN,
            )
            stats = time_avg_force(fpath, frac=AVG_FRAC)
            stats.update(alpha=cfg["conv_alpha"], elapsed=elapsed, nx=lvl["nx"], ny=lvl["ny"],
                         dx=lvl["dx"], tag=lvl["tag"], npts=npts)
            results["convergence"][name].append(stats)
            print(f"  {lvl['tag']:8s} dx={lvl['dx']:.3f} nx={lvl['nx']:4d}  "
                  f"Cl={stats['cl_mean']:+.4f}±{stats['cl_std']:.4f}  "
                  f"Cd={stats['cd_mean']:+.4f}±{stats['cd_std']:.4f}  ({elapsed:.1f}s)", flush=True)
            results_path.write_text(json.dumps(results, indent=2))

    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
