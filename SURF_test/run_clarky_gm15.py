"""Driver: run py/ibpm.py on the airfoils/ClarkY and airfoils/GM15 airfoils,
same methodology as run_all_airfoils.py (SD7003/SD8000, now in airfoils/):
polar sweep (5 alphas) + grid convergence study (dx=0.04/0.02/0.01) at
dx=0.02 production resolution, same domain/dt/nsteps.

ClarkY (Re=60700) and GM15 (Re=40600) come from the same UIUC LSAT
Volume 1/3 clean-tabulated .DRG/.LFT format as SD7003/SD8000 -- see
SURF_test/airfoils/README.md for the mentor question this answers (does a
genuinely lower-Re and/or more mainstream airfoil avoid SD7003's broadband
vorticity speckle) and provenance of the extracted data.

Usage:  python3 SURF_test/run_clarky_gm15.py
Output: SURF_test/airfoils/{ClarkY,GM15}/_run_data/, SURF_test/batch_results_clarky_gm15.json
"""
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
    "ClarkY": dict(
        dat=REPO / "SURF_test" / "airfoils" / "ClarkY" / "clarky.dat.txt",
        Re=60700,
        polar_alphas=[-2.07, 1.02, 4.16, 7.19, 10.26],
        conv_alpha=-0.45,
    ),
    "GM15": dict(
        dat=REPO / "SURF_test" / "airfoils" / "GM15" / "gm15.dat.txt",
        Re=40600,
        polar_alphas=[-3.62, 0.66, 4.61, 8.75, 9.94],
        conv_alpha=0.66,
    ),
}

CONV_LEVELS = [
    dict(tag="coarse", dx=0.04, nx=150, ny=75, dt=DT, nsteps=NSTEPS),
    dict(tag="medium", dx=0.02, nx=300, ny=150, dt=DT, nsteps=NSTEPS),
    # fine level uses dt=0.005/nsteps=6000 from the start (same physical
    # time, t=30) -- SD7003/README.md documents this dx=0.01 diverges to
    # NaN at dt=0.01 in both implementations; applying the fix from the
    # start here rather than as a rerun-after-NaN follow-up.
    dict(tag="fine",   dx=0.01, nx=600, ny=300, dt=0.005, nsteps=6000),
]
POLAR_LEVEL = dict(dx=0.02, nx=300, ny=150)


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
    results_path = REPO / "SURF_test" / "batch_results_clarky_gm15.json"
    if results_path.exists():
        results = json.loads(results_path.read_text())

    for name, cfg in CASES.items():
        print(f"=== {name}: polar sweep at dx={POLAR_LEVEL['dx']} ===", flush=True)
        geom_path, npts, perim = geom_for_dx(name, cfg["dat"], POLAR_LEVEL["dx"])
        print(f"  geometry: {npts} points, perimeter={perim:.4f}", flush=True)
        results["polar"].setdefault(name, [])
        done_alphas = {r["alpha"] for r in results["polar"][name]}
        for alpha in cfg["polar_alphas"]:
            if alpha in done_alphas:
                continue
            t0 = time.time()
            outdir = REPO / "SURF_test" / "airfoils" / name / "_run_data" / f"polar_a{alpha:+.2f}"
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
        results["convergence"].setdefault(name, [])
        done_tags = {r["tag"] for r in results["convergence"][name]}
        for lvl in CONV_LEVELS:
            if lvl["tag"] in done_tags:
                continue
            geom_path, npts, perim = geom_for_dx(name, cfg["dat"], lvl["dx"])
            outdir = REPO / "SURF_test" / "airfoils" / name / "_run_data" / f"conv_{lvl['tag']}"
            fpath, elapsed = run_case(
                geom=geom_path, name="run", outdir=outdir, alpha=cfg["conv_alpha"],
                nx=lvl["nx"], ny=lvl["ny"], Re=cfg["Re"],
                dt=lvl["dt"], nsteps=lvl["nsteps"], **DOMAIN,
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
