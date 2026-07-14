"""
run_sd_gridconv_extended.py

Extends run_all_airfoils.py / run_all_airfoils_cpp.py's SD7003/SD8000
grid-convergence sweep (dx=0.04/0.02/0.01, at each airfoil's conv_alpha)
with finer levels (dx=0.005, then dx=0.0025 if warranted), following the
same convention already used for NACA0012 (../low_re/NACA0012/run_gridconv.py):
dt = dx/2 for dx < 0.04, nsteps scaled to hold t_final=30 fixed.

Runs ONE (impl, airfoil) combination per invocation, so the caller can
launch py/SD7003, py/SD8000, cpp/SD7003, cpp/SD8000 as four separate
background processes on separate CPU cores simultaneously instead of
sequentially.

Usage: python3 SURF_test/run_sd_gridconv_extended.py <py|cpp> <SD7003|SD8000> [dx1 dx2 ...]
       (defaults to dx=0.005 if none given)
Output: SURF_test/airfoils/LSAT-<name>/_run_data{,_cpp}/conv_dx<dx>/
        SURF_test/batch_results_extended_<impl>.json
"""
import json
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO / "SURF_test"))
from make_airfoil_raw import make_raw_for_dx  # noqa: E402

GEOMDIR = REPO / "SURF_test" / "geom"
CPP_BIN = REPO / "build" / "ibpm"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
T_FINAL = 30.0
AVG_FRAC = 0.6

CASES = {
    "SD7003": dict(
        dat=REPO / "SURF_test" / "airfoils" / "LSAT-SD7003" / "sd7003.dat.txt",
        Re=61100, conv_alpha=-0.09,
    ),
    "SD8000": dict(
        dat=REPO / "SURF_test" / "airfoils" / "LSAT-SD8000" / "sd8000.dat.txt",
        Re=60800, conv_alpha=-0.81,
    ),
}


def geom_for_dx(name, dat, dx):
    raw_path = GEOMDIR / f"{name.lower()}_dx{dx:.4f}.txt"
    geom_path = GEOMDIR / f"{name.lower()}_dx{dx:.4f}.geom"
    if not geom_path.exists():
        n, perim = make_raw_for_dx(str(dat), dx, str(raw_path))
        geom_path.write_text(f"body {name}\n  raw {raw_path}\n  center 0.25 0.0\nend\n")
        print(f"  generated {geom_path.name}: {n} points, perimeter={perim:.4f}", flush=True)
    return geom_path


def dt_nsteps_for(dx):
    dt = dx / 2  # matches ../low_re/NACA0012/run_gridconv.py's rule for dx<0.04
    nsteps = int(round(T_FINAL / dt))
    return dt, nsteps


def time_avg_force(force_path, frac=AVG_FRAC):
    import numpy as np
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    n = len(d)
    seg = d[int(n * (1 - frac)):]
    return dict(cd_mean=float(seg[:, 2].mean()), cl_mean=float(seg[:, 3].mean()),
                cd_std=float(seg[:, 2].std()), cl_std=float(seg[:, 3].std()),
                t_final=float(d[-1, 1]), n_steps=int(n))


def run_one(impl, name, dx):
    cfg = CASES[name]
    geom_path = geom_for_dx(name, cfg["dat"], dx)
    dt, nsteps = dt_nsteps_for(dx)
    nx = int(round(DOMAIN["length"] / dx))
    ny = int(round(3.0 / dx))
    subdir = "_run_data_cpp" if impl == "cpp" else "_run_data"
    outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / subdir / f"conv_dx{dx}"
    fpath = outdir / ("run.force" if impl == "cpp" else "run.force")
    if fpath.exists():
        stats = time_avg_force(fpath)
        print(f"[{impl}/{name}] dx={dx}: already present, Cd={stats['cd_mean']:+.5f} "
              f"Cl={stats['cl_mean']:+.5f} (skipping)", flush=True)
        return stats

    outdir.mkdir(parents=True, exist_ok=True)
    cmd_prefix = [str(CPP_BIN)] if impl == "cpp" else [sys.executable, "-u", str(RUNNER)]
    cmd = cmd_prefix + [
        "-geom", str(geom_path), "-name", "run", "-outdir", str(outdir),
        "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(cfg["conv_alpha"]), "-Re", str(cfg["Re"]),
        "-dt", str(dt), "-nsteps", str(nsteps), "-tecplot", "0", "-restart", "0", "-force", "1",
    ]
    log_path = outdir / "run_log.txt"
    print(f"[{impl}/{name}] dx={dx}: nx={nx} ny={ny} dt={dt} nsteps={nsteps} "
          f"alpha={cfg['conv_alpha']} Re={cfg['Re']} -> {outdir.relative_to(REPO)}", flush=True)
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"run failed ({outdir}): see {log_path}")
    stats = time_avg_force(fpath)
    stats.update(dx=dx, elapsed=elapsed, nx=nx, ny=ny, alpha=cfg["conv_alpha"])
    print(f"[{impl}/{name}] dx={dx}: Cd={stats['cd_mean']:+.5f} Cl={stats['cl_mean']:+.5f} "
          f"({elapsed:.0f}s)", flush=True)
    return stats


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    impl, name = sys.argv[1], sys.argv[2]
    assert impl in ("py", "cpp"), impl
    assert name in CASES, name
    if impl == "cpp" and not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr)
        sys.exit(1)
    dxs = [float(a) for a in sys.argv[3:]] if len(sys.argv) > 3 else [0.005]

    # NOTE: one file PER (impl, name), not shared across airfoils -- an
    # earlier version shared batch_results_extended_<impl>.json across both
    # SD7003 and SD8000 processes running concurrently, and the last writer
    # silently clobbered the other's entry (confirmed: SD8000's dx=0.005
    # result vanished from batch_results_extended_py.json after both
    # finished). Recovered that time by recomputing straight from run.force;
    # avoided going forward by giving each (impl, name) its own file.
    results_path = REPO / "SURF_test" / f"batch_results_extended_{impl}_{name}.json"
    results = json.loads(results_path.read_text()) if results_path.exists() else {}

    for dx in dxs:
        stats = run_one(impl, name, dx)
        results[str(dx)] = stats
        results_path.write_text(json.dumps(results, indent=2))

    print(f"[{impl}/{name}] done, wrote {results_path}", flush=True)


if __name__ == "__main__":
    main()
