"""
run_naca0012_polar.py

Validates py/ibpm.py and C++ build/ibpm on NACA0012 at genuinely low
(hundreds) Reynolds number, against a NON-LSAT reference: published
low-Reynolds-number CFD benchmark drag coefficients (no wind-tunnel data
exists this low -- see SURF_test/airfoils/Lockard-NACA0012/README.md).

Reference drag anchors (alpha=0, where lift is zero by NACA0012 symmetry):
  Re=500:  Cd = 0.1762 (Lockard et al.), 0.1759 (Wu et al.), 0.178 (Nita et al. LBM)
  Re=1000: Cd = 0.119  (Di Ilio et al. HLBM & XFOIL; Kurtulus ~0.12)

Runs a small angle-of-attack polar at Re=500 (0,2,4,6,8,10 deg) so both
lift AND drag curves are produced, plus a Re=1000 alpha=0 point for a
second independent drag anchor. Same convention as the other airfoil
polars in this suite: dx=0.02 (nx=300, ny=150), dt=0.01, nsteps=3000
(t=30), Cl/Cd time-averaged over the last 60% of the run.

Usage: python3 SURF_test/run_naca0012_polar.py
Output: SURF_test/airfoils/Lockard-NACA0012/_run_data{,_cpp}/
        SURF_test/airfoils/Lockard-NACA0012/naca0012_polar_results.json
"""
import json
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
CPP_BIN = REPO / "build" / "ibpm"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTBASE = REPO / "SURF_test" / "airfoils" / "Lockard-NACA0012"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

DT = 0.01
NSTEPS = 3000
AVG_FRAC = 0.6
# (Re, alpha) cases: Re=500 polar for lift+drag curves, plus Re=1000 a=0 anchor
CASES = [(500, a) for a in (0, 2, 4, 6, 8, 10)] + [(1000, 0)]


def time_avg_force(force_path, frac=AVG_FRAC):
    import numpy as np
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    n = len(d)
    seg = d[int(n * (1 - frac)):]
    return dict(cl_mean=float(seg[:, 3].mean()), cl_std=float(seg[:, 3].std()),
                cd_mean=float(seg[:, 2].mean()), cd_std=float(seg[:, 2].std()),
                n_steps=int(n))


def run_one(cmd_prefix, outdir, Re, alpha):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = cmd_prefix + [
        "-geom", str(GEOM), "-name", "run", "-outdir", str(outdir),
        "-nx", "300", "-ny", "150", "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(alpha), "-Re", str(Re),
        "-dt", str(DT), "-nsteps", str(NSTEPS), "-tecplot", "0", "-restart", "0", "-force", "1",
    ]
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"run failed ({outdir}): see {log_path}")
    return outdir / "run.force", time.time() - t0


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    results_path = OUTBASE / "naca0012_polar_results.json"
    results = json.loads(results_path.read_text()) if results_path.exists() else {"py": [], "cpp": []}

    def do(impl, cmd_prefix, subdir):
        done = {(r["Re"], r["alpha"]) for r in results[impl]}
        for Re, alpha in CASES:
            if (Re, alpha) in done:
                continue
            outdir = OUTBASE / subdir / f"Re{Re}_a{alpha:+03.0f}"
            fpath, elapsed = run_one(cmd_prefix, outdir, Re, alpha)
            stats = time_avg_force(fpath)
            stats.update(Re=Re, alpha=alpha, elapsed=elapsed)
            results[impl].append(stats)
            print(f"[{impl}] Re={Re} a={alpha:+.0f}  Cl={stats['cl_mean']:+.4f}  "
                  f"Cd={stats['cd_mean']:+.4f}  ({elapsed:.0f}s)", flush=True)
            results_path.parent.mkdir(parents=True, exist_ok=True)
            results_path.write_text(json.dumps(results, indent=2))

    if which in ("cpp", "both"):
        if not CPP_BIN.exists():
            print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr); sys.exit(1)
        do("cpp", [str(CPP_BIN)], "_run_data_cpp")
    if which in ("py", "both"):
        do("py", [sys.executable, "-u", str(RUNNER)], "_run_data")
    print("done")


if __name__ == "__main__":
    main()
