"""
run_alpha_sweep.py

Test 2 of the LE vorticity investigation (see ../README.md).

The existing polar sweep (../../run_naca0012_polar.py) runs alpha=0,2,4,6,
8,10 at Re=500, dx=0.02, but with -restart 0 (no field snapshots), so the
existing force-coefficient polar can't say anything about whether the
leading-edge vorticity speck's TOP/BOTTOM asymmetry grows with alpha (as
expected if it tracks the stagnation point migrating off the geometric
nose).

Reruns alpha=0,2,8,10 (skipping 4,6 to keep compute down; 0 and a small
and large angle bracket the trend) at Re=500, dx=0.02 -- identical
settings to the polar sweep -- but with restart snapshots enabled at the
same physical times (t=0,5,...,30) used throughout this suite's
flow_evolution.png convention.

Usage: python3 SURF_test/low_re/NACA0012/leading_edge_investigation/run_alpha_sweep.py [a1 a2 ...]
Output: SURF_test/low_re/NACA0012/leading_edge_investigation/_run_data/alpha_a<alpha>_snap/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTBASE = REPO / "SURF_test" / "low_re" / "NACA0012" / "leading_edge_investigation" / "_run_data"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 500
DT = 0.01
NSTEPS = 3000
RESTART = 500  # t=5 spacing, matches flow_evolution.png STEPS convention

ALPHAS = [0, 2, 8, 10]


def run_one(outdir, alpha):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-u", str(RUNNER),
           "-geom", str(GEOM), "-name", "flow", "-outdir", str(outdir),
           "-nx", "300", "-ny", "150", "-ngrid", "1",
           "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
           "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(alpha), "-Re", str(RE),
           "-dt", str(DT), "-nsteps", str(NSTEPS), "-tecplot", "0",
           "-restart", str(RESTART), "-force", "1"]
    log_path = outdir / "run_log.txt"
    print(f"  alpha={alpha}: -> {outdir.relative_to(REPO)}", flush=True)
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"run failed ({outdir}): see {log_path}")
    print(f"    done in {elapsed:.0f}s", flush=True)


def main():
    alphas = [float(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else ALPHAS
    for alpha in alphas:
        outdir = OUTBASE / f"alpha_a{alpha:+03.0f}_snap"
        if list(outdir.glob("flow*.bin")):
            print(f"alpha={alpha}: snapshots already present in {outdir.relative_to(REPO)}, skipping")
            continue
        run_one(outdir, alpha)
    print("done")


if __name__ == "__main__":
    main()
