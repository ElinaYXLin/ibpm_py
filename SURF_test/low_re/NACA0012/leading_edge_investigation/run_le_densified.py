"""
run_le_densified.py

Test 3 of the LE vorticity investigation (see ../README.md): run the
LE-densified boundary geometry (../make_le_densified_geom.py) at the SAME
background grid dx=0.02, Re=500, alpha=0 as ../../run_gridconv.py's
dx=0.02 baseline point, with restart snapshots enabled so the resulting
vorticity field can be compared directly against that baseline's peak.

Usage: python3 SURF_test/low_re/NACA0012/leading_edge_investigation/run_le_densified.py
Output: SURF_test/low_re/NACA0012/leading_edge_investigation/_run_data/le_densified_snap/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM_DENSE = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200_LEdense.geom"
GEOM_UNIFORM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTBASE = REPO / "SURF_test" / "low_re" / "NACA0012" / "leading_edge_investigation" / "_run_data"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 500
ALPHA = 0.0
DT = 0.01
NSTEPS = 3000
RESTART = 500


def run_one(geom, outdir, label):
    if not geom.exists():
        print(f"ERROR: {geom} not found", file=sys.stderr)
        sys.exit(1)
    if list(outdir.glob("flow*.bin")):
        print(f"{label}: snapshots already present in {outdir.relative_to(REPO)}, skipping")
        return
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-u", str(RUNNER),
           "-geom", str(geom), "-name", "flow", "-outdir", str(outdir),
           "-nx", "300", "-ny", "150", "-ngrid", "1",
           "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
           "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(ALPHA), "-Re", str(RE),
           "-dt", str(DT), "-nsteps", str(NSTEPS), "-tecplot", "0",
           "-restart", str(RESTART), "-force", "1"]
    log_path = outdir / "run_log.txt"
    print(f"running {label} -> {outdir.relative_to(REPO)}", flush=True)
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"run failed: see {log_path}")
    print(f"  done in {elapsed:.0f}s", flush=True)


def main():
    # LE-densified boundary points, dx=0.02 background grid (the new case)
    run_one(GEOM_DENSE, OUTBASE / "le_densified_snap", "LE-densified")
    # uniform-boundary dx=0.02 baseline, WITH snapshots (the existing
    # ../../_run_data/gridconv_dx0.02 point used -restart 0, no snapshots)
    run_one(GEOM_UNIFORM, OUTBASE / "le_uniform_baseline_snap", "uniform baseline")


if __name__ == "__main__":
    main()
