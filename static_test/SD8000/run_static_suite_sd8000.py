"""
run_static_suite_sd8000.py

Runs the SD8000 coarse-grid (dx=0.04) case 5x with py_static/ibpm.py and 5x
with build_static/ibpm (cpp_static), in sequence -- see
static_test/SD8000/README.md for why this specific case (not NACA0012's
Re=500 flow evolution) is the "hardest" setup in SURF_test/airfoils/: it's
the one documented case (SURF_test/gen_port_fidelity_diagnostic.py,
SURF_test/airfoils/7-chaos_sensitivity/README.md) where this solver is
numerically unstable and chaotically sensitive to last-bit roundoff, to the
point of a catastrophic blow-up partway through the run.

Parameters match SURF_test/run_all_airfoils.py's CONV_LEVELS[0] ("coarse")
for SD8000 exactly: dx=0.04, nx=150, ny=75, ngrid=1, domain length=6/
xoffset=-2/yoffset=-1.5, Re=60800, alpha=-0.81 (conv_alpha), dt=0.01,
nsteps=3000. restart=250 added (not in the original conv_coarse run, which
used restart=0/no snapshots) so this check also has flowfield snapshots
for the vorticity-diff figures, same convention as static_test/'s NACA0012
check.

Usage: python3 static_test/SD8000/run_static_suite_sd8000.py
Output: static_test/SD8000/py_run{1..5}/flowfield/, .../cpp_run{1..5}/flowfield/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM = REPO / "SURF_test" / "geom" / "sd8000_dx0.0400.geom"
CPP_BIN = REPO / "build_static" / "ibpm"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"
OUTBASE = REPO / "static_test" / "SD8000"

DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 60800
ALPHA = -0.81
DT = 0.01
NSTEPS = 3000
RESTART = 250
N_RUNS = 5


def run(cmd_prefix, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = cmd_prefix + [
        "-geom", str(GEOM), "-name", "flow", "-outdir", str(outdir),
        "-nx", "150", "-ny", "75", "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(ALPHA), "-Re", str(RE),
        "-dt", str(DT), "-nsteps", str(NSTEPS), "-tecplot", "0",
        "-restart", str(RESTART), "-force", "1",
    ]
    log_path = outdir / "run_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    # NOTE: this case is expected to be numerically unstable and may blow
    # up to NaN/Inf partway through -- ibpm's own returncode is still 0 in
    # that case (it doesn't detect divergence itself), so we don't treat a
    # bad returncode as expected; a nonzero returncode here is a real failure.
    return proc.returncode == 0, elapsed


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    results = []
    if which in ("py", "both"):
        for i in range(1, N_RUNS + 1):
            outdir = OUTBASE / f"py_run{i}" / "flowfield"
            ok, elapsed = run([sys.executable, "-u", str(PY_RUNNER)], outdir)
            msg = f"py run {i}: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)"
            print(msg, flush=True)
            results.append(msg)
            if not ok:
                sys.exit(1)
    if which in ("cpp", "both"):
        if not CPP_BIN.exists():
            print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr)
            sys.exit(1)
        for i in range(1, N_RUNS + 1):
            outdir = OUTBASE / f"cpp_run{i}" / "flowfield"
            ok, elapsed = run([str(CPP_BIN)], outdir)
            msg = f"cpp run {i}: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)"
            print(msg, flush=True)
            results.append(msg)
            if not ok:
                sys.exit(1)
    (OUTBASE / "run_timing_summary.txt").write_text("\n".join(results) + "\n")


if __name__ == "__main__":
    main()
