"""
run_static_suite.py

Runs the standard NACA0012 flow-evolution case (SURF_test/low_re/run_naca0012.py's
exact parameters: dx=0.02, nx=300, ny=150, ngrid=1, domain length=6/xoffset=-2/
yoffset=-1.5, Re=500, alpha=5 deg, dt=0.01, nsteps=3000, restart=500 -> 7
snapshots t=0,5,...,30) 5 times with py_static/ibpm.py and 5 times with
build_static/ibpm (cpp_static), in sequence, to check that fixing the DST
planner flag (FFTW_EXHAUSTIVE -> FFTW_ESTIMATE|FFTW_UNALIGNED, see
py_static/_fftw_dst_shim.c and cpp_static/EllipticSolver2d.cc) makes each
implementation's own output bit-identical run-to-run.

Usage: python3 static_test/run_static_suite.py
Output: static_test/py_run{1..5}/flowfield/, static_test/cpp_run{1..5}/flowfield/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
CPP_BIN = REPO / "build_static" / "ibpm"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"
OUTBASE = REPO / "static_test"

DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 500
ALPHA = 5.0
DT = 0.01
NSTEPS = 3000
RESTART = 500
N_RUNS = 5


def run(cmd_prefix, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = cmd_prefix + [
        "-geom", str(GEOM), "-name", "flow", "-outdir", str(outdir),
        "-nx", "300", "-ny", "150", "-ngrid", "1",
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
