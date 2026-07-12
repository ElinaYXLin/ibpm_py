"""
run_naca0012.py

NACA0012 ("easy": thin, symmetric, no camber, minimal laminar-separation-
bubble risk) flowfield case at Re=500 -- genuinely in the "hundreds" Re
range the mentor asked about, not just SD7003/SD8000/ClarkY/GM15's
Re~40-61k. Runs BOTH py/ibpm.py and C++ build/ibpm, same convention as
every other flowfield case in this suite (dx=0.02, dt=0.01, nsteps=3000,
restart=500 -> 7 snapshots t=0,5,...,30).

No experimental wind-tunnel Cl/Cd data exists at this Re (confirmed by
search -- see SURF_test/low_re/README.md's "Why no experimental
validation data" section); this is a Python-vs-C++ fidelity + qualitative
flow-field check only, not a polar validation against a reference dataset
the way SURF_test/airfoils/ is.

Usage: python3 SURF_test/low_re/run_naca0012.py
Output: SURF_test/low_re/NACA0012/_run_data{,_cpp}/flowfield/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
CPP_BIN = REPO / "build" / "ibpm"
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUTDIR = REPO / "SURF_test" / "low_re" / "NACA0012"

DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 500
ALPHA = 5.0
DT = 0.01
NSTEPS = 3000
RESTART = 500


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
    if which in ("cpp", "both"):
        if not CPP_BIN.exists():
            print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr)
            sys.exit(1)
        ok, elapsed = run([str(CPP_BIN)], OUTDIR / "_run_data_cpp" / "flowfield")
        print(f"NACA0012 C++: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)
    if which in ("py", "both"):
        ok, elapsed = run([sys.executable, "-u", str(RUNNER)], OUTDIR / "_run_data" / "flowfield")
        print(f"NACA0012 py: {'OK' if ok else 'FAILED'} ({elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
