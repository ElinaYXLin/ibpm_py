"""C++ counterpart of run_flowfield.py: runs build/ibpm on the same
SD7003/SD8000 flowfield case (same geometry, alpha, grid, dt, nsteps,
restart cadence) so gen_flowfield_figs_v2.py can show a Python-vs-C++
vorticity-field comparison.

Usage: python3 SURF_test/run_flowfield_cpp.py
Output: SURF_test/airfoils/LSAT-{SD7003,SD8000}/_run_data_cpp/flowfield/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
CPP_BIN = REPO / "build" / "ibpm"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

CASES = {
    "SD7003": dict(geom=GEOMDIR / "sd7003_dx0.0200.geom", Re=61100, alpha=4.60),
    "SD8000": dict(geom=GEOMDIR / "sd8000_dx0.0200.geom", Re=60800, alpha=5.36),
}


def main():
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found. Build it first: cd build && make", file=sys.stderr)
        sys.exit(1)

    for name, cfg in CASES.items():
        outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_cpp" / "flowfield"
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = [str(CPP_BIN),
               "-geom", str(cfg["geom"]), "-name", "flow", "-outdir", str(outdir),
               "-nx", "300", "-ny", "150", "-ngrid", "1",
               "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
               "-yoffset", str(DOMAIN["yoffset"]),
               "-alpha", str(cfg["alpha"]), "-Re", str(cfg["Re"]), "-dt", "0.01",
               "-nsteps", "3000", "-tecplot", "0", "-restart", "250", "-force", "1"]
        log_path = outdir / "flow_log.txt"
        t0 = time.time()
        with open(log_path, "w") as logf:
            proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
        elapsed = time.time() - t0
        if proc.returncode != 0:
            raise RuntimeError(f"C++ flowfield run failed ({name}): see {log_path}")
        print(name, "flowfield run (C++) done,", f"{elapsed:.1f}s")

    print("done")


if __name__ == "__main__":
    main()
