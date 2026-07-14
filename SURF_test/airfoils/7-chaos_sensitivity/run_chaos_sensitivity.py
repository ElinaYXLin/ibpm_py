"""
run_chaos_sensitivity.py

Launches, in parallel (one OS process per run -- true multi-core
parallelism, not threads), all the new simulation runs behind this
directory's investigation:

1. SD7003 + SD8000 coarse-grid (dx=0.04) `conv_coarse` case (see
   ../2-c++included/port_fidelity_diagnostic.py for the original
   3000-step version), extended from 3000 -> 4000 steps, for BOTH
   py/ibpm.py and C++ build/ibpm, with a fine restart interval (every 25
   steps) so vorticity snapshots near any blow-up can be extracted
   afterward. (4 runs)
2. A tiny-Re-perturbation ensemble, C++ only, at SD8000's coarse config
   (the case already confirmed to blow up predictably near step 3000),
   varying Re by relative amounts of only 5e-8 to 1e-5 -- to test whether
   the blow-up step is chaotically sensitive to that (vs. fixed/
   deterministic). (16 runs)

The original coarse-grid runs used `-restart 0` (no checkpoints), so
"extending" them requires a fresh 0 -> 4000 rerun of each -- same
Re/alpha/grid as the original, nothing else changed.

Usage: python3 SURF_test/airfoils/7-chaos_sensitivity/run_chaos_sensitivity.py
(run from the repository root; requires build/ibpm. Takes ~70-90s on a
10-core machine; runs are launched concurrently via subprocess.Popen.)
"""
from __future__ import annotations

import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
CPP_BIN = REPO / "build" / "ibpm"
PY_RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
OUT = REPO / "SURF_test" / "airfoils" / "7-chaos_sensitivity" / "_run_data"
OUT.mkdir(parents=True, exist_ok=True)

COMMON = ["-tecplot", "0", "-force", "1", "-ngrid", "1", "-length", "6.0",
          "-xoffset", "-2.0", "-yoffset", "-1.5", "-nx", "150", "-ny", "75", "-dt", "0.01"]

# Same alpha/Re/dx as ../2-c++included/port_fidelity_diagnostic.py's
# conv_coarse case for each airfoil.
AIRFOILS = {
    "SD7003": dict(geom=REPO / "SURF_test" / "geom" / "sd7003_dx0.0400.geom", alpha=-0.09, Re=61100.0),
    "SD8000": dict(geom=REPO / "SURF_test" / "geom" / "sd8000_dx0.0400.geom", alpha=-0.81, Re=60800.0),
}

# Relative Re perturbations for the ensemble (SD8000 config only) -- the
# 6th-8th significant digit of Re=60800.
REL_PERTURBATIONS = [-1e-5, -5e-6, -2e-6, -1e-6, -5e-7, -2e-7, -1e-7, -5e-8,
                      5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5]


def build_cmd(impl, geom, alpha, Re, nsteps, restart, outdir):
    prefix = [str(CPP_BIN)] if impl == "cpp" else [sys.executable, "-u", str(PY_RUNNER)]
    return prefix + [
        "-geom", str(geom), "-name", "run", "-outdir", str(outdir),
        "-alpha", str(alpha), "-Re", str(Re), "-nsteps", str(nsteps), "-restart", str(restart),
        *COMMON,
    ]


def launch(cmd, outdir, log_name="run_log.txt"):
    outdir.mkdir(parents=True, exist_ok=True)
    log = open(outdir / log_name, "w")
    proc = subprocess.Popen(cmd, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
    return proc, log


def main():
    jobs = []  # (label, proc, logfile)

    # ---- 1. Extend SD7003 + SD8000 conv_coarse to 4000 steps, restart every 25 ----
    for name, cfg in AIRFOILS.items():
        for impl in ("py", "cpp"):
            outdir = OUT / f"{name}_ext4000_{impl}"
            cmd = build_cmd(impl, cfg["geom"], cfg["alpha"], cfg["Re"], 4000, 25, outdir)
            proc, log = launch(cmd, outdir)
            jobs.append((f"{name}-ext4000-{impl}", proc, log))
            print(f"launched {name}-ext4000-{impl}", flush=True)

    # ---- 2. Perturbation ensemble: C++ only, SD8000 coarse config, tiny Re shifts ----
    base_re = AIRFOILS["SD8000"]["Re"]
    for i, rel in enumerate(REL_PERTURBATIONS):
        re_val = base_re * (1.0 + rel)
        outdir = OUT / f"SD8000_perturb_{i:02d}"
        cmd = build_cmd("cpp", AIRFOILS["SD8000"]["geom"], AIRFOILS["SD8000"]["alpha"],
                         re_val, 4000, 0, outdir)
        proc, log = launch(cmd, outdir)
        jobs.append((f"perturb-{i:02d}-Re{re_val:.6f}", proc, log))
        print(f"launched perturb-{i:02d} Re={re_val:.6f} (rel={rel:+.1e})", flush=True)

    print(f"\n{len(jobs)} jobs launched, waiting for completion...", flush=True)
    t0 = time.time()
    remaining = list(jobs)
    while remaining:
        time.sleep(2)
        still = []
        for label, proc, log in remaining:
            if proc.poll() is None:
                still.append((label, proc, log))
            else:
                log.close()
                status = "OK" if proc.returncode == 0 else f"FAILED(rc={proc.returncode})"
                print(f"  [{time.time()-t0:6.1f}s] {label}: {status}", flush=True)
        remaining = still
    print(f"\nALL DONE in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
