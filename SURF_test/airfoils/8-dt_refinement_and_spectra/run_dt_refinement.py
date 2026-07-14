"""
run_dt_refinement.py

Follow-up to ../7-chaos_sensitivity/'s finding that SD8000's coarse grid
(dx=0.04) blows up (goes to NaN) within a few thousand steps, and that
the exact blow-up step is scrambled by imperceptible (1e-8 to 1e-5
relative) Re perturbations -- taken there as evidence of "chaos." That
conclusion has a gap: a NaN blow-up is the signature of a numerical
(discretization/CFL) instability, not of genuine chaotic physics (which
stays bounded forever, just unpredictable). The two are easy to conflate
because BOTH make the blow-up time sensitive to tiny perturbations.

This script runs the single most direct test to tell them apart: hold
dx=0.04 fixed and refine dt. If the blow-up is a temporal (CFL-type)
numerical instability, refining dt should delay or remove it; if it's
genuinely physical, the blow-up step (in physical TIME, not iteration
count) should be roughly dt-independent.

dt = 0.01 (the original), 0.005, 0.0025, each run 5 times (C++ only,
identical command each time) to separate genuine dt-dependence from the
run-to-run FFTW_EXHAUSTIVE replanning noise already documented in
../7-chaos_sensitivity/README.md (FFTW re-times its candidate plans at
each fresh process start and can select a different, numerically
non-identical-but-equivalent plan -- this alone perturbs the trajectory
at the round-off level, enough to reshuffle blow-up timing near a
marginal instability, so a single run per dt can't be trusted).

All runs held to the same physical end time t_final=40 (matching
../7-chaos_sensitivity's 4000-step, dt=0.01 extension), nsteps scaled per
dt to reach it. restart=0 (force traces only -- the companion spectrum/
energy analysis in compute_spectrum_and_energy.py reuses the field
snapshots ../7-chaos_sensitivity already saved, no new snapshots needed
here).

Usage: python3 SURF_test/airfoils/8-dt_refinement_and_spectra/run_dt_refinement.py
(run from repo root; requires build/ibpm; launches all 15 runs concurrently
via subprocess.Popen, same pattern as ../7-chaos_sensitivity/run_chaos_sensitivity.py)
"""
from __future__ import annotations

import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
CPP_BIN = REPO / "build" / "ibpm"
OUT = REPO / "SURF_test" / "airfoils" / "8-dt_refinement_and_spectra" / "_run_data"
OUT.mkdir(parents=True, exist_ok=True)

GEOM = REPO / "SURF_test" / "geom" / "sd8000_dx0.0400.geom"
ALPHA = -0.81
RE = 60800.0
T_FINAL = 40.0
N_REPEATS = 5
DTS = [0.01, 0.005, 0.0025]

COMMON = ["-tecplot", "0", "-force", "1", "-ngrid", "1", "-length", "6.0",
          "-xoffset", "-2.0", "-yoffset", "-1.5", "-nx", "150", "-ny", "75",
          "-restart", "0"]


def build_cmd(dt, nsteps, outdir):
    return [str(CPP_BIN), "-geom", str(GEOM), "-name", "run", "-outdir", str(outdir),
            "-alpha", str(ALPHA), "-Re", str(RE), "-dt", str(dt), "-nsteps", str(nsteps),
            *COMMON]


def launch(cmd, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    log = open(outdir / "run_log.txt", "w")
    proc = subprocess.Popen(cmd, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
    return proc, log


def main():
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found.", file=sys.stderr)
        sys.exit(1)

    jobs = []
    for dt in DTS:
        nsteps = int(round(T_FINAL / dt))
        for rep in range(N_REPEATS):
            outdir = OUT / f"dt{dt}_rep{rep:02d}"
            cmd = build_cmd(dt, nsteps, outdir)
            proc, log = launch(cmd, outdir)
            jobs.append((f"dt={dt} rep{rep:02d} (nsteps={nsteps})", proc, log))
            print(f"launched dt={dt} rep{rep:02d}: nsteps={nsteps} -> {outdir.relative_to(REPO)}", flush=True)

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
