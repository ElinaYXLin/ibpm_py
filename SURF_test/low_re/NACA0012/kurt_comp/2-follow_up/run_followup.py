"""
run_followup.py

"Cheap steady runs" follow-up to 1-paper_based/: targeted, py_static-only,
steady-only tests (no pitching, no cpp_static -- py/cpp agreement was
already established byte-for-byte-ish in 1-paper_based/, so these physics/
configuration sensitivity studies don't need to re-verify that) probing
WHICH ibpm configuration choices explain three of the anomalies flagged in
1-paper_based/README.md:

  D. ngrid sweep (1/2/3), alpha=0..5: does the far-field/blockage treatment
     (single- vs multi-domain) explain the ~14%-high lift-curve slope?
  E. domain-size sweep (larger uniform domain, ngrid=1), alpha=0..5: same
     question, isolated from the multi-domain method itself -- just less
     blockage from a bigger box.
  F. grid-alignment test, alpha=0 only: does shifting the domain by half a
     cell (dx/2) change/flip the small nonzero Cl(0) -- a fingerprint of a
     discrete-IB-on-Cartesian-grid symmetry limit rather than a geometry
     asymmetry (already checked: the resampled boundary itself is symmetric
     to ~1e-5).
  G. dx refinement, alpha=0 only: does the Cl(0) offset shrink at dx=0.01,
     confirming it's a resolution-limited (not fixed) effect?

All at Re=1000, steady, matching 1-paper_based's baseline (D/E/F at
nx=300,ny=150,dt=0.01,nsteps=3000; G at nx=600,ny=300,dt=0.005,nsteps=6000).
D/E/F reuse 1-paper_based's ngrid=1/length=6/yoffset=-1.5/dx=0.02 alpha=0..5
steady runs as their own baseline point (no rerun needed there).

Usage: python3 SURF_test/low_re/NACA0012/kurt_comp/2-follow_up/run_followup.py [njobs]
Output: 2-follow_up/runs/<test>/<case>/
"""
import pathlib
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
FOLLOWUP = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "2-follow_up"
RUNS = FOLLOWUP / "runs"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"
GEOM_020 = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"
GEOM_010 = REPO / "SURF_test" / "geom" / "naca0012_dx0.0100.geom"

RE = 1000
BASE_DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
ANGLES_D_E = [0, 1, 2, 3, 4, 5]


def is_done(outdir, nsteps):
    f = outdir / "flow.force"
    if not f.exists():
        return False
    try:
        last = None
        with open(f) as fh:
            for line in fh:
                if line.strip():
                    last = line
        return last is not None and int(last.split()[0]) >= nsteps
    except Exception:
        return False


def build_jobs():
    jobs = []

    # D. ngrid sweep: ngrid=2,3 (ngrid=1 already exists in 1-paper_based).
    # Grid.py asserts ngrid==1 or (nx%4==0 and ny%4==0) for the multi-domain
    # scheme -- ny=150 (1-paper_based's baseline) fails that (150%4=2), so
    # this uses ny=152 (a 1.3% taller domain: height 3.0 -> 3.04) to keep
    # dx=0.02 EXACTLY while satisfying the constraint. nx=300 already
    # satisfies it unchanged.
    ngrid_domain = dict(length=6, xoffset=-2, yoffset=-1.52)
    for ngrid in (2, 3):
        for a in ANGLES_D_E:
            outdir = RUNS / "ngrid_sweep" / f"ngrid{ngrid}_a{a:02d}"
            jobs.append(dict(
                outdir=outdir, geom=GEOM_020, nx=300, ny=152, ngrid=ngrid,
                domain=ngrid_domain, alpha=a, dt=0.01, nsteps=3000,
            ))

    # E. domain-size sweep: same dx=0.02, ~1.67x taller/wider uniform domain
    big_domain = dict(length=10, xoffset=-3.5, yoffset=-2.5)
    for a in ANGLES_D_E:
        outdir = RUNS / "domain_sweep" / f"large_a{a:02d}"
        jobs.append(dict(
            outdir=outdir, geom=GEOM_020, nx=500, ny=250, ngrid=1,
            domain=big_domain, alpha=a, dt=0.01, nsteps=3000,
        ))

    # F. grid-alignment test: shift the domain by half a cell (dx/2=0.01)
    shifted_domain = dict(length=6, xoffset=-2, yoffset=-1.5 + 0.01)
    outdir = RUNS / "grid_alignment" / "yshift_a00"
    jobs.append(dict(
        outdir=outdir, geom=GEOM_020, nx=300, ny=150, ngrid=1,
        domain=shifted_domain, alpha=0, dt=0.01, nsteps=3000,
    ))

    # G. dx refinement: alpha=0 only, at dx=0.01 (this repo's kurt_comp
    # dx=0.01 grid-check sweep was explicitly skipped for the paper
    # comparison, so this point doesn't exist yet anywhere)
    outdir = RUNS / "dx_refine" / "dx0.010_a00"
    jobs.append(dict(
        outdir=outdir, geom=GEOM_010, nx=600, ny=300, ngrid=1,
        domain=BASE_DOMAIN, alpha=0, dt=0.005, nsteps=6000,
    ))

    return jobs


def run_one(job):
    outdir = job["outdir"]
    if is_done(outdir, job["nsteps"]):
        return (str(outdir.relative_to(RUNS)), "skip", 0.0)
    outdir.mkdir(parents=True, exist_ok=True)
    d = job["domain"]
    cmd = [sys.executable, "-u", str(PY_RUNNER),
           "-geom", str(job["geom"]), "-name", "flow", "-outdir", str(outdir),
           "-nx", str(job["nx"]), "-ny", str(job["ny"]), "-ngrid", str(job["ngrid"]),
           "-length", str(d["length"]), "-xoffset", str(d["xoffset"]),
           "-yoffset", str(d["yoffset"]), "-alpha", str(job["alpha"]), "-Re", str(RE),
           "-dt", str(job["dt"]), "-nsteps", str(job["nsteps"]),
           "-tecplot", "0", "-restart", "0", "-force", "1"]
    t0 = time.time()
    with open(outdir / "run_log.txt", "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    dt_ = time.time() - t0
    status = "ok" if (proc.returncode == 0 and is_done(outdir, job["nsteps"])) else "FAIL"
    return (str(outdir.relative_to(RUNS)), status, dt_)


def main():
    njobs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    jobs = build_jobs()
    print(f"{len(jobs)} jobs, {njobs}-way parallel")
    t0 = time.time()
    done = fail = skip = 0
    with ProcessPoolExecutor(max_workers=njobs) as ex:
        futs = [ex.submit(run_one, j) for j in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            name, status, dt_ = fut.result()
            if status == "ok":
                done += 1
            elif status == "skip":
                skip += 1
            else:
                fail += 1
                print(f"  FAIL: {name}")
            if i % 5 == 0 or i == len(jobs):
                print(f"[{i}/{len(jobs)}] done={done} skip={skip} fail={fail} "
                      f"elapsed={time.time()-t0:.0f}s", flush=True)
    print(f"DONE. ran={done} skipped={skip} failed={fail} in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
