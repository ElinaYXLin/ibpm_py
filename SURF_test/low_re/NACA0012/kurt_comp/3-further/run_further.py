"""
run_further.py

New simulations for the 7 further tests proposed as follow-up to
2-follow_up/README.md's "Open questions" (1a/1b/2a/3a need no new runs --
see analyze_further.py). Steady, py_static only, Re=1000 baseline (see
1-paper_based/README.md for why py-only is enough here).

  2b. dx refinement of shedding frequency at 4 representative angles
      (15, 20, 30, 40 -- inside/near each Strouhal plateau/transition),
      at dx=0.01; plus ONE angle (20, the low-St plateau) also at
      dx=0.005 -- this single point is the expensive "long" job (~3h).
  2c. ngrid=2,3 sweep of shedding frequency at the same 4 angles, at the
      baseline dx=0.02 grid (ny=152 for the ngrid>1 divisible-by-4
      requirement, same fix as 2-follow_up's ngrid_sweep).
  3b. Initial-condition ensemble at 3 adjacent "jagged" angles (25, 28,
      30 deg): 5 "source" runs (steady runs at the neighboring angles
      24/26/27/29/31, WITH a final restart snapshot -- 1-paper_based only
      kept snapshots for 0/9/12 deg, so these don't exist yet), 6
      "approach-direction" runs (continue from a neighbor's developed
      state into the target angle, from below and from above), and 9
      "perturbation" runs (impulsive start, Re nudged by a small relative
      amount, three angles x three perturbation levels).

  2d. Does refining dx AND increasing ngrid together do better than either
      alone? Combines 2b's dx=0.01 with 2c's ngrid=2,3 (same 4 angles,
      same ngrid-domain-height fix as 2c, scaled to dx=0.01: ny=304
      instead of 152, keeping the same 3.04c physical domain height and
      the ngrid>1 divisible-by-4 requirement). dx=0.005 excluded per the
      mentor's/user's own call -- too slow to be worth it here.

Usage:
  python3 run_further.py short [njobs]   -- everything except the dx=0.005 point
  python3 run_further.py long            -- just the dx=0.005 point (serial, 1 job)
  python3 run_further.py test2d [njobs]  -- 2d: dx=0.01 x ngrid=2,3 (8 jobs)
Output: 3-further/runs/<test>/<case>/
"""
import pathlib
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
FURTHER = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "3-further"
RUNS = FURTHER / "runs"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"
PAPER_RUNS = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "1-paper_based" / "runs" / "dx0.020"

GEOM = {
    "0.02": REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom",
    "0.01": REPO / "SURF_test" / "geom" / "naca0012_dx0.0100.geom",
    "0.005": REPO / "SURF_test" / "geom" / "naca0012_dx0.0050.geom",
}
RE = 1000
BASE_DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
REP_ANGLES = [15, 20, 30, 40]  # 2b/2c representative angles
JAGGED_ANGLES = [25, 28, 30]   # 3b target angles
NEIGHBORS = {25: (24, 26), 28: (27, 29), 30: (29, 31)}
PERTURBATIONS = [-1e-2, -1e-3, 1e-2]  # relative Re nudges


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


def run_case(outdir, geom, nx, ny, ngrid, domain, alpha, Re, dt, nsteps,
             restart=0, ic=None, resettime=False):
    if is_done(outdir, nsteps):
        return (str(outdir.relative_to(RUNS)), "skip", 0.0)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-u", str(PY_RUNNER),
           "-geom", str(geom), "-name", "flow", "-outdir", str(outdir),
           "-nx", str(nx), "-ny", str(ny), "-ngrid", str(ngrid),
           "-length", str(domain["length"]), "-xoffset", str(domain["xoffset"]),
           "-yoffset", str(domain["yoffset"]), "-alpha", str(alpha), "-Re", str(Re),
           "-dt", str(dt), "-nsteps", str(nsteps), "-tecplot", "0",
           "-restart", str(restart), "-force", "1"]
    if ic is not None:
        cmd += ["-ic", str(ic)]
    if resettime:
        cmd += ["-resettime", "1"]
    t0 = time.time()
    with open(outdir / "run_log.txt", "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    dt_ = time.time() - t0
    status = "ok" if (proc.returncode == 0 and is_done(outdir, nsteps)) else "FAIL"
    return (str(outdir.relative_to(RUNS)), status, dt_)


def build_short_jobs():
    jobs = []

    # 2b: dx=0.01 at 4 representative angles
    for a in REP_ANGLES:
        outdir = RUNS / "dx_refine" / f"dx0.010_a{a:02d}"
        jobs.append(dict(outdir=outdir, geom=GEOM["0.01"], nx=600, ny=300, ngrid=1,
                          domain=BASE_DOMAIN, alpha=a, Re=RE, dt=0.005, nsteps=6000))

    # 2c: ngrid=2,3 at same 4 angles (ny=152 domain, matching 2-follow_up's fix)
    ngrid_domain = dict(length=6, xoffset=-2, yoffset=-1.52)
    for ngrid in (2, 3):
        for a in REP_ANGLES:
            outdir = RUNS / "ngrid_sweep" / f"ngrid{ngrid}_a{a:02d}"
            jobs.append(dict(outdir=outdir, geom=GEOM["0.02"], nx=300, ny=152, ngrid=ngrid,
                              domain=ngrid_domain, alpha=a, Re=RE, dt=0.01, nsteps=3000))

    # 3b source runs: neighbor angles, restart only at the final step
    neighbor_angles = sorted(set(n for pair in NEIGHBORS.values() for n in pair))
    for a in neighbor_angles:
        outdir = RUNS / "ic_ensemble" / "source" / f"steady_a{a:02d}"
        jobs.append(dict(outdir=outdir, geom=GEOM["0.02"], nx=300, ny=150, ngrid=1,
                          domain=BASE_DOMAIN, alpha=a, Re=RE, dt=0.01, nsteps=3000,
                          restart=3000))

    # 3b approach-direction runs: queued after source runs by dependency,
    # but run_one resolves the ic path lazily at call time, so it's safe to
    # submit everything to the pool at once as long as sources finish first
    # within the same pool -- to guarantee that ordering simply, run these
    # in a SEPARATE, later pool stage (see main()).

    return jobs


def build_approach_jobs():
    jobs = []
    for target, (below, above) in NEIGHBORS.items():
        for direction, neighbor in (("from_below", below), ("from_above", above)):
            src_bin = RUNS / "ic_ensemble" / "source" / f"steady_a{neighbor:02d}" / "flow03000.bin"
            outdir = RUNS / "ic_ensemble" / "approach" / f"a{target:02d}_{direction}"
            jobs.append(dict(outdir=outdir, geom=GEOM["0.02"], nx=300, ny=150, ngrid=1,
                              domain=BASE_DOMAIN, alpha=target, Re=RE, dt=0.01, nsteps=3000,
                              ic=src_bin, resettime=True))
    return jobs


def build_perturbation_jobs():
    jobs = []
    for a in JAGGED_ANGLES:
        for rel in PERTURBATIONS:
            re_val = RE * (1 + rel)
            tag = f"rel{rel:+.0e}".replace("+", "p").replace("-", "m")
            outdir = RUNS / "ic_ensemble" / "perturb" / f"a{a:02d}_{tag}"
            jobs.append(dict(outdir=outdir, geom=GEOM["0.02"], nx=300, ny=150, ngrid=1,
                              domain=BASE_DOMAIN, alpha=a, Re=re_val, dt=0.01, nsteps=3000))
    return jobs


def build_test2d_jobs():
    """2d: dx=0.01 combined with ngrid=2,3 at the same 4 representative
    angles. Domain height widened from BASE_DOMAIN's 3.0c to 3.04c (same
    fix 2c used, here at 2x resolution: ny=304 instead of 152) purely to
    satisfy the ngrid>1 divisible-by-4 requirement -- not a physical
    change vs. the dx=0.01/ngrid=1 baseline already on disk from 2b,
    which used yoffset=-1.5 (3.00c) instead; the 0.04c difference in
    domain height is negligible next to the resolution/far-field
    variables actually being tested."""
    jobs = []
    ngrid_domain = dict(length=6, xoffset=-2, yoffset=-1.52)
    for ngrid in (2, 3):
        for a in REP_ANGLES:
            outdir = RUNS / "dx_ngrid_sweep" / f"dx0.010_ngrid{ngrid}_a{a:02d}"
            jobs.append(dict(outdir=outdir, geom=GEOM["0.01"], nx=600, ny=304, ngrid=ngrid,
                              domain=ngrid_domain, alpha=a, Re=RE, dt=0.005, nsteps=6000))
    return jobs


def build_long_job():
    outdir = RUNS / "dx_refine" / "dx0.005_a20"
    return dict(outdir=outdir, geom=GEOM["0.005"], nx=1200, ny=600, ngrid=1,
                domain=BASE_DOMAIN, alpha=20, Re=RE, dt=0.0025, nsteps=12000)


def run_pool(jobs, njobs, label):
    print(f"{label}: {len(jobs)} jobs, {njobs}-way parallel")
    t0 = time.time()
    done = fail = skip = 0
    with ProcessPoolExecutor(max_workers=njobs) as ex:
        futs = [ex.submit(run_case, **j) for j in jobs]
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
    print(f"{label} DONE. ran={done} skipped={skip} failed={fail} in {time.time()-t0:.0f}s")
    return fail == 0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "short"
    njobs = int(sys.argv[2]) if len(sys.argv) > 2 else 8

    if mode == "long":
        ok = run_pool([build_long_job()], 1, "2b (dx=0.005)")
        sys.exit(0 if ok else 1)

    if mode == "test2d":
        ok = run_pool(build_test2d_jobs(), njobs, "2d (dx=0.01 x ngrid=2,3)")
        sys.exit(0 if ok else 1)

    # short mode: 2b(dx=0.01) + 2c(ngrid) + 3b(source) in one pool, THEN
    # 3b(approach) in a second pool (needs the source runs' .bin files to
    # exist first), THEN 3b(perturb) (independent, could run anytime, but
    # kept last to keep output easy to read)
    stage1 = build_short_jobs()
    ok1 = run_pool(stage1, njobs, "stage 1 (2b dx=0.01, 2c ngrid, 3b sources)")

    stage2 = build_approach_jobs()
    ok2 = run_pool(stage2, njobs, "stage 2 (3b approach-direction)")

    stage3 = build_perturbation_jobs()
    ok3 = run_pool(stage3, njobs, "stage 3 (3b perturbations)")

    sys.exit(0 if (ok1 and ok2 and ok3) else 1)


if __name__ == "__main__":
    main()
