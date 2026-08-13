"""
run_ngrid_sweep.py

Table 5: does far-field domain extent (ngrid) affect the LE/TE field-max
peak, with every other setting held at this folder's default (dx=0.02,
ds=dx, alpha=0, Re=1000, naca0012_baseline shape)? Nothing in
`5-leading_edge` or `../` (6-edges_further) tested ngrid on this exact
LE/TE-peak question before -- `../` Groups B-D vary phase/shape/density
at ngrid=1 only, and `3-further`'s ngrid sweeps are a different quantity
(shedding Strouhal, at angles 15/20/30/40, not the LE/TE peak at alpha=0).

ngrid=1 is NOT simply "reuse the existing baseline" here: ngrid>1
requires ny%4==0 (`py_static/grid.py`'s own assert), so ny is bumped
from 150 to 152 (yoffset -1.5 -> -1.52, a 0.04c/1.3% domain-height
change) for ALL FOUR ngrid values, including ngrid=1 -- otherwise "vary
only ngrid" wouldn't be true, since ngrid=1 would sit on a very slightly
different domain than 2/3/4. All four are launched fresh from this
identical NGRID_DOMAIN so the comparison isolates ngrid alone.

Both py_static and cpp_static, 8 runs total. dx=0.02 (nx=300, ny=152),
dt=0.01, nsteps=3000 (t=30, this folder's standard duration), alpha=0,
Re=1000, naca0012_baseline geometry (ds=dx, no densification).

Usage: python3 run_ngrid_sweep.py [njobs]
Output: runs/ngrid_sweep/ngrid{1,2,3,4}_{py,cpp}/
"""
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as c  # noqa: E402

RUNS = HERE / "runs" / "ngrid_sweep"
GEOM = c.BASE_GEOM_DX002  # naca0012_dx0.0200.geom, ds=dx baseline (same shape as Table 1-4's baseline)

DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.52)
NX, NY = 300, 152
DT, NSTEPS = 0.01, 3000
ALPHA, RE = 0.0, 1000.0


def build_jobs():
    jobs = []
    for ngrid in (1, 2, 3, 4):
        for impl in ("py", "cpp"):
            outdir = RUNS / f"ngrid{ngrid}_{impl}"
            jobs.append(dict(impl=impl, geom=GEOM, outdir=outdir, nx=NX, ny=NY,
                              dt=DT, nsteps=NSTEPS, alpha=ALPHA, re=RE,
                              domain=DOMAIN, ngrid=ngrid, restart=NSTEPS))
    return jobs


def run_one(job):
    ok, elapsed, status = c.run_case(**job)
    return (job["outdir"].name, status, elapsed)


def main():
    njobs = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    jobs = build_jobs()
    print(f"{len(jobs)} jobs, {njobs}-way parallel", flush=True)
    t0 = time.time()
    done = fail = skip = 0
    with ProcessPoolExecutor(max_workers=njobs) as ex:
        futs = [ex.submit(run_one, j) for j in jobs]
        for fut in as_completed(futs):
            name, status, elapsed = fut.result()
            if status == "skip":
                skip += 1
            elif status == "ok":
                done += 1
            else:
                fail += 1
                print(f"  FAIL: {name}", flush=True)
            print(f"  {name}: {status} ({elapsed:.0f}s) "
                  f"[{done+fail+skip}/{len(jobs)}, elapsed={time.time()-t0:.0f}s]", flush=True)
    print(f"DONE. ran={done} skipped={skip} failed={fail} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
