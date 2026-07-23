"""
run_all_6.py

Single orchestrator for every new simulation this folder needs: both
py_static and cpp_static, launched via a bounded process pool (this
machine has 10 cores). Resumable -- a job whose target step count is
already on disk is skipped.

Job groups:
  recon2  -- LE-only dense @ Re=1000 (2 impls); LE+TE dense @ Re=500 (2 impls)
  B1      -- NACA0012 baseline, grid phase-shifted (7 shifts x 2 impls)
  B2      -- NACA0006/0012/0018, phase-equalized to a common LE sub-cell
             phase (built from testA2's audit -- see phase_shift_for())
  C1      -- NACA0006/0018 at dx=0.01, dx=0.005 (2 shapes x 2 dx x 2 impls)
  C2      -- thickness family at dx=0.02 (6 shapes x 2 impls)
  D1      -- NACA0012 LE point-density levels (4 factors x 2 impls;
             factor=4 reuses recon2's LEonly_dense run)
  E1      -- front-only / TE-only decoupled variants (3 shapes x 2 impls)
  E2      -- common-TE resweep (2 shapes x 2 impls)

Usage: python3 run_all_6.py [njobs]
Output: runs/<group>/<case>[_cpp]/
"""
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

import common as c

ALPHA0 = 0.0
RE1000 = 1000.0
RE500 = 500.0
BASE_DOMAIN = dict(length=6.0, xoffset=-2.0, yoffset=-1.5)
IMPLS = ("py", "cpp")


def outdir(group, name, impl):
    return c.RUNS / group / (name if impl == "py" else f"{name}_cpp")


def job(group, name, geom, nx, ny, dt, nsteps, impl, re=RE1000, alpha=ALPHA0, domain=None):
    od = outdir(group, name, impl)
    # restart=nsteps -> exactly one snapshot, at the final (developed) step;
    # restart=0 (the default) writes NO restart file at all, ever, which
    # silently makes every run look "failed" to is_done()'s .bin check.
    ok, elapsed, status = c.run_case(impl, geom, od, nx, ny, dt, nsteps, alpha=alpha, re=re,
                                       domain=domain or BASE_DOMAIN, restart=nsteps)
    return f"{group}/{name} {impl}", status, elapsed


def build_jobs():
    jobs = []

    # ---------------- recon2 ----------------
    le_only = c.GEOMDIR / "naca0012_dx0.0200_LEonly_dense.geom"
    lte_dense = c.KURT5 / "geom" / "naca0012_dx0.0200_LTEdense.geom"
    for impl in IMPLS:
        jobs.append(("recon2", "LEonly_dense_Re1000", le_only, 300, 150, 0.01, 3000, impl, RE1000, ALPHA0, None))
        jobs.append(("recon2", "LTEdense_Re500", lte_dense, 300, 150, 0.01, 3000, impl, RE500, ALPHA0, None))

    # ---------------- B1: phase sweep, NACA0012 baseline ----------------
    dx = 0.02
    shifts = [(0.25, 0), (0.5, 0), (0.75, 0), (0, 0.25), (0, 0.5), (0, 0.75), (0.5, 0.5)]
    for fx, fy in shifts:
        name = f"phase_x{fx:.2f}_y{fy:.2f}"
        domain = dict(length=6.0, xoffset=-2.0 + fx * dx, yoffset=-1.5 + fy * dx)
        for impl in IMPLS:
            jobs.append(("B1_phase", name, c.BASE_GEOM_DX002, 300, 150, 0.01, 3000, impl, RE1000, ALPHA0, domain))

    # ---------------- B2: phase-equalized 0006/0012/0018 ----------------
    # computed from testA2's audit (phase_x of each shape's LE, relative to
    # the standard xoffset=-2.0 grid) -- shift each shape's OWN grid so its
    # LE lands at the SAME target phase (use naca0012's own native phase as
    # the common target, since it needs no shift)
    import csv as _csv
    audit = {}
    audit_path = c.DATA / "testA2_geometry_audit.csv"
    if audit_path.exists():
        with open(audit_path) as f:
            for row in _csv.DictReader(f):
                audit[row["shape"]] = row
    if audit:
        target_phase_x = float(audit["naca0012"]["phase_x"])
        b2_shapes = {
            "naca0006": c.KURT5 / "geom" / "naca0006_dx0.0200.geom",
            "naca0012": c.BASE_GEOM_DX002,
            "naca0018": c.KURT5 / "geom" / "naca0018_dx0.0200.geom",
        }
        for shape, geom in b2_shapes.items():
            phase_x = float(audit[shape]["phase_x"])
            # shift xoffset so this shape's LE phase matches the target
            dphase = target_phase_x - phase_x
            domain = dict(length=6.0, xoffset=-2.0 + dphase * dx, yoffset=-1.5)
            for impl in IMPLS:
                jobs.append(("B2_phase_equalized", f"{shape}_equalized", geom, 300, 150, 0.01, 3000,
                             impl, RE1000, ALPHA0, domain))

    # ---------------- C1: per-shape refinement ----------------
    for shape in ("naca0006", "naca0018"):
        for dxv, nx, ny, dt, nsteps in [(0.01, 600, 300, 0.005, 6000), (0.005, 1200, 600, 0.0025, 12000)]:
            geom = c.GEOMDIR / f"{shape}_dx{dxv:.4f}.geom"
            for impl in IMPLS:
                jobs.append(("C1_shape_refine", f"{shape}_dx{dxv:.4f}", geom, nx, ny, dt, nsteps,
                             impl, RE1000, ALPHA0, None))

    # ---------------- C2: thickness family @ dx=0.02 ----------------
    for name in ("naca0004", "naca0008", "naca0010", "naca0014", "naca0016", "naca0020"):
        geom = c.GEOMDIR / f"{name}_dx0.0200.geom"
        for impl in IMPLS:
            jobs.append(("C2_thickness_family", name, geom, 300, 150, 0.01, 3000, impl, RE1000, ALPHA0, None))

    # ---------------- D1: LE point-density levels ----------------
    for factor in (0.5, 2, 8, 16):
        geom = c.GEOMDIR / f"naca0012_dx0.0200_LEdensity{factor:g}x.geom"
        for impl in IMPLS:
            jobs.append(("D1_point_density", f"LEdensity{factor:g}x", geom, 300, 150, 0.01, 3000,
                         impl, RE1000, ALPHA0, None))

    # ---------------- E1: decoupled front/TE variants ----------------
    for name in ("naca0012_dx0.0200_frontsharp_TEnative", "naca0012_dx0.0200_frontblunt_TEnative",
                 "naca0012_dx0.0200_TEsharp_frontnative"):
        geom = c.GEOMDIR / f"{name}.geom"
        for impl in IMPLS:
            jobs.append(("E1_decouple", name, geom, 300, 150, 0.01, 3000, impl, RE1000, ALPHA0, None))

    # ---------------- E2: common-TE resweep ----------------
    for name in ("naca0006_dx0.0200_roundTE", "naca0018_dx0.0200_roundTE"):
        geom = c.GEOMDIR / f"{name}.geom"
        for impl in IMPLS:
            jobs.append(("E2_common_TE", name, geom, 300, 150, 0.01, 3000, impl, RE1000, ALPHA0, None))

    return jobs


def main():
    njobs = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    jobs = build_jobs()
    print(f"{len(jobs)} jobs, {njobs}-way parallel", flush=True)
    done = fail = skip = 0
    with ProcessPoolExecutor(max_workers=njobs) as ex:
        futs = [ex.submit(job, *j) for j in jobs]
        for fut in as_completed(futs):
            name, status, elapsed = fut.result()
            if status == "skip":
                skip += 1
            elif status == "ok":
                done += 1
            else:
                fail += 1
                print(f"  FAIL: {name}", flush=True)
            n = done + fail + skip
            print(f"[{n}/{len(jobs)}] {name}: {status} ({elapsed:.0f}s)", flush=True)
    print(f"DONE. ran={done} skipped={skip} failed={fail}", flush=True)


if __name__ == "__main__":
    main()
