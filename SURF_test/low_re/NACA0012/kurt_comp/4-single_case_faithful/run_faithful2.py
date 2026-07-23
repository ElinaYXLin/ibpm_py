"""
run_faithful2.py

Second attempt at matching Kurtulus (2019)'s actual numerical setup:
steady (non-oscillating) NACA0012, alpha_0=12deg, Re=1000.

What changed vs. old/run_faithful.py: that run kept dx=0.02 uniform
everywhere (this repo's production resolution) because matching the
paper's near-wall spacing (0.0015c) uniformly across the paper's full
34c x 30c domain would need ~20,000+ cells per direction -- estimated at
over a year of wall-clock time, infeasible.

This run instead uses the solver's multi-domain nesting (ngrid>1, each
level doubling dx and physical extent outward, same nx x ny cell count
at every level -- see py/grid.py) to get dx=0.0015c matched AT THE BODY,
while still reaching a 34c x 30c-scale far field via nesting.

The finest box is deliberately short-and-wide (rectangular, ny < nx):
this solver's grids don't have to be square (dx=length/nx is the only
spacing, shared by x and y, but nx and ny are independent), and the
per-level cost scales with nx*ny, so a thin airfoil doesn't need as
much vertical margin as chordwise margin -- shrinking ny cuts cost
directly without shrinking the region that actually needs resolving.

  - Finest level (lev=0): dx=0.0015c exactly, box x in [-0.3, 1.704]
    (2.004c: 0.3c upstream of the leading edge at x=0, 0.704c downstream
    of the trailing edge at x=1 -- enough to hold the near-wake vortex
    formation region before it crosses to the coarser level), y in
    [-0.3, 0.3] (0.6c, symmetric about the chord line -- ~5x the
    airfoil's own half-thickness of margin above/below, tighter than a
    square box but chosen deliberately to cut cost; this is the one
    parameter most likely to need revisiting if the near-body vortices
    turn out to be clipped by this margin).
  - ngrid=7 (one more than a square-box version would need, since a
    shorter finest box must nest one level further to reach the same
    far-field extent): outermost (lev=6) box is 2.004 * 2^6 = 128.3c (x)
    x 0.6 * 2^6 = 38.4c (y) -- exceeds the paper's 34c x 30c far field
    in both directions, comfortably (x is generously oversized, a side
    effect of prioritizing y-extent; extra levels are cheap, ~1/7 more
    elliptic-solve cost each, so this isn't worth trimming).
  - No xshift/yshift: nesting is centered on the finest box's own center
    (~0.7c, not exactly the leading edge), giving deep, roughly
    symmetric margins upstream and downstream at the outer levels.

Timestep: dt must shrink with the finest dx for CFL stability (the
scheme's advection term is explicit). Baseline (old run): dt=0.01,
dx=0.02 -> Courant number U*dt/dx = 0.5. Preserving that same Courant
number at dx=0.0015: dt = 0.5*0.0015/1.0 = 0.00075.

Duration: same non-dimensional t=146 as the old run (t_dim=100s,
c=0.1m, U_inf=0.146m/s -- see old/run_faithful.py's derivation),
averaging over the last 50% (t=73-146). nsteps = 146/0.00075 = 194667
(t=146.00025, close enough).

Estimated wall-clock: ~3.3 days per implementation (see chat discussion;
extrapolated from the old run's measured 4.58h/14600 steps at
1700x1500/ngrid=1, scaled for cell count, ngrid, log-factor, and the
~13.3x more steps from the smaller dt -- not independently benchmarked
at this resolution, so treat as order-of-magnitude, +/-~2x).

Usage: python3 run_faithful2.py [py|cpp|both]
Output: 4-single_case_faithful/runs/<impl>/
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "4-single_case_faithful"
RUNS = HERE / "runs"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"
CPP_BIN = REPO / "build_static" / "ibpm"
GEOM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0015.geom"

# Finest level: dx = LENGTH_X/NX = 2.004/1336 = 0.0015 exactly.
# y-extent = NY*dx = 400*0.0015 = 0.6c. Short-and-wide box (ny < nx):
# the airfoil is thin, so less vertical margin is needed than chordwise
# margin -- see module docstring.
LENGTH_X, XOFFSET, YOFFSET = 2.004, -0.3, -0.3
NX, NY, NGRID = 1336, 400, 7
RE = 1000
ALPHA = 12.0
DT = 0.00075
NSTEPS = 194667  # t = 194667*0.00075 = 146.00025, matching the paper's steady-run duration
RESTART = NSTEPS // 10  # ~10 snapshots evenly through the run, same convention as old run


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


def run_one(impl):
    outdir = RUNS / impl
    if is_done(outdir, NSTEPS):
        print(f"{impl}: already done, skipping")
        return True
    outdir.mkdir(parents=True, exist_ok=True)
    common = [
        "-geom", str(GEOM), "-name", "flow", "-outdir", str(outdir),
        "-nx", str(NX), "-ny", str(NY), "-ngrid", str(NGRID),
        "-length", str(LENGTH_X), "-xoffset", str(XOFFSET), "-yoffset", str(YOFFSET),
        "-alpha", str(ALPHA), "-Re", str(RE), "-dt", str(DT), "-nsteps", str(NSTEPS),
        "-tecplot", "0", "-restart", str(RESTART), "-force", "1",
    ]
    if impl == "cpp":
        cmd = [str(CPP_BIN)] + common
    else:
        cmd = [sys.executable, "-u", str(PY_RUNNER)] + common
    print(f"{impl}: launching ({NX}x{NY} grid, ngrid={NGRID}, {NSTEPS} steps) ...", flush=True)
    t0 = time.time()
    with open(outdir / "run_log.txt", "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    dt_ = time.time() - t0
    ok = proc.returncode == 0 and is_done(outdir, NSTEPS)
    print(f"{impl}: {'OK' if ok else 'FAILED'} in {dt_:.0f}s ({dt_/3600:.2f}h)")
    return ok


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    results = {}
    if which in ("py", "both"):
        results["py"] = run_one("py")
    if which in ("cpp", "both"):
        results["cpp"] = run_one("cpp")
    print("SUMMARY:", results)
    sys.exit(0 if all(results.values()) else 1)


if __name__ == "__main__":
    main()
