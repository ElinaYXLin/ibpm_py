"""
run_faithful.py

One case, matched to Kurtulus (2019)'s actual numerical setup as closely
as this solver's architecture allows: steady (non-oscillating) NACA0012,
alpha_0=12deg, Re=1000.

What's matched exactly:
  - Domain: x in [-15, 19] (34c: "far field boundary is located at 15c
    upstream and 19c downstream from the leading edge", paper p.3-4),
    y in [-15, 15] (30c -- paper doesn't state the lateral extent
    explicitly; inferred as the same 15c radius as upstream, standard for
    a C-type mesh's semicircular front/side boundary).
  - Duration and averaging window: paper runs steady cases to t=100s,
    dimensional, averaging over 50s<=t<=100s. Converting with this
    solver's own non-dimensionalization (c=0.1m, U_inf=0.146m/s, same
    conversion already used for the pitch frequencies elsewhere in
    kurt_comp): t_nondim = t_dim * U_inf/c = 100 * 0.146/0.1 = 146.
    nsteps=14600 at dt=0.01, averaging over the last 50% (t=73-146) --
    same FRACTION as before, ~5x the ABSOLUTE window.

What's deliberately NOT matched, and why: the paper's mesh is unstructured
near the body (first wall-normal cell 0.0015c) grading into a C-type
structured far-field -- not a single uniform dx. This solver is uniform
Cartesian only; matching 0.0015c across a 34c domain would need ~20,000+
cells per direction (days, not a day). Kept dx=0.02 (already this repo's
production resolution for NACA0012) as the stated, honest compromise --
matching what IS directly matchable (domain extent, run duration) while
being explicit about what isn't.

Usage: python3 run_faithful.py [py|cpp|both]
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
GEOM = REPO / "SURF_test" / "geom" / "naca0012_dx0.0200.geom"

# Grid.py uses a single dx = length/nx shared by both axes (no independent
# y-spacing) -- x in [xoffset, xoffset+length] = [-15, 19] (34c, "length"
# is the x-extent only), y in [yoffset, yoffset+ny*dx] = [-15, 15] (30c,
# since ny*dx = 1500*0.02 = 30). Both give dx=dy=0.02.
LENGTH_X, XOFFSET, YOFFSET = 34, -15, -15
NX, NY = 1700, 1500  # dx = 34/1700 = 0.02; y-extent = 1500*0.02 = 30 -> y in [-15,15]
RE = 1000
ALPHA = 12.0
DT = 0.01
NSTEPS = 14600  # t=146, matching the paper's steady run duration
RESTART = 1460  # 10 snapshots evenly through the run (paper reports at t=100s "instantaneous",
                 # mean over t=73-146 -- restart cadence just needs to bracket that window)


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
        "-nx", str(NX), "-ny", str(NY), "-ngrid", "1",
        "-length", "34", "-xoffset", "-15", "-yoffset", "-15",
        "-alpha", str(ALPHA), "-Re", str(RE), "-dt", str(DT), "-nsteps", str(NSTEPS),
        "-tecplot", "0", "-restart", str(RESTART), "-force", "1",
    ]
    if impl == "cpp":
        cmd = [str(CPP_BIN)] + common
    else:
        cmd = [sys.executable, "-u", str(PY_RUNNER)] + common
    print(f"{impl}: launching ({NX}x{NY} grid, {NSTEPS} steps) ...", flush=True)
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
