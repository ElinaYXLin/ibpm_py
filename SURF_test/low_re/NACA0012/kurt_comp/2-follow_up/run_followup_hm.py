"""
run_followup_hm.py

New simulations for Tests H-M (see README.md's "H-M: why is Cd's excess
bigger than Cl's?" section) -- all follow-ups to fig13_14_hysteresis.png's
unresolved Cd-vs-paper discrepancy. Tests H and M need no new runs (H is a
cross-check of already-existing mean-coefficient data; M cross-checks the
paper's own two tables against each other). Tests I/J/K/L each need ONE
new f=4Hz, alpha0=0 pitching run per implementation (py_static AND
cpp_static, per mentor request -- see 5-leading_edge's Group 2/3 for the
same convention), varying exactly one knob from ../1-paper_based's
existing f4hz_{py,cpp}_a00 baseline:

  I. ngrid=2,3 (domain/blockage -- the same knob that fixed Cl's slope in
     Test D, here checked against the DYNAMIC Cd hysteresis for the first
     time)
  J. LE+TE boundary points refined to ds=dx/4 (the naca0012_dx0.0200_LTEdense
     geometry from ../5-leading_edge/, made pitching by make_hm_geoms.py)
  K. dt refined 0.005->0.0025 (same duration, double the steps)
  L. Re nudged +1% (1000->1010)

Usage: python3 run_followup_hm.py [py|cpp|both]   (default: both)
Output: runs/hm/<test>/ (py), runs/hm/<test>_cpp/ (cpp)
"""
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
HERE = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "2-follow_up"
RUNS = HERE / "runs" / "hm"
GEOMDIR = HERE / "geom"
PY_RUNNER = REPO / "static_test" / "run_ibpm_case_static.py"
CPP_BIN = REPO / "build_static" / "ibpm"
# NOTE: must use the PITCHING (f4hz) geometry variant for every one of these
# runs -- naca0012_dx0.0200.geom (no "motion pitchplunge" line) is the plain
# STEADY geometry and produces a non-oscillating run, which silently breaks
# the hysteresis-amplitude analysis (constant Cl/Cd -> zero peak-to-peak).
GEOM_020_F4HZ = (REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" /
                  "1-paper_based" / "geom" / "naca0012_dx0.0200_f4hz.geom")
LTEDENSE_F4HZ_GEOM = GEOMDIR / "naca0012_dx0.0200_LTEdense_f4hz.geom"

BASE_DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)
RE = 1000
ALPHA0 = 0.0


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


def outdir_for(name, impl):
    # py kept its original (pre-cpp-request) unsuffixed dir name so the
    # already-completed/in-progress py runs aren't disturbed; cpp is new,
    # so it gets an explicit suffix -- same convention as 5-leading_edge's
    # test2/test3 scripts.
    return RUNS / (name if impl == "py" else f"{name}_cpp")


def run_one(name, geom, nx, ny, ngrid, domain, dt, nsteps, re=RE, impl="py"):
    outdir = outdir_for(name, impl)
    if is_done(outdir, nsteps):
        print(f"{name} {impl}: already done", flush=True)
        return
    outdir.mkdir(parents=True, exist_ok=True)
    common = [
        "-geom", str(geom), "-name", "flow", "-outdir", str(outdir),
        "-nx", str(nx), "-ny", str(ny), "-ngrid", str(ngrid),
        "-length", str(domain["length"]), "-xoffset", str(domain["xoffset"]),
        "-yoffset", str(domain["yoffset"]), "-alpha", str(ALPHA0), "-Re", str(re),
        "-dt", str(dt), "-nsteps", str(nsteps), "-tecplot", "0",
        "-restart", "0", "-force", "1",
    ]
    cmd = [str(CPP_BIN)] + common if impl == "cpp" else [sys.executable, "-u", str(PY_RUNNER)] + common
    print(f"{name} {impl}: launching (nx={nx} ny={ny} ngrid={ngrid} dt={dt} nsteps={nsteps} Re={re})", flush=True)
    t0 = time.time()
    with open(outdir / "run_log.txt", "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    ok = proc.returncode == 0 and is_done(outdir, nsteps)
    print(f"  {'OK' if ok else 'FAILED'} in {elapsed:.0f}s", flush=True)


def main(impls=("py", "cpp")):
    ngrid_domain = dict(length=6, xoffset=-2, yoffset=-1.52)
    for impl in impls:
        # I: ngrid=2,3 (same ny=152 divisibility fix as Test D)
        for ngrid in (2, 3):
            run_one(f"I_ngrid{ngrid}", GEOM_020_F4HZ, 300, 152, ngrid, ngrid_domain,
                     dt=0.005, nsteps=6000, impl=impl)

        # J: LE+TE boundary points refined (background grid unchanged)
        run_one("J_LTEdense", LTEDENSE_F4HZ_GEOM, 300, 150, 1, BASE_DOMAIN,
                 dt=0.005, nsteps=6000, impl=impl)

        # K: dt refined 0.005 -> 0.0025, same duration (t=30)
        run_one("K_dt0.0025", GEOM_020_F4HZ, 300, 150, 1, BASE_DOMAIN,
                 dt=0.0025, nsteps=12000, impl=impl)

        # L: Re nudged +1%
        run_one("L_Re1010", GEOM_020_F4HZ, 300, 150, 1, BASE_DOMAIN,
                 dt=0.005, nsteps=6000, re=RE * 1.01, impl=impl)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "both"
    impls = {"py": ("py",), "cpp": ("cpp",), "both": ("py", "cpp")}[mode]
    main(impls)
