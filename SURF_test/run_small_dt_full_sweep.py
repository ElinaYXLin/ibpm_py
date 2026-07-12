"""
run_small_dt_full_sweep.py

Extends SURF_test/run_small_dt_validation.py's single-case check into the
FULL diagram set that lives in SD7003/SD8000's 2-c++included/ -- polar
sweep (5 alphas), grid convergence (3 dx levels), and the flowfield
snapshot case -- but all run at dt=0.001 (10x smaller than the dt=0.01 used
for 2-c++included/), for both py/ibpm.py (now native-FFTW3,
FFTW_EXHAUSTIVE -- see py/elliptic_solver_2d.py) and C++ build/ibpm, so
gen_small_dt_full_report.py can produce the same diagram types --
polar_comparison.png, drag_polar.png, grid_convergence.png,
flow_evolution.png -- for 3-small_dt/ that 2-c++included/ has for the
normal dt.

Time-averaging window is shorter than 2-c++included/'s (t=6 instead of
t=30): at dt=0.001 this is already 6000 steps per run (vs. 3000 at the
normal dt=0.01 for the same physical time), and FFTW_EXHAUSTIVE's planning
cost (paid once per elliptic-solver instance, 4 per Python run) is real --
see run_small_dt_validation.py's docstring. t=6 is well past the initial
impulsive-start transient (compare flow_evolution.png's t=0..30 montage:
the wake is already established within the first few time units) while
keeping this sweep's total wall-clock time tractable.

Usage:  python3 SURF_test/run_small_dt_full_sweep.py
Output: SURF_test/airfoils/LSAT-{SD7003,SD8000}/_run_data_smalldt_full/       (Python)
        SURF_test/airfoils/LSAT-{SD7003,SD8000}/_run_data_smalldt_full_cpp/   (C++)
        SURF_test/small_dt_full_results.json
"""
import json
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
GEOMDIR = REPO / "SURF_test" / "geom"
CPP_BIN = REPO / "build" / "ibpm"
DOMAIN = dict(length=6, xoffset=-2, yoffset=-1.5)

DT = 0.001            # 10x smaller than this suite's usual dt=0.01
NSTEPS = 6000           # t = 6.0
AVG_FRAC = 0.6           # average over the last 60% of the run (t=2.4-6.0), same convention as elsewhere
RESTART_EVERY_FLOW = 1000  # 7 snapshots (0,1000,...,6000) for flow_evolution.png

CASES = {
    "SD7003": dict(Re=61100, polar_alphas=[-2.92, -0.09, 1.66, 4.60, 7.72],
                    conv_alpha=-0.09, flow_alpha=4.60),
    "SD8000": dict(Re=60800, polar_alphas=[-3.88, -0.81, 2.29, 5.36, 8.36],
                    conv_alpha=-0.81, flow_alpha=5.36),
}
CONV_LEVELS = [
    dict(tag="coarse", dx=0.04, nx=150, ny=75),
    dict(tag="medium", dx=0.02, nx=300, ny=150),
    dict(tag="fine", dx=0.01, nx=600, ny=300),
]
POLAR_LEVEL = dict(dx=0.02, nx=300, ny=150)

RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"
RUNNER.write_text(
    "import sys, types, pathlib\n"
    "repo_root = pathlib.Path(%r)\n" % str(REPO) +
    "sys.path.insert(0, str(repo_root))\n"
    "pkg = types.ModuleType('py')\n"
    "pkg.__path__ = [str(repo_root / 'py')]\n"
    "sys.modules['py'] = pkg\n"
    "from py.ibpm import main\n"
    "sys.exit(main(['py.ibpm'] + sys.argv[1:]))\n"
)


def _args(geom, name, outdir, alpha, Re, nx, ny, restart=0):
    return [
        "-geom", str(geom), "-name", name, "-outdir", str(outdir),
        "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]), "-alpha", str(alpha), "-Re", str(Re),
        "-dt", str(DT), "-nsteps", str(NSTEPS), "-tecplot", "0",
        "-restart", str(restart), "-force", "1",
    ]


def run_case(cmd_prefix, geom, name, outdir, alpha, Re, nx, ny, restart=0):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = cmd_prefix + _args(geom, name, outdir, alpha, Re, nx, ny, restart)
    log_path = outdir / f"{name}_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"run failed ({outdir}): see {log_path}")
    return elapsed


def time_avg_force(force_path, frac=AVG_FRAC):
    import numpy as np
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    n = len(d)
    i0 = int(n * (1 - frac))
    seg = d[i0:]
    return dict(cl_mean=float(seg[:, 3].mean()), cl_std=float(seg[:, 3].std()),
                cd_mean=float(seg[:, 2].mean()), cd_std=float(seg[:, 2].std()))


def geom_for_dx(name, dx):
    p = GEOMDIR / f"{name.lower()}_dx{dx:.4f}.geom"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing -- run SURF_test/run_all_airfoils.py first")
    return p


def main():
    if not CPP_BIN.exists():
        print(f"ERROR: {CPP_BIN} not found. Build it first: cd build && make", file=sys.stderr)
        sys.exit(1)

    results_path = REPO / "SURF_test" / "small_dt_full_results.json"
    results = {"polar": {}, "convergence": {}}
    if results_path.exists():
        results = json.loads(results_path.read_text())

    def save():
        results_path.write_text(json.dumps(results, indent=2))

    for name, cfg in CASES.items():
        results["polar"].setdefault(name, {"py": [], "cpp": []})
        results["convergence"].setdefault(name, {"py": [], "cpp": []})
        done_py_alphas = {r["alpha"] for r in results["polar"][name]["py"]}
        done_cpp_alphas = {r["alpha"] for r in results["polar"][name]["cpp"]}

        geom = geom_for_dx(name, POLAR_LEVEL["dx"])

        # ---- polar sweep ----
        for alpha in cfg["polar_alphas"]:
            if alpha not in done_py_alphas:
                outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt_full" / f"polar_a{alpha:+.2f}"
                t = run_case([sys.executable, "-u", str(RUNNER)], geom, "run", outdir, alpha,
                             cfg["Re"], POLAR_LEVEL["nx"], POLAR_LEVEL["ny"])
                stats = time_avg_force(outdir / "run.force")
                stats.update(alpha=alpha, elapsed=t)
                results["polar"][name]["py"].append(stats)
                print(f"[py]  {name} polar alpha={alpha:+.2f}  Cl={stats['cl_mean']:+.4f}  "
                      f"Cd={stats['cd_mean']:+.4f}  ({t:.0f}s)", flush=True)
                save()
            if alpha not in done_cpp_alphas:
                outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt_full_cpp" / f"polar_a{alpha:+.2f}"
                t = run_case([str(CPP_BIN)], geom, "run", outdir, alpha,
                             cfg["Re"], POLAR_LEVEL["nx"], POLAR_LEVEL["ny"])
                stats = time_avg_force(outdir / "run.force")
                stats.update(alpha=alpha, elapsed=t)
                results["polar"][name]["cpp"].append(stats)
                print(f"[cpp] {name} polar alpha={alpha:+.2f}  Cl={stats['cl_mean']:+.4f}  "
                      f"Cd={stats['cd_mean']:+.4f}  ({t:.0f}s)", flush=True)
                save()

        # ---- grid convergence ----
        done_py_tags = {r["tag"] for r in results["convergence"][name]["py"]}
        done_cpp_tags = {r["tag"] for r in results["convergence"][name]["cpp"]}
        for lvl in CONV_LEVELS:
            geom_c = geom_for_dx(name, lvl["dx"])
            if lvl["tag"] not in done_py_tags:
                outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt_full" / f"conv_{lvl['tag']}"
                t = run_case([sys.executable, "-u", str(RUNNER)], geom_c, "run", outdir,
                             cfg["conv_alpha"], cfg["Re"], lvl["nx"], lvl["ny"])
                stats = time_avg_force(outdir / "run.force")
                stats.update(dx=lvl["dx"], tag=lvl["tag"], alpha=cfg["conv_alpha"], elapsed=t)
                results["convergence"][name]["py"].append(stats)
                print(f"[py]  {name} conv {lvl['tag']} dx={lvl['dx']}  Cl={stats['cl_mean']:+.4f}  "
                      f"Cd={stats['cd_mean']:+.4f}  ({t:.0f}s)", flush=True)
                save()
            if lvl["tag"] not in done_cpp_tags:
                outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt_full_cpp" / f"conv_{lvl['tag']}"
                t = run_case([str(CPP_BIN)], geom_c, "run", outdir,
                             cfg["conv_alpha"], cfg["Re"], lvl["nx"], lvl["ny"])
                stats = time_avg_force(outdir / "run.force")
                stats.update(dx=lvl["dx"], tag=lvl["tag"], alpha=cfg["conv_alpha"], elapsed=t)
                results["convergence"][name]["cpp"].append(stats)
                print(f"[cpp] {name} conv {lvl['tag']} dx={lvl['dx']}  Cl={stats['cl_mean']:+.4f}  "
                      f"Cd={stats['cd_mean']:+.4f}  ({t:.0f}s)", flush=True)
                save()

        # ---- flowfield (vorticity snapshots) ----
        flow_py_dir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt_full" / "flowfield"
        flow_cpp_dir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}" / "_run_data_smalldt_full_cpp" / "flowfield"
        if not (flow_py_dir / "run.force").exists():
            t = run_case([sys.executable, "-u", str(RUNNER)], geom, "run", flow_py_dir,
                         cfg["flow_alpha"], cfg["Re"], POLAR_LEVEL["nx"], POLAR_LEVEL["ny"],
                         restart=RESTART_EVERY_FLOW)
            print(f"[py]  {name} flowfield done ({t:.0f}s)", flush=True)
        if not (flow_cpp_dir / "run.force").exists():
            t = run_case([str(CPP_BIN)], geom, "run", flow_cpp_dir,
                         cfg["flow_alpha"], cfg["Re"], POLAR_LEVEL["nx"], POLAR_LEVEL["ny"],
                         restart=RESTART_EVERY_FLOW)
            print(f"[cpp] {name} flowfield done ({t:.0f}s)", flush=True)

    print("ALL DONE")


if __name__ == "__main__":
    main()
