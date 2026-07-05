"""Driver: run py/ibpm.py on SD7003/SD8000 airfoils for a range of angle
of attack (polar comparison) and grid resolution (grid convergence study),
recording time-averaged Cl/Cd from each run's .force file.
"""
import sys, types, pathlib, subprocess, csv, time, json

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
RUNNER = REPO / "SURF_test" / "run_ibpm_case.py"

# Write out a tiny standalone runner (avoids re-creating the sys.modules
# shadow-fix workaround inline in every subprocess.run call).
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


def run_case(geom, name, outdir, alpha, nx, ny, length, xoffset, yoffset, dt, nsteps, Re,
             restart=0, force_every=1, extra_args=None):
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(RUNNER),
           "-geom", str(geom), "-name", name, "-outdir", str(outdir),
           "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
           "-length", str(length), "-xoffset", str(xoffset), "-yoffset", str(yoffset),
           "-alpha", str(alpha), "-Re", str(Re), "-dt", str(dt), "-nsteps", str(nsteps),
           "-tecplot", "0", "-restart", str(restart), "-force", str(force_every)]
    if extra_args:
        cmd += extra_args
    log_path = outdir / f"{name}_log.txt"
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"run failed ({name}): see {log_path}")
    return outdir / f"{name}.force", elapsed


def time_avg_force(force_path, frac=0.4):
    """Average Cd, Cl over the last `frac` of the run (skip transient)."""
    import numpy as np
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    n = len(d)
    i0 = int(n * (1 - frac))
    seg = d[i0:]
    t, cd, cl = d[:, 1], d[:, 2], d[:, 3]
    cd_mean, cl_mean = seg[:, 2].mean(), seg[:, 3].mean()
    cd_std, cl_std = seg[:, 2].std(), seg[:, 3].std()
    return dict(cd_mean=float(cd_mean), cl_mean=float(cl_mean),
                cd_std=float(cd_std), cl_std=float(cl_std),
                t_final=float(t[-1]), n_steps=int(n))


if __name__ == "__main__":
    print("driver module ready")
