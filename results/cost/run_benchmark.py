# run_benchmark.py
#
# Computational cost analysis: runs each *implemented* backend in
# backends.BACKENDS across a range of grid resolutions, on the same
# cylinder-at-Re=100 case, and records:
#   - wall-clock runtime (total, and split into model-construction /
#     solver-factorization / timestepping phases, using the same
#     stdout-marker technique as results/generate_validation_report.py --
#     see that script's module docstring for why "total wall time" alone
#     is a misleading cost comparison for this codebase)
#   - CPU time (user+sys seconds), via resource.getrusage(RUSAGE_CHILDREN)
#     (authoritative -- computed at process exit, no polling/race issues)
#   - peak RSS memory, both via getrusage (authoritative) and via a
#     psutil polling thread (used only for the RSS/CPU% *time series*,
#     which getrusage cannot provide)
#
# Usage: python3 results/cost/run_benchmark.py
# (run from the repository root; requires build/ibpm to exist -- see
# build/Makefile's docstring/README for how it was reconstructed)

from __future__ import annotations

import json
import os
import pathlib
import resource
import subprocess
import sys
import threading
import time
from typing import Dict, List, Tuple

import psutil

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from backends import BACKENDS, REPO_ROOT, make_circle_geom

COST_DIR = pathlib.Path(__file__).resolve().parent
GEOM_DIR = COST_DIR / "geom"
RAW_DIR = COST_DIR / "raw"
GEOM_DIR.mkdir(exist_ok=True)
RAW_DIR.mkdir(exist_ok=True)

MARKER_SOLVER_START = "solver for projection step"
MARKER_INTEGRATION_START = "Integrating for"
MARKER_STEP_PREFIX = "step "

# Grid resolutions to benchmark (nx=ny; single-domain, ngrid=1). dt=0.01 is
# used at nx=400 instead of the default 0.02 because results/README.md
# already documented that this exact case (400x400, dt=0.02) blows up to
# NaN in BOTH C++ and Python -- a genuine CFL limit, not implementation-
# specific -- so 0.01 (stable) is used there instead, exactly as that
# report did.
GRID_SIZES = [100, 200, 300, 400]
NSTEPS = 150
SAMPLE_INTERVAL = 0.05  # seconds, psutil polling cadence for time-series


def dt_for(nx: int) -> float:
    return 0.01 if nx >= 400 else 0.02


class ResourceMonitor:
    """Polls a subprocess's (and its children's) RSS/CPU% on a background
    thread while it runs, for time-series plots. Not used for the
    authoritative summary numbers (see getrusage in run_one)."""

    def __init__(self, pid: int, interval: float = SAMPLE_INTERVAL):
        self.pid = pid
        self.interval = interval
        self.samples: List[Tuple[float, float, float]] = []  # (t, rss_bytes, cpu_percent)
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        # NOTE: psutil.Process.cpu_percent() tracks elapsed CPU time *since
        # the previous call on that same Process object* -- it always
        # returns 0.0 on a freshly-constructed Process (no prior baseline).
        # Re-instantiating psutil.Process(pid) every poll (as an earlier
        # version of this code did) therefore always reports 0.0. Process
        # objects are cached here instead, one per (pid), reused across
        # polls, with new children picked up (and dead ones dropped) each
        # iteration.
        self._proc_cache: Dict[int, psutil.Process] = {}

    def _refresh_procs(self):
        try:
            root = psutil.Process(self.pid)
        except psutil.NoSuchProcess:
            self._proc_cache = {}
            return []
        live_pids = {self.pid}
        try:
            live_pids |= {c.pid for c in root.children(recursive=True)}
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        for pid in list(self._proc_cache):
            if pid not in live_pids:
                del self._proc_cache[pid]
        for pid in live_pids:
            if pid not in self._proc_cache:
                try:
                    self._proc_cache[pid] = psutil.Process(pid)
                    self._proc_cache[pid].cpu_percent()  # prime baseline
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        return list(self._proc_cache.values())

    def _run(self):
        t0 = time.perf_counter()
        self._refresh_procs()  # prime baselines before the first real sample
        while not self._stop.is_set():
            time.sleep(self.interval)
            procs = self._refresh_procs()
            if not procs:
                break
            rss = 0
            cpu = 0.0
            for p in procs:
                try:
                    rss += p.memory_info().rss
                    cpu += p.cpu_percent()
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            self.samples.append((time.perf_counter() - t0, rss, cpu))

    def start(self):
        self._thread.start()

    def stop(self):
        self._stop.set()
        self._thread.join(timeout=2 * self.interval)


def run_one(cmd: List[str], outdir: pathlib.Path) -> Dict:
    outdir.mkdir(parents=True, exist_ok=True)
    ru_before = resource.getrusage(resource.RUSAGE_CHILDREN)

    proc = subprocess.Popen(cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, bufsize=1)
    monitor = ResourceMonitor(proc.pid)
    monitor.start()

    t0 = time.perf_counter()
    markers: Dict[str, float] = {}
    last_step_t = 0.0
    log_lines = []
    for line in proc.stdout:
        t = time.perf_counter() - t0
        log_lines.append(f"{t:8.3f}s  {line}")
        if MARKER_SOLVER_START in line and "solver_start" not in markers:
            markers["solver_start"] = t
        if MARKER_INTEGRATION_START in line and "integration_start" not in markers:
            markers["integration_start"] = t
        if line.strip().startswith(MARKER_STEP_PREFIX):
            last_step_t = t
    proc.wait()
    total_wall = time.perf_counter() - t0
    monitor.stop()

    ru_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    (outdir / "timestamped_log.txt").write_text("".join(log_lines))
    if proc.returncode != 0:
        raise RuntimeError(f"run failed (exit {proc.returncode}); see {outdir}/timestamped_log.txt")

    markers["integration_end"] = last_step_t if last_step_t > 0 else total_wall
    model_t = markers.get("solver_start", 0.0)
    setup_t = markers.get("integration_start", total_wall) - model_t
    steps_t = markers["integration_end"] - markers.get("integration_start", 0.0)

    # ru_maxrss is bytes on macOS/Darwin, KB on Linux -- normalize to bytes.
    maxrss_units = 1 if sys.platform == "darwin" else 1024
    peak_rss_bytes = (ru_after.ru_maxrss - 0) * maxrss_units  # ru_maxrss is already a peak, not cumulative
    cpu_user = ru_after.ru_utime - ru_before.ru_utime
    cpu_sys = ru_after.ru_stime - ru_before.ru_stime

    ts_peak_rss = max((s[1] for s in monitor.samples), default=0.0)

    return dict(
        wall_time=total_wall,
        phase_model=model_t,
        phase_setup=setup_t,
        phase_steps=steps_t,
        cpu_user=cpu_user,
        cpu_sys=cpu_sys,
        cpu_total=cpu_user + cpu_sys,
        peak_rss_bytes=peak_rss_bytes,
        peak_rss_bytes_sampled=ts_peak_rss,
        timeseries=[{"t": t, "rss_bytes": rss, "cpu_percent": cpu} for t, rss, cpu in monitor.samples],
    )


def main():
    results = []
    results_path = RAW_DIR / "cost_results.json"

    for backend in BACKENDS:
        if not backend.implemented:
            print(f"=== skipping backend '{backend.name}': {backend.note} ===", flush=True)
            continue
        for nx in GRID_SIZES:
            ny = nx
            dt = dt_for(nx)
            dx = 4.0 / nx
            geom_path = GEOM_DIR / f"cylinder_dx{dx:.5f}.geom"
            npts = make_circle_geom(dx, geom_path)

            name = f"{backend.name}_nx{nx}"
            outdir = RAW_DIR / name
            cmd = backend.build_cmd(geom_path, outdir, "run", nx, ny, dt, NSTEPS)

            print(f"=== {backend.name}: nx=ny={nx} dx={dx:.4f} dt={dt} npts={npts} ===", flush=True)
            t0 = time.time()
            r = run_one(cmd, outdir)
            elapsed = time.time() - t0
            r.update(backend=backend.name, label=backend.label, nx=nx, ny=ny, dx=dx, dt=dt,
                     nsteps=NSTEPS, npts=npts)
            results.append(r)
            print(f"    wall={r['wall_time']:.2f}s  (model={r['phase_model']:.2f}s "
                  f"setup={r['phase_setup']:.2f}s steps={r['phase_steps']:.2f}s)  "
                  f"cpu={r['cpu_total']:.2f}s  peakRSS={r['peak_rss_bytes']/1e6:.1f}MB  "
                  f"[{elapsed:.1f}s elapsed]", flush=True)

            results_path.write_text(json.dumps(results, indent=2))

    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
