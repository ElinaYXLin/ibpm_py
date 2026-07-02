#!/usr/bin/env python3
"""
generate_validation_report.py

Validates the Python port (py/) against the reference C++ implementation
(src/, compiled as build/ibpm) by running BOTH on every valid example
geometry in this repository and comparing:

  - Drag/lift coefficients (Cd, Cl) over a full timestepping run
  - Absolute force error vs. time (a dense, per-timestep error signal)
  - Vorticity field snapshots at t=0 and at the final step
  - Wall-clock runtime, broken down into three phases (see "A note on
    runtime", below) -- NOT just a single total-time number.

Every number and figure this script produces comes directly from freshly
re-running the two actual simulators (subprocess calls to build/ibpm and
`python3 -m py.ibpm`) and reading their binary/text output -- nothing here
is hand-entered or AI-generated.

Which geometries are run
-------------------------
This repository contains four `.geom` files. Three are used here; one is
excluded, with the reason recorded in EXCLUDED_CASES and printed at
runtime (not silently skipped):

  - examples/cylinder.geom          (160-point circle; the standard example)
  - ibpm.geom                       (314-point circle; a finer boundary,
                                      originally used to match a Fortran
                                      reference case)
  - benchmarking/cylinder2Pa.geom   (160-point circle with an explicit
                                      "motion fixed 0 0 0" command -- same
                                      body as cylinder.geom, but exercises
                                      RigidBody's motion-parsing code path,
                                      which cylinder.geom does not)
  - benchmarking/cylinder2PaPlunge.geom -- EXCLUDED. Its "motion PitchPlunge
    0 0 0.5 0.2" line only supplies 4 of PitchPlunge's 6 required
    parameters (a malformed file, not a porting artifact), AND PitchPlunge
    motion itself is not yet ported to Python (see py/fixed_position.py's
    module docstring for why). Geometry.load() returns False for this file
    in the Python port; it cannot currently be run for comparison.

A note on runtime
------------------
A naive "total wall-clock time" comparison is misleading for this specific
codebase: src/EllipticSolver2d.cc constructs its FFTW plan with
`FFTW_EXHAUSTIVE`, which does an expensive one-time *search* over FFT
algorithms the first time a transform of a given size is requested (this
is a property of the C++ build, unrelated to numerical correctness). That
one-time cost is paid once, during NavierStokesModel construction, and
dominates a cold-start run of this size -- it is NOT representative of the
per-timestep computational cost, and it makes the *whole-run* time
essentially incomparable between C++ and Python (whose scipy.fft has no
such planning phase at all). This script instead times three phases
separately, using the same progress-message markers both programs print,
so the misleading part (one-time FFTW planning) is visible on its own
axis instead of contaminating the "real" per-step comparison:

  1. "model construction"  -- process start -> first ProjectionSolver
     construction begins (dominated by FFTW_EXHAUSTIVE in C++)
  2. "solver factorization" -- building/factoring the projection operator
     for each timestepping substep
  3. "timestepping"        -- the actual N-step time-integration loop
     (the number that should be compared to judge computational
     efficiency of the ported numerics)

Usage:
    python3 results/generate_validation_report.py

Requirements:
    - numpy, matplotlib (see py/requirements.txt; matplotlib is only
      needed for this script, not for py/ itself)
    - a compiled C++ reference binary at build/ibpm (see the top-level
      README's "Installation" section: `cd build && make`). If it is not
      present, the script exits early with an explanation rather than
      fabricating comparison data.

Output (written under results/):
    figures/<case>/force_coefficients_vs_time.png
    figures/<case>/force_error_vs_time.png
    figures/<case>/vorticity_field_comparison.png
    figures/<case>/vorticity_parity_plot.png
    figures/<case>/flow_evolution_python.png
    figures/all_cases_force_comparison.png      (stitched, all cases)
    figures/all_cases_vorticity_comparison.png  (stitched, all cases)
    figures/runtime_phase_breakdown.png         (stitched, all cases)
    validation_metrics.csv
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))  # so `from py.state import State` resolves below

RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
RUN_DATA_DIR = RESULTS_DIR / "_run_data"  # raw run output, gitignored

CPP_BIN = REPO_ROOT / "build" / "ibpm"

# Default timestepping parameters; overridable per-case (see CASES below)
# since stability depends on grid resolution (a finer grid needs a smaller
# dt -- this is a genuine CFL/stability property of the explicit-viscous
# fractional-step scheme, not a porting concern; see the "ibpm_geom" case).
DT = 0.02
NSTEPS = 250
RE = 100.0
SCHEME = "rk3"
SNAPSHOT_STEPS = [0, 200]  # final snapshot kept at step 200 (not 250) so
# every case only needs 2 restart files (cheaper) while still showing a
# well-developed wake; -restart 200 below writes exactly these two.

CASES = [
    dict(name="cylinder", geom="examples/cylinder.geom", nx=200, ny=200,
         label="cylinder.geom (160 pts)", dt=DT, nsteps=NSTEPS),
    dict(name="ibpm_geom", geom="ibpm.geom", nx=400, ny=400,
         label="ibpm.geom (314 pts)",
         # dt=0.02 (the default used by the other cases) blows up to NaN
         # around t=0.3 on this finer 400x400 grid, IDENTICALLY in both
         # C++ and Python (verified: both diverge at the same step, with
         # matching diverging values) -- a genuine CFL/stability limit of
         # this explicit scheme at this resolution, not a port discrepancy.
         # dt=0.01 is stable for the full 250-step run.
         dt=0.01, nsteps=NSTEPS),
    dict(name="cylinder2pa", geom="benchmarking/cylinder2Pa.geom", nx=200, ny=200,
         label="cylinder2Pa.geom (160 pts, 'motion fixed')", dt=DT, nsteps=NSTEPS),
]

EXCLUDED_CASES = [
    dict(name="cylinder2paplunge", geom="benchmarking/cylinder2PaPlunge.geom",
         reason="'motion PitchPlunge 0 0 0.5 0.2' supplies only 4 of the 6 "
                "parameters PitchPlunge requires (a malformed file), and "
                "PitchPlunge motion is not yet ported to Python -- "
                "Geometry.load() returns False for this file in py/. "
                "Not run."),
]

# Markers (exact substrings) both programs print, used to split a run into
# phases. See "A note on runtime" above.
MARKER_SOLVER_START = "solver for projection step"  # "Using Cholesky/ConjugateGradient ..."
MARKER_INTEGRATION_START = "Integrating for"
MARKER_STEP_PREFIX = "step "


def _grid_args(nx: int, ny: int, dt: float, nsteps: int) -> List[str]:
    return [
        "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
        "-length", "4.0", "-xoffset", "-2.0", "-yoffset", "-2.0",
        "-Re", str(RE), "-dt", str(dt), "-nsteps", str(nsteps), "-scheme", SCHEME,
        "-tecplot", "0", "-restart", str(SNAPSHOT_STEPS[-1]), "-force", "1", "-energy", "0",
    ]


def _run_timed(cmd: List[str], outdir: Path, env: Dict[str, str]) -> Tuple[float, Dict[str, float]]:
    """Run `cmd`, streaming stdout+stderr with per-line timestamps (relative
    to process start), and return (total_wall_time, phase_markers).

    phase_markers maps marker name -> first timestamp (seconds) at which a
    line containing that marker text appeared.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.Popen(cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, bufsize=1, env=env)
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
    total = time.perf_counter() - t0
    markers["integration_end"] = last_step_t if last_step_t > 0 else total
    (outdir / "timestamped_log.txt").write_text("".join(log_lines))
    if proc.returncode != 0:
        raise RuntimeError(f"run failed (exit {proc.returncode}); see {outdir}/timestamped_log.txt")
    return total, markers


def run_cpp(geom: Path, nx: int, ny: int, dt: float, nsteps: int, outdir: Path) -> Tuple[float, Dict[str, float]]:
    cmd = [str(CPP_BIN), "-geom", str(geom), "-outdir", str(outdir), "-name", "ibpm",
           *_grid_args(nx, ny, dt, nsteps)]
    return _run_timed(cmd, outdir, env=dict(os.environ))


def run_python(geom: Path, nx: int, ny: int, dt: float, nsteps: int, outdir: Path) -> Tuple[float, Dict[str, float]]:
    cmd = [sys.executable, "-u", "-m", "py.ibpm", "-geom", str(geom), "-outdir", str(outdir),
           "-name", "ibpm", *_grid_args(nx, ny, dt, nsteps)]
    env = dict(os.environ, PYTHONUNBUFFERED="1")
    return _run_timed(cmd, outdir, env=env)


def load_omega(bin_path: Path) -> np.ndarray:
    """Load the finest-level vorticity field from an IBPM restart file,
    using the already-validated py.state.State reader (same binary format
    for both the C++ and Python outputs)."""
    from py.state import State

    s = State(filename=str(bin_path))
    return s.omega._data[0].copy()  # shape (nx-1, ny-1), finest grid level


def load_force(force_path: Path) -> np.ndarray:
    return np.loadtxt(force_path)


def phase_durations(markers: Dict[str, float]) -> Tuple[float, float, float]:
    model = markers["solver_start"]
    setup = markers["integration_start"] - markers["solver_start"]
    steps = markers["integration_end"] - markers["integration_start"]
    return model, setup, steps


def field_errors(cpp_f: np.ndarray, py_f: np.ndarray) -> Tuple[float, float, float, float]:
    diff = np.abs(py_f - cpp_f)
    peak = np.abs(cpp_f).max()
    if peak == 0:
        return diff.max(), diff.mean(), 0.0, 0.0
    return diff.max(), diff.mean(), diff.max() / peak, diff.mean() / peak


def run_case(case: dict) -> dict:
    name = case["name"]
    geom = REPO_ROOT / case["geom"]
    nx, ny = case["nx"], case["ny"]
    dt = case.get("dt", DT)
    nsteps = case.get("nsteps", NSTEPS)
    case_dir = FIGURES_DIR / name
    case_dir.mkdir(parents=True, exist_ok=True)
    cpp_dir = RUN_DATA_DIR / name / "cpp"
    py_dir = RUN_DATA_DIR / name / "python"

    print(f"[{name}] running C++ reference ({geom.name}, {nx}x{ny}, dt={dt}, nsteps={nsteps})...")
    cpp_total, cpp_markers = run_cpp(geom, nx, ny, dt, nsteps, cpp_dir)
    print(f"[{name}]   done in {cpp_total:.2f} s")

    print(f"[{name}] running Python port...")
    py_total, py_markers = run_python(geom, nx, ny, dt, nsteps, py_dir)
    print(f"[{name}]   done in {py_total:.2f} s")

    cpp_force = load_force(cpp_dir / "ibpm.force")
    py_force = load_force(py_dir / "ibpm.force")
    assert cpp_force.shape == py_force.shape, f"[{name}] step count mismatch"

    t = cpp_force[:, 1]
    cd_cpp, cl_cpp = cpp_force[:, 2], cpp_force[:, 3]
    cd_py, cl_py = py_force[:, 2], py_force[:, 3]
    cd_abs_err = np.abs(cd_py - cd_cpp)
    cl_abs_err = np.abs(cl_py - cl_cpp)
    nonzero = cd_cpp != 0
    cd_rel_err = cd_abs_err[nonzero] / np.abs(cd_cpp[nonzero])

    xs = -2.0 + np.arange(1, nx) * (4.0 / nx)
    ys = -2.0 + np.arange(1, ny) * (4.0 / ny)
    X, Y = np.meshgrid(xs, ys, indexing="ij")

    omega_cpp = {s: load_omega(cpp_dir / f"ibpm{s:05d}.bin") for s in SNAPSHOT_STEPS}
    omega_py = {s: load_omega(py_dir / f"ibpm{s:05d}.bin") for s in SNAPSHOT_STEPS}

    # ---- per-case figures ----
    _fig_force(case_dir, case, t, cd_cpp, cl_cpp, cd_py, cl_py, cd_abs_err, cl_abs_err)
    vmax = max(np.abs(omega_cpp[s]).max() for s in SNAPSHOT_STEPS) or 1.0
    _fig_vorticity_grid(case_dir, case, X, Y, omega_cpp, omega_py, vmax)
    _fig_parity(case_dir, case, omega_cpp[SNAPSHOT_STEPS[-1]], omega_py[SNAPSHOT_STEPS[-1]], nx, ny)
    _fig_flow_evolution(case_dir, case, X, Y, omega_py, vmax)

    return dict(
        case=case, dt=dt, nsteps=nsteps, t=t, cd_cpp=cd_cpp, cl_cpp=cl_cpp, cd_py=cd_py, cl_py=cl_py,
        cd_rel_err=cd_rel_err, cl_abs_err=cl_abs_err,
        omega_cpp=omega_cpp, omega_py=omega_py, X=X, Y=Y, vmax=vmax,
        cpp_total=cpp_total, py_total=py_total,
        cpp_phases=phase_durations(cpp_markers), py_phases=phase_durations(py_markers),
    )


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _fig_force(case_dir, case, t, cd_cpp, cl_cpp, cd_py, cl_py, cd_abs_err, cl_abs_err) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    axes[0].plot(t, cd_cpp, "-", color="tab:blue", lw=2.5, alpha=0.5, label="C++ (src/)")
    axes[0].plot(t, cd_py, "--", color="tab:orange", lw=1.5, label="Python (py/)")
    axes[0].set_ylabel(r"$C_d$")
    axes[0].legend()
    axes[0].set_title(f"{case['label']}: force coefficients vs. time")
    axes[1].plot(t, cl_cpp, "-", color="tab:blue", lw=2.5, alpha=0.5, label="C++ (src/)")
    axes[1].plot(t, cl_py, "--", color="tab:orange", lw=1.5, label="Python (py/)")
    axes[1].set_ylabel(r"$C_l$")
    axes[1].set_xlabel("time")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "force_coefficients_vs_time.png", dpi=150)
    plt.close(fig)

    floor = 1e-18
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    axes[0].semilogy(t, np.maximum(cd_abs_err, floor), color="tab:red")
    axes[0].set_ylabel(r"$|C_{d,\mathrm{py}} - C_{d,\mathrm{cpp}}|$")
    axes[0].set_title(f"{case['label']}: Python-vs-C++ absolute force error vs. time")
    if cd_abs_err.max() <= floor:
        axes[0].text(0.5, 0.5, "error is exactly 0.0 at every timestep\n(line at floor is a plotting artifact)",
                     transform=axes[0].transAxes, ha="center", va="center", fontsize=9, color="dimgray")
    axes[1].semilogy(t, np.maximum(cl_abs_err, floor), color="tab:red")
    axes[1].set_ylabel(r"$|C_{l,\mathrm{py}} - C_{l,\mathrm{cpp}}|$")
    axes[1].set_xlabel("time")
    fig.tight_layout()
    fig.savefig(case_dir / "force_error_vs_time.png", dpi=150)
    plt.close(fig)


def _fig_vorticity_grid(case_dir, case, X, Y, omega_cpp, omega_py, vmax) -> None:
    dt = case.get("dt", DT)
    fig, axes = plt.subplots(len(SNAPSHOT_STEPS), 3, figsize=(12, 3.3 * len(SNAPSHOT_STEPS)))
    if len(SNAPSHOT_STEPS) == 1:
        axes = axes[None, :]
    for row, step in enumerate(SNAPSHOT_STEPS):
        time_val = step * dt
        diff = omega_py[step] - omega_cpp[step]
        panels = [
            (omega_cpp[step], f"C++, t={time_val:g}", "RdBu_r", -vmax, vmax),
            (omega_py[step], f"Python, t={time_val:g}", "RdBu_r", -vmax, vmax),
            (diff, f"Python - C++ (max |diff| = {np.abs(diff).max():.1e})", "PuOr", None, None),
        ]
        for col, (field, title, cmap, vmin_, vmax_) in enumerate(panels):
            ax = axes[row, col]
            im = ax.contourf(X, Y, field, levels=41, cmap=cmap, vmin=vmin_, vmax=vmax_)
            ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, color="k", lw=1.2))
            ax.set_aspect("equal")
            ax.set_xlim(-2, 3)
            ax.set_ylim(-2, 2)
            ax.set_title(title, fontsize=10)
            fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle(f"{case['label']}: vorticity field, C++ vs. Python vs. difference")
    fig.tight_layout()
    fig.savefig(case_dir / "vorticity_field_comparison.png", dpi=150)
    plt.close(fig)


def _fig_parity(case_dir, case, cpp_field, py_field, nx, ny) -> None:
    a, b = cpp_field.ravel(), py_field.ravel()
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(a, b, s=2, alpha=0.3, color="tab:blue")
    lims = [min(a.min(), b.min()), max(a.max(), b.max())]
    ax.plot(lims, lims, "k--", lw=1, label="y = x")
    ax.set_xlabel(r"C++ vorticity $\omega$")
    ax.set_ylabel(r"Python vorticity $\omega$")
    r2 = np.corrcoef(a, b)[0, 1] ** 2 if np.ptp(a) > 0 else 1.0
    ax.set_title(f"{case['label']}: parity plot, all interior points\n"
                 f"({(nx - 1) * (ny - 1):,} points, $R^2$ = {r2:.10f})")
    ax.legend()
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(case_dir / "vorticity_parity_plot.png", dpi=150)
    plt.close(fig)


def _fig_flow_evolution(case_dir, case, X, Y, omega_py, vmax) -> None:
    dt = case.get("dt", DT)
    fig, axes = plt.subplots(1, len(SNAPSHOT_STEPS), figsize=(4.2 * len(SNAPSHOT_STEPS), 4.5),
                              constrained_layout=True)
    if len(SNAPSHOT_STEPS) == 1:
        axes = [axes]
    for ax, step in zip(axes, SNAPSHOT_STEPS):
        im = ax.contourf(X, Y, omega_py[step], levels=41, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, color="k", lw=1.2))
        ax.set_aspect("equal")
        ax.set_xlim(-2, 3)
        ax.set_ylim(-2, 2)
        ax.set_title(f"t = {step * dt:g}")
    fig.suptitle(f"{case['label']}: flow development (Python port, Re={RE:g})")
    fig.colorbar(im, ax=axes, shrink=0.7, label="vorticity")
    fig.savefig(case_dir / "flow_evolution_python.png", dpi=150)
    plt.close(fig)


def _fig_stitched_force(results: List[dict]) -> None:
    fig, axes = plt.subplots(len(results), 1, figsize=(7, 3 * len(results)), sharex=True)
    if len(results) == 1:
        axes = [axes]
    for ax, r in zip(axes, results):
        ax.plot(r["t"], r["cd_cpp"], "-", color="tab:blue", lw=2.5, alpha=0.5, label="C++")
        ax.plot(r["t"], r["cd_py"], "--", color="tab:orange", lw=1.5, label="Python")
        ax.set_ylabel(r"$C_d$")
        ax.set_title(r["case"]["label"], fontsize=10)
        ax.legend(fontsize=8)
    axes[-1].set_xlabel("time")
    fig.suptitle("All cases: drag coefficient, C++ vs. Python")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "all_cases_force_comparison.png", dpi=150)
    plt.close(fig)


def _fig_stitched_vorticity(results: List[dict]) -> None:
    """Stitch each case's final-snapshot {C++, Python, diff} row into one
    combined figure -- one row per case."""
    final_step = SNAPSHOT_STEPS[-1]
    fig, axes = plt.subplots(len(results), 3, figsize=(12, 3.3 * len(results)))
    if len(results) == 1:
        axes = axes[None, :]
    for row, r in zip(range(len(results)), results):
        vmax = r["vmax"]
        dt = r["dt"]
        cpp_f = r["omega_cpp"][final_step]
        py_f = r["omega_py"][final_step]
        diff = py_f - cpp_f
        panels = [
            (cpp_f, f"C++, {r['case']['label']}, t={final_step * dt:g}", "RdBu_r", -vmax, vmax),
            (py_f, f"Python, {r['case']['label']}, t={final_step * dt:g}", "RdBu_r", -vmax, vmax),
            (diff, f"diff (max |diff| = {np.abs(diff).max():.1e})", "PuOr", None, None),
        ]
        for col, (field, title, cmap, vmin_, vmax_) in enumerate(panels):
            ax = axes[row, col]
            im = ax.contourf(r["X"], r["Y"], field, levels=41, cmap=cmap, vmin=vmin_, vmax=vmax_)
            ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, color="k", lw=1.2))
            ax.set_aspect("equal")
            ax.set_xlim(-2, 3)
            ax.set_ylim(-2, 2)
            ax.set_title(title, fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("All cases: vorticity field at final snapshot, C++ vs. Python vs. difference\n"
                 "(note: per-case final time differs -- see subplot titles -- since dt is tuned per case)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "all_cases_vorticity_comparison.png", dpi=150)
    plt.close(fig)


def _fig_runtime_phases(results: List[dict]) -> None:
    """Stacked bar chart: for each case, one C++ bar and one Python bar,
    each stacked into the 3 phases described in the module docstring."""
    labels = ["model\nconstruction", "solver\nfactorization", "timestepping\n(250 steps)"]
    colors = ["#c9c9c9", "#8ecae6", "#219ebc"]
    n = len(results)
    fig, ax = plt.subplots(figsize=(2.2 * n + 2, 5.5))
    x = np.arange(n)
    width = 0.35
    for impl, offset, hatch in [("cpp_phases", -width / 2, None), ("py_phases", width / 2, "//")]:
        bottoms = np.zeros(n)
        for phase_idx, (label, color) in enumerate(zip(labels, colors)):
            vals = np.array([r[impl][phase_idx] for r in results])
            ax.bar(x + offset, vals, width, bottom=bottoms, color=color, hatch=hatch,
                   edgecolor="black", linewidth=0.5,
                   label=label if offset < 0 else None)
            bottoms += vals
    # totals as text
    for i, r in enumerate(results):
        cpp_total = sum(r["cpp_phases"])
        py_total = sum(r["py_phases"])
        ax.text(i - width / 2, cpp_total, f"{cpp_total:.1f}s", ha="center", va="bottom", fontsize=8)
        ax.text(i + width / 2, py_total, f"{py_total:.1f}s", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([r["case"]["label"] for r in results], fontsize=8)
    ax.set_ylabel("wall-clock time (s)")
    ax.set_title("Runtime by phase (solid = C++, hatched = Python)\n"
                 "'model construction' is dominated by C++'s one-time FFTW_EXHAUSTIVE planning --\n"
                 "compare 'timestepping' for actual per-step computational cost")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "runtime_phase_breakdown.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    if not CPP_BIN.exists():
        print(
            f"ERROR: {CPP_BIN} not found.\n"
            "Build the C++ reference first (see the top-level README:\n"
            "  cd build && make\n"
            ")\nThis script deliberately does not fabricate comparison data "
            "when the reference binary is unavailable.",
            file=sys.stderr,
        )
        sys.exit(1)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("Excluded cases (not run; see script docstring for why):")
    for c in EXCLUDED_CASES:
        print(f"  - {c['name']} ({c['geom']}): {c['reason']}")
    print()

    results = [run_case(case) for case in CASES]

    _fig_stitched_force(results)
    _fig_stitched_vorticity(results)
    _fig_runtime_phases(results)

    metrics_path = RESULTS_DIR / "validation_metrics.csv"
    with open(metrics_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "quantity", "error_type", "max_error", "mean_error",
                    "max_relative_to_peak", "mean_relative_to_peak"])
        for r in results:
            name = r["case"]["name"]
            w.writerow([name, "Cd (drag coefficient)", "relative (|py-cpp|/|cpp|, t>0)",
                        f"{r['cd_rel_err'].max():.3e}", f"{r['cd_rel_err'].mean():.3e}", "", ""])
            w.writerow([name, "Cl (lift coefficient)", "absolute (reference ~0 by symmetry)",
                        f"{r['cl_abs_err'].max():.3e}", f"{r['cl_abs_err'].mean():.3e}", "", ""])
            dt, nsteps = r["dt"], r["nsteps"]
            for step in SNAPSHOT_STEPS:
                amax, amean, rmax, rmean = field_errors(r["omega_cpp"][step], r["omega_py"][step])
                w.writerow([name, f"vorticity field, t={step * dt:g}", "absolute",
                            f"{amax:.3e}", f"{amean:.3e}", f"{rmax:.3e}", f"{rmean:.3e}"])
            model_c, setup_c, steps_c = r["cpp_phases"]
            model_p, setup_p, steps_p = r["py_phases"]
            w.writerow([name, "runtime, C++: model construction", "wall-clock seconds", f"{model_c:.3f}", "", "", ""])
            w.writerow([name, "runtime, C++: solver factorization", "wall-clock seconds", f"{setup_c:.3f}", "", "", ""])
            w.writerow([name, f"runtime, C++: timestepping ({nsteps} steps)", "wall-clock seconds", f"{steps_c:.3f}", "", "", ""])
            w.writerow([name, "runtime, Python: model construction", "wall-clock seconds", f"{model_p:.3f}", "", "", ""])
            w.writerow([name, "runtime, Python: solver factorization", "wall-clock seconds", f"{setup_p:.3f}", "", "", ""])
            w.writerow([name, f"runtime, Python: timestepping ({nsteps} steps)", "wall-clock seconds", f"{steps_p:.3f}", "", "", ""])
            w.writerow([name, "runtime ratio, Python/C++: timestepping only", "dimensionless",
                        f"{steps_p / steps_c:.3f}", "", "", ""])
            w.writerow([name, "ms/step, C++: timestepping", "milliseconds", f"{1000*steps_c/nsteps:.3f}", "", "", ""])
            w.writerow([name, "ms/step, Python: timestepping", "milliseconds", f"{1000*steps_p/nsteps:.3f}", "", "", ""])
        for c in EXCLUDED_CASES:
            w.writerow([c["name"], "EXCLUDED", c["reason"], "", "", "", ""])

    print(f"\nWrote per-case figures under {FIGURES_DIR}/<case>/")
    print(f"Wrote stitched figures to {FIGURES_DIR}")
    print(f"Wrote metrics to {metrics_path}\n")
    for r in results:
        name = r["case"]["name"]
        nsteps = r["nsteps"]
        model_c, setup_c, steps_c = r["cpp_phases"]
        model_p, setup_p, steps_p = r["py_phases"]
        print(f"[{name}] C++    total={sum(r['cpp_phases']):6.2f}s  "
              f"(model={model_c:5.2f}s  factorize={setup_c:5.2f}s  steps={steps_c:5.2f}s -> {1000*steps_c/nsteps:.2f} ms/step)")
        print(f"[{name}] Python total={sum(r['py_phases']):6.2f}s  "
              f"(model={model_p:5.2f}s  factorize={setup_p:5.2f}s  steps={steps_p:5.2f}s -> {1000*steps_p/nsteps:.2f} ms/step)")
        print(f"[{name}] Cd max relative error: {r['cd_rel_err'].max():.3e}   "
              f"Cl max absolute error: {r['cl_abs_err'].max():.3e}")


if __name__ == "__main__":
    main()
