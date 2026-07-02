#!/usr/bin/env python3
"""
generate_validation_report.py

Validates the Python port (py/) against the reference C++ implementation
(src/, compiled as build/ibpm) by running BOTH on the standard
examples/cylinder.geom test case (the default configuration used by
examples/ibpm.cmd: 200x200 grid, Re=100, dt=0.02, 250 steps) and comparing:

  - Drag/lift coefficients (Cd, Cl) over the full 250-step run
  - Absolute force error vs. time (a dense, per-timestep error signal)
  - Vorticity field snapshots at t=0, 2, 4 (steps 0, 100, 200)
  - Wall-clock runtime (cold start, i.e. including one-time Cholesky
    factorization of the projection operator)

Every number and figure this script produces comes directly from freshly
re-running the two actual simulators (subprocess calls to build/ibpm and
`python3 -m py.ibpm`) and reading their binary/text output -- nothing here
is hand-entered or AI-generated.

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
    figures/force_coefficients_vs_time.png
    figures/force_error_vs_time.png
    figures/vorticity_field_comparison.png
    figures/vorticity_parity_plot.png
    figures/runtime_comparison.png
    figures/flow_evolution_python.png
    validation_metrics.csv
"""

from __future__ import annotations

import csv
import subprocess
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
RUN_DATA_DIR = RESULTS_DIR / "_run_data"  # raw run output, gitignored

CPP_BIN = REPO_ROOT / "build" / "ibpm"
GEOM_FILE = REPO_ROOT / "examples" / "cylinder.geom"

# Configuration matching examples/ibpm.cmd's defaults, made explicit here
# so this script is self-documenting and unaffected if the C++/Python
# defaults ever change.
RUN_ARGS = [
    "-nx", "200", "-ny", "200", "-ngrid", "1",
    "-length", "4.0", "-xoffset", "-2.0", "-yoffset", "-2.0",
    "-Re", "100", "-dt", "0.02", "-nsteps", "250", "-scheme", "rk3",
    "-tecplot", "100", "-restart", "100", "-force", "1", "-energy", "0",
]
NX, NY, LENGTH, XOFFSET, YOFFSET = 200, 200, 4.0, -2.0, -2.0
DX = LENGTH / NX
SNAPSHOT_STEPS = [0, 100, 200]  # matches -restart 100 above


def run_cpp(outdir: Path) -> float:
    """Run the compiled C++ reference binary; return wall-clock seconds."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(CPP_BIN), "-geom", str(GEOM_FILE), "-outdir", str(outdir), "-name", "ibpm", *RUN_ARGS]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    (outdir / "stdout.log").write_text(proc.stdout)
    (outdir / "stderr.log").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"C++ run failed (see {outdir}/stderr.log)")
    return elapsed


def run_python(outdir: Path) -> float:
    """Run the Python port as a subprocess (fair wall-clock comparison,
    including interpreter startup); return wall-clock seconds."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "py.ibpm", "-geom", str(GEOM_FILE), "-outdir", str(outdir), "-name", "ibpm", *RUN_ARGS]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    (outdir / "stdout.log").write_text(proc.stdout)
    (outdir / "stderr.log").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"Python run failed (see {outdir}/stderr.log)")
    return elapsed


def load_omega(bin_path: Path) -> np.ndarray:
    """Load the finest-level vorticity field from an IBPM restart file,
    using the already-validated py.state.State reader (same binary
    format for both the C++ and Python outputs)."""
    sys.path.insert(0, str(REPO_ROOT))
    from py.state import State  # local import: keeps `py` off sys.path for --help etc.

    s = State(filename=str(bin_path))
    return s.omega._data[0].copy()  # shape (nx-1, ny-1), finest grid level


def load_force(force_path: Path) -> np.ndarray:
    """Load a .force file: columns [timestep, time, Cd, Cl]."""
    return np.loadtxt(force_path)


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
    cpp_dir = RUN_DATA_DIR / "cpp"
    py_dir = RUN_DATA_DIR / "python"

    print("Running C++ reference (build/ibpm)...")
    cpp_time = run_cpp(cpp_dir)
    print(f"  done in {cpp_time:.2f} s")

    print("Running Python port (python3 -m py.ibpm)...")
    py_time = run_python(py_dir)
    print(f"  done in {py_time:.2f} s")

    # ------------------------------------------------------------------
    # Force (Cd, Cl) comparison
    # ------------------------------------------------------------------
    cpp_force = load_force(cpp_dir / "ibpm.force")
    py_force = load_force(py_dir / "ibpm.force")
    assert cpp_force.shape == py_force.shape, "step count mismatch between runs"

    t = cpp_force[:, 1]
    cd_cpp, cl_cpp = cpp_force[:, 2], cpp_force[:, 3]
    cd_py, cl_py = py_force[:, 2], py_force[:, 3]

    cd_abs_err = np.abs(cd_py - cd_cpp)
    cl_abs_err = np.abs(cl_py - cl_cpp)
    # Cd is O(1) throughout (after the initial impulsive-start transient at
    # t=0, where Cd is exactly 0 in both and relative error is undefined);
    # report relative error over t>0.
    nonzero = cd_cpp != 0
    cd_rel_err = cd_abs_err[nonzero] / np.abs(cd_cpp[nonzero])

    # --- Figure 1: Cd/Cl vs time, Python overlaid on C++ ---
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    axes[0].plot(t, cd_cpp, "-", color="tab:blue", lw=2.5, alpha=0.5, label="C++ (src/)")
    axes[0].plot(t, cd_py, "--", color="tab:orange", lw=1.5, label="Python (py/)")
    axes[0].set_ylabel(r"$C_d$ (drag coefficient)")
    axes[0].legend()
    axes[0].set_title(f"Cylinder, Re=100, {NX}x{NY} grid: force coefficients vs. time")

    axes[1].plot(t, cl_cpp, "-", color="tab:blue", lw=2.5, alpha=0.5, label="C++ (src/)")
    axes[1].plot(t, cl_py, "--", color="tab:orange", lw=1.5, label="Python (py/)")
    axes[1].set_ylabel(r"$C_l$ (lift coefficient)")
    axes[1].set_xlabel("time")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "force_coefficients_vs_time.png", dpi=150)
    plt.close(fig)

    # --- Figure 2: force error vs time (dense, per-timestep signal) ---
    # NOTE: a floor of 1e-18 is applied only so exact-zero errors are
    # visible on a log-scale axis; it is not a measured error value.
    _FLOOR = 1e-18
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    axes[0].semilogy(t, np.maximum(cd_abs_err, _FLOOR), color="tab:red")
    axes[0].set_ylabel(r"$|C_{d,\mathrm{py}} - C_{d,\mathrm{cpp}}|$")
    axes[0].set_title("Python-vs-C++ absolute force error vs. time")
    if cd_abs_err.max() <= _FLOOR:
        axes[0].text(0.5, 0.5, "error is exactly 0.0 at every timestep\n(line at floor is a plotting artifact)",
                     transform=axes[0].transAxes, ha="center", va="center", fontsize=9, color="dimgray")
    axes[1].semilogy(t, np.maximum(cl_abs_err, _FLOOR), color="tab:red")
    axes[1].set_ylabel(r"$|C_{l,\mathrm{py}} - C_{l,\mathrm{cpp}}|$")
    axes[1].set_xlabel("time")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "force_error_vs_time.png", dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Vorticity field comparison at 3 snapshots
    # ------------------------------------------------------------------
    xs = XOFFSET + np.arange(1, NX) * DX  # interior-node coordinates (see Grid.getXEdge)
    ys = YOFFSET + np.arange(1, NY) * DX
    X, Y = np.meshgrid(xs, ys, indexing="ij")

    omega_cpp = {}
    omega_py = {}
    for step in SNAPSHOT_STEPS:
        tag = f"{step:05d}"
        omega_cpp[step] = load_omega(cpp_dir / f"ibpm{tag}.bin")
        omega_py[step] = load_omega(py_dir / f"ibpm{tag}.bin")

    vmax = max(np.abs(omega_cpp[s]).max() for s in SNAPSHOT_STEPS)
    fig, axes = plt.subplots(len(SNAPSHOT_STEPS), 3, figsize=(12, 3.3 * len(SNAPSHOT_STEPS)))
    for row, step in enumerate(SNAPSHOT_STEPS):
        time_val = step * 0.02
        diff = omega_py[step] - omega_cpp[step]
        panels = [
            (omega_cpp[step], f"C++, t={time_val:g}", "RdBu_r", -vmax, vmax),
            (omega_py[step], f"Python, t={time_val:g}", "RdBu_r", -vmax, vmax),
            (diff, f"Python - C++ (max |diff| = {np.abs(diff).max():.1e})", "PuOr", None, None),
        ]
        for col, (field, title, cmap, vmin_, vmax_) in enumerate(panels):
            ax = axes[row, col]
            im = ax.contourf(X, Y, field, levels=41, cmap=cmap, vmin=vmin_, vmax=vmax_)
            circle = plt.Circle((0, 0), 0.5, fill=False, color="k", lw=1.2)
            ax.add_patch(circle)
            ax.set_aspect("equal")
            ax.set_xlim(-2, 3)
            ax.set_ylim(-2, 2)
            ax.set_title(title, fontsize=10)
            fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("Vorticity field: C++ vs. Python vs. difference, at 3 snapshots")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "vorticity_field_comparison.png", dpi=150)
    plt.close(fig)

    # --- Figure: parity plot (Python vs C++ vorticity, all interior points, final snapshot) ---
    final_step = SNAPSHOT_STEPS[-1]
    a = omega_cpp[final_step].ravel()
    b = omega_py[final_step].ravel()
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(a, b, s=2, alpha=0.3, color="tab:blue")
    lims = [min(a.min(), b.min()), max(a.max(), b.max())]
    ax.plot(lims, lims, "k--", lw=1, label="y = x")
    ax.set_xlabel(r"C++ vorticity $\omega$")
    ax.set_ylabel(r"Python vorticity $\omega$")
    ax.set_title(f"Parity plot: all interior points, t={final_step * 0.02:g}\n"
                 f"({(NX - 1) * (NY - 1):,} points, $R^2$ = {np.corrcoef(a, b)[0, 1]**2:.10f})")
    ax.legend()
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "vorticity_parity_plot.png", dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Runtime comparison
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5, 4))
    bars = ax.bar(["C++\n(build/ibpm)", "Python\n(py/ibpm.py)"], [cpp_time, py_time],
                   color=["tab:blue", "tab:orange"])
    for rect, val in zip(bars, [cpp_time, py_time]):
        ax.text(rect.get_x() + rect.get_width() / 2, val, f"{val:.2f} s",
                ha="center", va="bottom")
    ax.set_ylabel("wall-clock time (s)")
    ax.set_title(f"Total runtime, {NX}x{NY} grid, 250 steps\n(cold start, incl. Cholesky factorization)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "runtime_comparison.png", dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Bonus: flow evolution (Python only) -- physical visualization
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, len(SNAPSHOT_STEPS), figsize=(4.2 * len(SNAPSHOT_STEPS), 4.5),
                              constrained_layout=True)
    for ax, step in zip(axes, SNAPSHOT_STEPS):
        im = ax.contourf(X, Y, omega_py[step], levels=41, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        circle = plt.Circle((0, 0), 0.5, fill=False, color="k", lw=1.2)
        ax.add_patch(circle)
        ax.set_aspect("equal")
        ax.set_xlim(-2, 3)
        ax.set_ylim(-2, 2)
        ax.set_title(f"t = {step * 0.02:g}")
    fig.suptitle("Flow development behind an impulsively-started cylinder (Python port, Re=100)")
    fig.colorbar(im, ax=axes, shrink=0.7, label="vorticity")
    fig.savefig(FIGURES_DIR / "flow_evolution_python.png", dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Summary metrics table
    # ------------------------------------------------------------------
    def field_errors(step: int) -> "tuple[float, float, float, float]":
        cpp_f = omega_cpp[step]
        diff = np.abs(omega_py[step] - cpp_f)
        peak = np.abs(cpp_f).max()
        if peak == 0:
            # t=0: zero initial condition, field is exactly zero in both --
            # 0/0 is undefined, but "no error" is the correct reading here
            # (diff.max() is itself 0.0 in this case).
            return diff.max(), diff.mean(), 0.0, 0.0
        return diff.max(), diff.mean(), diff.max() / peak, diff.mean() / peak

    metrics_path = RESULTS_DIR / "validation_metrics.csv"
    with open(metrics_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "error_type", "max_error", "mean_error",
                    "max_relative_to_peak", "mean_relative_to_peak"])
        w.writerow(["Cd (drag coefficient)", "relative (|py-cpp|/|cpp|, t>0)",
                    f"{cd_rel_err.max():.3e}", f"{cd_rel_err.mean():.3e}", "", ""])
        w.writerow(["Cl (lift coefficient)", "absolute (reference ~0 by symmetry)",
                    f"{cl_abs_err.max():.3e}", f"{cl_abs_err.mean():.3e}", "", ""])
        for step in SNAPSHOT_STEPS:
            amax, amean, rmax, rmean = field_errors(step)
            w.writerow([f"vorticity field, t={step * 0.02:g}", "absolute",
                        f"{amax:.3e}", f"{amean:.3e}", f"{rmax:.3e}", f"{rmean:.3e}"])
        w.writerow(["runtime, C++ (build/ibpm)", "wall-clock seconds", f"{cpp_time:.3f}", "", "", ""])
        w.writerow(["runtime, Python (py/ibpm.py)", "wall-clock seconds", f"{py_time:.3f}", "", "", ""])
        w.writerow(["runtime, Python / C++ ratio", "dimensionless", f"{py_time / cpp_time:.3f}", "", "", ""])

    print(f"\nWrote figures to {FIGURES_DIR}")
    print(f"Wrote metrics to {metrics_path}")
    print(f"\nC++ time:    {cpp_time:.2f} s")
    print(f"Python time: {py_time:.2f} s")
    for step in SNAPSHOT_STEPS:
        amax, amean, rmax, rmean = field_errors(step)
        print(f"vorticity t={step*0.02:g}: max|diff|={amax:.3e}  (rel. to peak: {rmax:.3e})")
    print(f"Cd max relative error: {cd_rel_err.max():.3e}")
    print(f"Cl max absolute error: {cl_abs_err.max():.3e}")


if __name__ == "__main__":
    main()
