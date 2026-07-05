# gen_cost_report.py
#
# Reads results/cost/raw/cost_results.json (produced by run_benchmark.py)
# and produces the tables (CSV + Markdown) and figures (matplotlib PNGs)
# for the computational cost analysis. Nothing here is hand-drawn or
# AI-generated -- every number traces back to an actual `build/ibpm` /
# `python3 results/cost/run_ibpm_case.py` run, measured with
# resource.getrusage (CPU time, peak RSS) and psutil (RSS/CPU% time
# series), see run_benchmark.py.

from __future__ import annotations

import csv
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COST_DIR = pathlib.Path(__file__).resolve().parent
RAW_DIR = COST_DIR / "raw"
TABLES_DIR = COST_DIR / "tables"
TABLES_DIR.mkdir(exist_ok=True)

results = json.loads((RAW_DIR / "cost_results.json").read_text())

BACKEND_STYLE = {
    "cpp": dict(color="C1", marker="o", label="C++ (build/ibpm)"),
    "python": dict(color="C0", marker="s", label="Python (py/ibpm.py)"),
    "jax": dict(color="C2", marker="^", label="JAX (planned, not yet run)"),
}


def by_backend(backend):
    rows = [r for r in results if r["backend"] == backend]
    return sorted(rows, key=lambda r: r["nx"])


BACKENDS_PRESENT = sorted({r["backend"] for r in results})

# ======================================================================
# Table 1: full per-run summary (CSV + Markdown)
# ======================================================================
csv_path = TABLES_DIR / "cost_summary.csv"
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["backend", "nx", "ny", "dx", "dt", "nsteps", "npts",
                "wall_time_s", "phase_model_s", "phase_setup_s", "phase_steps_s",
                "ms_per_step", "cpu_user_s", "cpu_sys_s", "cpu_total_s",
                "peak_rss_MB", "cpu_efficiency_pct"])
    for r in sorted(results, key=lambda r: (r["backend"], r["nx"])):
        ms_per_step = 1000 * r["phase_steps"] / r["nsteps"]
        cpu_eff = 100 * r["cpu_total"] / r["wall_time"] if r["wall_time"] > 0 else float("nan")
        w.writerow([r["backend"], r["nx"], r["ny"], f"{r['dx']:.5f}", r["dt"], r["nsteps"], r["npts"],
                    f"{r['wall_time']:.3f}", f"{r['phase_model']:.3f}", f"{r['phase_setup']:.3f}",
                    f"{r['phase_steps']:.3f}", f"{ms_per_step:.3f}",
                    f"{r['cpu_user']:.3f}", f"{r['cpu_sys']:.3f}", f"{r['cpu_total']:.3f}",
                    f"{r['peak_rss_bytes']/1e6:.2f}", f"{cpu_eff:.1f}"])
print(f"wrote {csv_path}")

md_path = TABLES_DIR / "cost_summary.md"
with open(md_path, "w") as f:
    f.write("# Computational cost summary\n\n")
    nsteps_label = results[0]["nsteps"] if results else "N"
    f.write(f"Cylinder, Re=100, RK3, ngrid=1, domain [-2,2]x[-2,2], {nsteps_label} steps "
            "(dt=0.01 at nx=400, else 0.02 -- see run_benchmark.py for why).\n\n")
    f.write("| backend | nx=ny | dx | wall (s) | model (s) | setup (s) | "
            "steps (s) | ms/step | CPU (s) | peak RSS (MB) | CPU eff. (%) |\n")
    f.write("|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in sorted(results, key=lambda r: (r["backend"], r["nx"])):
        ms_per_step = 1000 * r["phase_steps"] / r["nsteps"]
        cpu_eff = 100 * r["cpu_total"] / r["wall_time"] if r["wall_time"] > 0 else float("nan")
        f.write(f"| {r['backend']} | {r['nx']} | {r['dx']:.4f} | {r['wall_time']:.2f} | "
                f"{r['phase_model']:.2f} | {r['phase_setup']:.2f} | {r['phase_steps']:.2f} | "
                f"{ms_per_step:.2f} | {r['cpu_total']:.2f} | {r['peak_rss_bytes']/1e6:.1f} | "
                f"{cpu_eff:.0f} |\n")
    f.write("\n_JAX row intentionally absent: not implemented yet, see backends.py._\n")
print(f"wrote {md_path}")

# Table 2: head-to-head ratio table (Python/C++), the single most useful summary
ratio_path = TABLES_DIR / "cost_ratio_python_vs_cpp.md"
with open(ratio_path, "w") as f:
    f.write("# Python / C++ cost ratio (>1 = Python costs more)\n\n")
    f.write("| nx=ny | wall time ratio | timestepping-only ratio | peak RSS ratio | CPU-seconds ratio |\n")
    f.write("|---|---|---|---|---|\n")
    cpp_rows = {r["nx"]: r for r in by_backend("cpp")}
    py_rows = {r["nx"]: r for r in by_backend("python")}
    for nx in sorted(set(cpp_rows) & set(py_rows)):
        c, p = cpp_rows[nx], py_rows[nx]
        f.write(f"| {nx} | {p['wall_time']/c['wall_time']:.2f}x | "
                f"{p['phase_steps']/c['phase_steps']:.2f}x | "
                f"{p['peak_rss_bytes']/c['peak_rss_bytes']:.2f}x | "
                f"{p['cpu_total']/c['cpu_total']:.2f}x |\n")
print(f"wrote {ratio_path}")

# ======================================================================
# Figure 1: wall-clock time vs grid size (total, log-log)
# ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for backend in BACKENDS_PRESENT:
    rows = by_backend(backend)
    nx = [r["nx"] for r in rows]
    wall = [r["wall_time"] for r in rows]
    steps = [r["phase_steps"] for r in rows]
    style = BACKEND_STYLE[backend]
    axes[0].plot(nx, wall, "-", marker=style["marker"], color=style["color"], label=style["label"])
    axes[1].plot(nx, steps, "-", marker=style["marker"], color=style["color"], label=style["label"])
for ax, title in zip(axes, ["Total wall-clock time (includes one-time setup)",
                             "Timestepping-only wall-clock time"]):
    ax.set_xlabel("grid resolution (nx = ny)")
    ax.set_ylabel("wall-clock time (s)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
fig.suptitle("Runtime vs. grid resolution: total (misleading -- see right panel) vs. "
             "timestepping-only cost")
fig.tight_layout()
fig.savefig(COST_DIR / "runtime_vs_gridsize.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote runtime_vs_gridsize.png")

# ======================================================================
# Figure 2: phase breakdown (stacked bars), same convention as
# results/generate_validation_report.py's runtime_phase_breakdown.png
# ======================================================================
fig, ax = plt.subplots(figsize=(10, 5.5))
labels_phase = ["model\nconstruction", "solver\nfactorization", f"timestepping\n({results[0]['nsteps']} steps)"]
colors = ["0.65", "#8ecae6", "#023047"]
nxs = sorted({r["nx"] for r in results})
width = 0.35
x = np.arange(len(nxs))
for backend, offset, hatch in [("cpp", -width / 2, None), ("python", width / 2, "//")]:
    rows_by_nx = {r["nx"]: r for r in by_backend(backend)}
    bottoms = np.zeros(len(nxs))
    for phase_idx, (key, color) in enumerate(zip(["phase_model", "phase_setup", "phase_steps"], colors)):
        vals = np.array([rows_by_nx[nx][key] if nx in rows_by_nx else 0.0 for nx in nxs])
        ax.bar(x + offset, vals, width, bottom=bottoms, color=color, hatch=hatch,
                edgecolor="black", linewidth=0.5,
                label=(labels_phase[phase_idx] if backend == "cpp" else None))
        bottoms += vals
ax.set_xticks(x)
ax.set_xticklabels([str(n) for n in nxs])
ax.set_xlabel("grid resolution (nx = ny)")
ax.set_ylabel("wall-clock time (s)")
ax.set_title("Runtime by phase (solid = C++, hatched = Python)\n"
             "'model construction' is dominated by FFTW planning in C++ (one-time cost); "
             "compare 'timestepping' for actual per-step cost")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(COST_DIR / "runtime_phase_breakdown.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote runtime_phase_breakdown.png")

# ======================================================================
# Figure 3: peak RAM vs grid size
# ======================================================================
fig, ax = plt.subplots(figsize=(7, 5))
for backend in BACKENDS_PRESENT:
    rows = by_backend(backend)
    nx = [r["nx"] for r in rows]
    rss_mb = [r["peak_rss_bytes"] / 1e6 for r in rows]
    style = BACKEND_STYLE[backend]
    ax.plot(nx, rss_mb, "-", marker=style["marker"], color=style["color"], label=style["label"])
ax.set_xlabel("grid resolution (nx = ny)")
ax.set_ylabel("peak resident set size (MB)")
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_title("Peak RAM usage vs. grid resolution")
ax.legend(fontsize=9)
ax.grid(alpha=0.3, which="both")
fig.savefig(COST_DIR / "peak_ram_vs_gridsize.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote peak_ram_vs_gridsize.png")

# ======================================================================
# Figure 4: CPU-seconds vs grid size (user+sys), plus CPU efficiency (%)
# ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for backend in BACKENDS_PRESENT:
    rows = by_backend(backend)
    nx = [r["nx"] for r in rows]
    cpu = [r["cpu_total"] for r in rows]
    eff = [100 * r["cpu_total"] / r["wall_time"] for r in rows]
    style = BACKEND_STYLE[backend]
    axes[0].plot(nx, cpu, "-", marker=style["marker"], color=style["color"], label=style["label"])
    axes[1].plot(nx, eff, "-", marker=style["marker"], color=style["color"], label=style["label"])
axes[0].set_ylabel("CPU time, user+sys (s)")
axes[0].set_yscale("log")
axes[0].set_title("Total CPU-seconds consumed")
axes[1].set_ylabel("CPU efficiency = 100 x CPU-seconds / wall-seconds (%)")
axes[1].axhline(100, color="k", ls=":", lw=1)
axes[1].set_title("CPU efficiency (100% = fully single-core-bound,\n>100% = multi-threaded, <100% = I/O/blocked)")
for ax in axes:
    ax.set_xlabel("grid resolution (nx = ny)")
    ax.set_xscale("log", base=2)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(COST_DIR / "cpu_usage_vs_gridsize.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote cpu_usage_vs_gridsize.png")

# ======================================================================
# Figure 5: time-series RSS and CPU% for the largest grid size (richest signal)
# ======================================================================
nx_focus = max(nxs)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for backend in BACKENDS_PRESENT:
    rows = [r for r in by_backend(backend) if r["nx"] == nx_focus]
    if not rows:
        continue
    r = rows[0]
    ts = r["timeseries"]
    if not ts:
        continue
    t = [s["t"] for s in ts]
    rss = [s["rss_bytes"] / 1e6 for s in ts]
    cpu = [s["cpu_percent"] for s in ts]
    style = BACKEND_STYLE[backend]
    axes[0].plot(t, rss, color=style["color"], label=style["label"])
    axes[1].plot(t, cpu, color=style["color"], label=style["label"])
axes[0].set_ylabel("RSS memory (MB)")
axes[0].set_title(f"Memory over time, nx=ny={nx_focus}")
axes[1].set_ylabel("CPU usage (%, 100% = 1 core saturated)")
axes[1].set_title(f"CPU usage over time, nx=ny={nx_focus}")
for ax in axes:
    ax.set_xlabel("wall-clock time since process start (s)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(COST_DIR / "timeseries_largest_case.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote timeseries_largest_case.png")

print("done")
