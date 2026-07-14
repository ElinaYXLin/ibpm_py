"""
gen_dt_refinement_figs.py

Figure for run_dt_refinement.py's decisive test: does refining dt (at
fixed dx=0.04) delay or remove SD8000's coarse-grid blow-up? 5 repeats
per dt isolate genuine dt-dependence from the run-to-run FFTW_EXHAUSTIVE
replanning noise documented in ../7-chaos_sensitivity/README.md.

Usage: python3 SURF_test/airfoils/8-dt_refinement_and_spectra/gen_dt_refinement_figs.py
Output: SURF_test/airfoils/8-dt_refinement_and_spectra/figures/dt_refinement_blowup.png
        SURF_test/airfoils/8-dt_refinement_and_spectra/data/dt_refinement_summary.txt
"""
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
HERE = REPO / "SURF_test" / "airfoils" / "8-dt_refinement_and_spectra"
RUN_DATA = HERE / "_run_data"
FIGS = HERE / "figures"; FIGS.mkdir(parents=True, exist_ok=True)
DATA = HERE / "data"; DATA.mkdir(parents=True, exist_ok=True)

DTS = [0.01, 0.005, 0.0025]
N_REPEATS = 5
THRESHOLD = 20.0
T_FINAL = 40.0


def find_blowup(force_path, dt, threshold=THRESHOLD):
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    step, cd, cl = d[:, 0], d[:, 2], d[:, 3]
    mask = (np.abs(cl) > threshold) | (np.abs(cd) > threshold)
    if not mask.any():
        return None
    return int(step[np.argmax(mask)]) * dt  # physical time


def main():
    results = {}  # dt -> list of blowup times (None if never)
    for dt in DTS:
        times = []
        for rep in range(N_REPEATS):
            fp = RUN_DATA / f"dt{dt}_rep{rep:02d}" / "run.force"
            times.append(find_blowup(fp, dt))
        results[dt] = times

    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = {0.01: "C0", 0.005: "C1", 0.0025: "C2"}
    means, stds = [], []
    for i, dt in enumerate(DTS):
        times = results[dt]
        finite = [t for t in times if t is not None]
        n_never = len(times) - len(finite)
        x = np.full(len(times), i) + np.random.uniform(-0.06, 0.06, len(times))
        y = [t if t is not None else T_FINAL for t in times]
        markers = ["x" if t is None else "o" for t in times]
        for xi, yi, m in zip(x, y, markers):
            ax.scatter(xi, yi, marker=m, s=70, color=colors[dt],
                       edgecolor="k", linewidth=0.6, zorder=3)
        if finite:
            mean_t, std_t = np.mean(finite), np.std(finite)
            means.append(mean_t); stds.append(std_t)
            ax.errorbar(i, mean_t, yerr=std_t, fmt="_", color="k", ms=20,
                        capsize=6, elinewidth=1.5, zorder=4)
        else:
            means.append(np.nan); stds.append(np.nan)
        if n_never:
            ax.annotate(f"{n_never}/{len(times)} never blew up", (i, -2.5), ha="center", fontsize=8,
                        annotation_clip=False, color="C2")

    ax.axhline(T_FINAL, color="gray", ls=":", lw=1, label=f"run length cap (t={T_FINAL})")
    ax.set_xticks(range(len(DTS))); ax.set_xticklabels([f"dt={dt}" for dt in DTS])
    ax.set_ylabel("physical time $t$ of blow-up ($|C_l|$ or $|C_d|$ > 20)")
    ax.set_xlabel("timestep (grid dx=0.04 fixed)")
    ax.set_title("SD8000 coarse grid: does refining dt delay the blow-up?\n"
                 "(o = blew up, x = ran clean to t=40; 5 repeats/dt isolate FFTW-replanning noise; "
                 "black bar = mean±std)")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_ylim(bottom=-5)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    outp = FIGS / "dt_refinement_blowup.png"
    fig.savefig(outp, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outp}")

    with open(DATA / "dt_refinement_summary.txt", "w") as f:
        f.write("SD8000 coarse grid (dx=0.04), C++, 5 repeats per dt\n")
        f.write("blow-up threshold: |Cl| or |Cd| > 20\n\n")
        for dt in DTS:
            times = results[dt]
            finite = [t for t in times if t is not None]
            f.write(f"dt={dt}: {[f'{t:.3f}' if t is not None else 'none' for t in times]}\n")
            if finite:
                f.write(f"  mean={np.mean(finite):.3f}  std={np.std(finite):.3f}  "
                        f"n_blowup={len(finite)}/{len(times)}\n")
            else:
                f.write(f"  no blow-ups in any of {len(times)} repeats\n")
    print(f"wrote {DATA / 'dt_refinement_summary.txt'}")
    print(open(DATA / "dt_refinement_summary.txt").read())


if __name__ == "__main__":
    main()
