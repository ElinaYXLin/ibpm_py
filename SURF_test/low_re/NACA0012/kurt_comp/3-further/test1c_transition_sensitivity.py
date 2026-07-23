"""
test1c_transition_sensitivity.py

Follow-up to a gap identified in Test 1b's justification: Test 1b claims
the 34-35deg thrust/lift/Strouhal jump is "explained by Question 2's
grid-resolution-driven mode transition," but Question 2's own grid/domain
sensitivity checks (tests 2b/2c) only ever probed 4 representative angles
INSIDE the plateaus (15, 20, 30, 40) -- none of them AT the 34-35 transition
itself. So the claim that the transition is grid/domain-sensitive was an
extrapolation from nearby points, never directly tested. This script tests
it directly: does the transition's LOCATION (which pair of adjacent
integer angles shows the jump) or its SIZE move under the same two knobs
(dx refinement, ngrid/domain refinement) already used for tests 2b/2c?

Steady, py_static only, Re=1000, same conventions as run_further.py.
New angles: 33,34,35,36 (straddling the transition found in Test 1b) at
dx=0.01 and at ngrid=2,3 (dx=0.02). dx=0.02/ngrid=1 baseline already exists
in ../1-paper_based/runs/dx0.020/steady_py_a{33..36}.

Usage:
  python3 test1c_transition_sensitivity.py run       # launch new runs
  python3 test1c_transition_sensitivity.py analyze   # analyze + plot
Output: runs/transition_check/, figures/test_1c_transition_sensitivity.png,
        data/test1c_transition_sensitivity.csv
"""
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_further as rf  # reuses GEOM, run_case, REPO, RUNS, PAPER_RUNS

TRANS_ANGLES = [33, 34, 35, 36]
FIGS = rf.FURTHER / "figures"
DATA = rf.FURTHER / "data"


def build_jobs():
    jobs = []
    for a in TRANS_ANGLES:
        outdir = rf.RUNS / "transition_check" / "dx_refine" / f"dx0.010_a{a:02d}"
        jobs.append(dict(outdir=outdir, geom=rf.GEOM["0.01"], nx=600, ny=300, ngrid=1,
                          domain=rf.BASE_DOMAIN, alpha=a, Re=rf.RE, dt=0.005, nsteps=6000))
    ngrid_domain = dict(length=6, xoffset=-2, yoffset=-1.52)
    for ngrid in (2, 3):
        for a in TRANS_ANGLES:
            outdir = rf.RUNS / "transition_check" / "ngrid_sweep" / f"ngrid{ngrid}_a{a:02d}"
            jobs.append(dict(outdir=outdir, geom=rf.GEOM["0.02"], nx=300, ny=152, ngrid=ngrid,
                              domain=ngrid_domain, alpha=a, Re=rf.RE, dt=0.01, nsteps=3000))
    return jobs


def do_run():
    jobs = build_jobs()
    print(f"{len(jobs)} jobs", flush=True)
    for j in jobs:
        name, status, elapsed = rf.run_case(**j)
        print(f"  {name}: {status} ({elapsed:.0f}s)", flush=True)


def _cl_st(run_dir):
    d = np.loadtxt(run_dir / "flow.force") if (run_dir / "flow.force").exists() else None
    if d is None:
        return None, None
    n = len(d)
    seg = d[int(n * 0.5):]
    t, cd, cl = seg[:, 1], seg[:, 2], seg[:, 3]
    cl_mean = float(np.mean(cl))
    cl0 = cl - cl_mean
    dt_ = np.mean(np.diff(t))
    win = np.hanning(len(cl0))
    amp = np.abs(np.fft.rfft(cl0 * win))
    freqs = np.fft.rfftfreq(len(cl0), d=dt_)
    band = (freqs > 0.1) & (freqs < 3.0)
    st = float(freqs[band][np.argmax(amp[band])]) if band.any() and amp[band].max() > 5 * amp[band].mean() else 0.0
    return cl_mean, st


def analyze():
    rows = []
    configs = [
        ("dx=0.02 (baseline)", lambda a: rf.PAPER_RUNS / f"steady_py_a{a:02d}"),
        ("dx=0.01", lambda a: rf.RUNS / "transition_check" / "dx_refine" / f"dx0.010_a{a:02d}"),
        ("ngrid=2", lambda a: rf.RUNS / "transition_check" / "ngrid_sweep" / f"ngrid2_a{a:02d}"),
        ("ngrid=3", lambda a: rf.RUNS / "transition_check" / "ngrid_sweep" / f"ngrid3_a{a:02d}"),
    ]
    for label, path_fn in configs:
        for a in TRANS_ANGLES:
            rd = path_fn(a)
            cl_mean, st = _cl_st(rd)
            rows.append(dict(config=label, alpha=a, cl_mean=cl_mean, strouhal=st))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"dx=0.02 (baseline)": "#c0392b", "dx=0.01": "#2980b9",
              "ngrid=2": "#16a085", "ngrid=3": "#8e44ad"}
    for label, _ in configs:
        sub = [r for r in rows if r["config"] == label and r["cl_mean"] is not None]
        sub.sort(key=lambda r: r["alpha"])
        if not sub:
            continue
        a = [r["alpha"] for r in sub]
        axes[0].plot(a, [r["cl_mean"] for r in sub], "o-", color=colors[label], label=label, lw=1.6)
        axes[1].plot(a, [r["strouhal"] for r in sub], "o-", color=colors[label], label=label, lw=1.6)
    axes[0].set_ylabel(r"mean $\overline{C_l}$"); axes[0].set_title("Does the Cl jump's location/size move?")
    axes[1].set_ylabel("Strouhal number St"); axes[1].set_title("Does the Strouhal jump's location/size move?")
    for ax in axes:
        ax.set_xlabel("mean angle of attack (deg)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        ax.set_xticks(TRANS_ANGLES)
    fig.suptitle("Test 1c: is the 34-35deg transition itself grid/domain-sensitive?\n"
                 "(directly testing what Test 1b's justification assumed but never checked)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    out = FIGS / "test_1c_transition_sensitivity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out.name}")

    import csv
    with open(DATA / "test1c_transition_sensitivity.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["config", "alpha", "cl_mean", "strouhal"])
        w.writeheader(); w.writerows(rows)

    print("\nTest 1c summary:")
    for label, _ in configs:
        sub = [r for r in rows if r["config"] == label]
        sub.sort(key=lambda r: r["alpha"])
        print(f"  {label}: " + ", ".join(
            f"a{r['alpha']}: Cl={r['cl_mean']:.3f} St={r['strouhal']:.3f}"
            if r["cl_mean"] is not None else f"a{r['alpha']}: (missing)" for r in sub))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        do_run()
    else:
        analyze()
