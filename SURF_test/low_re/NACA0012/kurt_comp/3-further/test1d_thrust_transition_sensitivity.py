"""
test1d_thrust_transition_sensitivity.py

Separate from Test 1c: does grid/domain refinement shrink IBPM's
extended thrust regime back toward the paper's reported cutoff?

`../1-paper_based`'s thrust_check found instantaneous Cd<0 (thrust) for
3<=alpha0<=~50deg at f=4Hz in IBPM, vs. the paper's reported 3-37/38deg
-- the single largest quantitative disagreement in that comparison
(onset matches almost exactly; the *upper* cutoff does not). Unlike
Test 1c (steady, non-pitching runs -- can't produce thrust, since the
effect is driven by the pitching motion itself, not steady shedding),
this test reruns the ACTUAL oscillating f4hz pitching motion (pitchplunge,
1deg amplitude, same motion line as ../1-paper_based/geom/
naca0012_dx0.0200_f4hz.geom) at dx=0.01 and ngrid=2/3, at alpha0=
35,40,45,50 -- 35 is inside the range both IBPM and the paper agree
shows thrust; 40/45/50 are inside IBPM's baseline thrust regime but
outside the paper's claimed cutoff, i.e. exactly the disputed range.

Same convention as Test 1c: py_static only, Re=1000, dx=0.02/ngrid=1
baseline reuses ../1-paper_based/runs/dx0.020/f4hz_py_a{35,40,45,50}
directly (zero new runs for that row); dx=0.01 and ngrid=2/3 are new
(cheap relative to the faithful2 runs elsewhere: same duration/dt as
the existing f4hz baseline, t=30, dt=0.005).

Usage:
  python3 test1d_thrust_transition_sensitivity.py run       # launch new runs
  python3 test1d_thrust_transition_sensitivity.py analyze   # analyze + plot
Output: runs/thrust_check/, figures/test_1d_thrust_transition_sensitivity.png,
        data/test1d_thrust_transition_sensitivity.csv
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_further as rf

ANGLES = [35, 40, 45, 50]
FIGS = rf.FURTHER / "figures"
DATA = rf.FURTHER / "data"
GEOM_DX002_F4HZ = rf.PAPER_RUNS.parent.parent / "geom" / "naca0012_dx0.0200_f4hz.geom"
GEOM_DX001_F4HZ = rf.FURTHER / "geom" / "naca0012_dx0.0100_f4hz.geom"
DT, NSTEPS = 0.005, 6000  # matches the existing f4hz baseline's temporal resolution exactly


def build_jobs():
    jobs = []
    for a in ANGLES:
        outdir = rf.RUNS / "thrust_check" / "dx_refine" / f"dx0.010_a{a:02d}"
        jobs.append(dict(outdir=outdir, geom=GEOM_DX001_F4HZ, nx=600, ny=300, ngrid=1,
                          domain=rf.BASE_DOMAIN, alpha=a, Re=rf.RE, dt=DT, nsteps=NSTEPS))
    ngrid_domain = dict(length=6, xoffset=-2, yoffset=-1.52)
    for ngrid in (2, 3):
        for a in ANGLES:
            outdir = rf.RUNS / "thrust_check" / "ngrid_sweep" / f"ngrid{ngrid}_a{a:02d}"
            jobs.append(dict(outdir=outdir, geom=GEOM_DX002_F4HZ, nx=300, ny=152, ngrid=ngrid,
                              domain=ngrid_domain, alpha=a, Re=rf.RE, dt=DT, nsteps=NSTEPS))
    return jobs


def do_run():
    jobs = build_jobs()
    print(f"{len(jobs)} jobs", flush=True)
    for j in jobs:
        name, status, elapsed = rf.run_case(**j)
        print(f"  {name}: {status} ({elapsed:.0f}s)", flush=True)


def _cd_min(run_dir):
    f = run_dir / "flow.force"
    if not f.exists():
        return None
    d = np.loadtxt(f)
    n = len(d)
    seg = d[int(n * 0.5):]
    return float(np.min(seg[:, 2]))


CONFIGS = [
    ("dx=0.02 (baseline)", lambda a: rf.PAPER_RUNS / f"f4hz_py_a{a:02d}"),
    ("dx=0.01", lambda a: rf.RUNS / "thrust_check" / "dx_refine" / f"dx0.010_a{a:02d}"),
    ("ngrid=2", lambda a: rf.RUNS / "thrust_check" / "ngrid_sweep" / f"ngrid2_a{a:02d}"),
    ("ngrid=3", lambda a: rf.RUNS / "thrust_check" / "ngrid_sweep" / f"ngrid3_a{a:02d}"),
]
COLORS = {"dx=0.02 (baseline)": "#c0392b", "dx=0.01": "#2980b9",
          "ngrid=2": "#16a085", "ngrid=3": "#8e44ad"}


def analyze():
    rows = []
    for label, path_fn in CONFIGS:
        for a in ANGLES:
            rd = path_fn(a)
            rows.append(dict(config=label, alpha=a, cd_min=_cd_min(rd)))

    with open(DATA / "test1d_thrust_transition_sensitivity.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["config", "alpha", "cd_min"])
        w.writeheader(); w.writerows(rows)
    print("wrote test1d_thrust_transition_sensitivity.csv")

    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    for label, _ in CONFIGS:
        sub = [r for r in rows if r["config"] == label and r["cd_min"] is not None]
        sub.sort(key=lambda r: r["alpha"])
        if not sub:
            continue
        a = [r["alpha"] for r in sub]
        ax.plot(a, [r["cd_min"] for r in sub], "o-", color=COLORS[label], label=label, lw=1.8, ms=7)
    ax.axhline(0, color="black", lw=1.0, ls="--", alpha=0.6)
    ax.axvline(37.5, color="gray", lw=1.0, ls=":", alpha=0.7)
    ax.text(37.7, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.05,
            "paper's\nreported\ncutoff", fontsize=7, va="top", color="gray")
    ax.set_xlabel("mean angle of attack alpha0 (deg)")
    ax.set_ylabel(r"minimum instantaneous $C_d$ over cycle (f=4Hz pitching)")
    ax.set_title("Does IBPM's extended thrust regime (min-Cd<0 past 37-38deg)\n"
                  "shrink under grid/domain refinement?")
    ax.set_xticks(ANGLES)
    ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle("Test 1d: thrust_check companion -- is the disagreement with the\n"
                 "paper's 37-38deg cutoff a resolution artifact?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = FIGS / "test_1d_thrust_transition_sensitivity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out.name}")

    print("\nTest 1d min-Cd (thrust) summary:")
    for label, _ in CONFIGS:
        sub = [r for r in rows if r["config"] == label]
        sub.sort(key=lambda r: r["alpha"])
        print(f"  {label}: " + ", ".join(
            f"a{r['alpha']}: min(Cd)={r['cd_min']:.4f}" if r["cd_min"] is not None
            else f"a{r['alpha']}: (missing)" for r in sub))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        do_run()
    else:
        analyze()
