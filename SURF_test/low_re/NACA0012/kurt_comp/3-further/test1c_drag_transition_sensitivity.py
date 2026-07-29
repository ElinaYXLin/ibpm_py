"""
test1c_drag_transition_sensitivity.py

Companion to test1c_transition_sensitivity.py: that script checked
whether the 34-35deg jump in mean Cl and Strouhal number moves or
vanishes under grid/domain refinement (dx=0.01, ngrid=2, ngrid=3 vs.
the dx=0.02/ngrid=1 baseline). This does the same check for mean drag
Cd, reusing the exact same runs (zero new runs) -- same 4 configs, same
alpha=33,34,35,36 straddling the transition.

(For minimum instantaneous Cd over an oscillating pitching cycle -- the
thrust_check-relevant quantity -- see test1d_thrust_transition_sensitivity.py,
a separate test: these are steady, non-pitching runs, so they can't
probe the pitching-motion-driven thrust phenomenon that test1d checks.)

Usage: python3 test1c_drag_transition_sensitivity.py
Output: figures/test_1c_drag_transition_sensitivity.png,
        data/test1c_drag_transition_sensitivity.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_further as rf
from test1c_transition_sensitivity import TRANS_ANGLES, FIGS, DATA

CONFIGS = [
    ("dx=0.02 (baseline)", lambda a: rf.PAPER_RUNS / f"steady_py_a{a:02d}"),
    ("dx=0.01", lambda a: rf.RUNS / "transition_check" / "dx_refine" / f"dx0.010_a{a:02d}"),
    ("ngrid=2", lambda a: rf.RUNS / "transition_check" / "ngrid_sweep" / f"ngrid2_a{a:02d}"),
    ("ngrid=3", lambda a: rf.RUNS / "transition_check" / "ngrid_sweep" / f"ngrid3_a{a:02d}"),
]
COLORS = {"dx=0.02 (baseline)": "#c0392b", "dx=0.01": "#2980b9",
          "ngrid=2": "#16a085", "ngrid=3": "#8e44ad"}


def _cd(run_dir):
    f = run_dir / "flow.force"
    if not f.exists():
        return None
    d = np.loadtxt(f)
    n = len(d)
    seg = d[int(n * 0.5):]
    cd = seg[:, 2]
    return float(np.mean(cd))


def main():
    rows = []
    for label, path_fn in CONFIGS:
        for a in TRANS_ANGLES:
            rd = path_fn(a)
            rows.append(dict(config=label, alpha=a, cd_mean=_cd(rd)))

    with open(DATA / "test1c_drag_transition_sensitivity.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["config", "alpha", "cd_mean"])
        w.writeheader(); w.writerows(rows)
    print("wrote test1c_drag_transition_sensitivity.csv")

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for label, _ in CONFIGS:
        sub = [r for r in rows if r["config"] == label and r["cd_mean"] is not None]
        sub.sort(key=lambda r: r["alpha"])
        if not sub:
            continue
        a = [r["alpha"] for r in sub]
        ax.plot(a, [r["cd_mean"] for r in sub], "o-", color=COLORS[label], label=label, lw=1.8, ms=7)
    ax.set_xlabel("mean angle of attack (deg)")
    ax.set_ylabel(r"mean $\overline{C_d}$")
    ax.set_title("Does the Cd jump's location/size move under grid/domain refinement?")
    ax.set_xticks(TRANS_ANGLES)
    ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle("Test 1c (drag companion): is the 34-35deg Cd transition\n"
                 "grid/domain-sensitive, like Cl and Strouhal?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = FIGS / "test_1c_drag_transition_sensitivity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out.name}")

    print("\nTest 1c drag summary:")
    for label, _ in CONFIGS:
        sub = [r for r in rows if r["config"] == label]
        sub.sort(key=lambda r: r["alpha"])
        print(f"  {label}: " + ", ".join(
            f"a{r['alpha']}: Cd={r['cd_mean']:.4f}" if r["cd_mean"] is not None
            else f"a{r['alpha']}: (missing)" for r in sub))


if __name__ == "__main__":
    main()
