"""
test1c_mindrag_transition_sensitivity.py

Minimum-Cd companion to test1c_drag_transition_sensitivity.py: that
script tracked the 34-35deg transition's effect on *mean* Cd. But
Kurtulus's own thrust reporting (Fig 1, and the thrust_check discussion
in ../1-paper_based/README.md) is about the minimum instantaneous Cd
over a cycle -- the most-negative/least-drag point, not the time-mean
-- so a mean-Cd comparison isn't actually the quantity she reported.
This reuses the exact same runs as test1c_drag_transition_sensitivity.py
(zero new runs) and just changes the reduction from mean(Cd) to min(Cd)
over the same developed (last-50%) segment.

**Important caveat, discovered by running this**: at alpha=33-36deg,
steady (non-pitching), min(Cd) never actually goes negative in any of
the 4 configs -- there is no real thrust here, only a smaller-than-mean
drag at the low point of the vortex-shedding cycle. Actual thrust
(Cd<0) in this repo only shows up for the oscillating f=4Hz pitching
motion, which is what test1d_thrust_transition_sensitivity.py tests
directly. This script is still useful as the min-Cd analog of test1c's
mean-Cd check (does the 34-35deg transition's location/size in min-Cd
move under refinement, same question test1c asks for the mean) -- it
just isn't a thrust/max-thrust measurement for this particular set of
runs, and is labeled accordingly rather than as "max thrust."

Usage: python3 test1c_mindrag_transition_sensitivity.py
Output: figures/test_1c_mindrag_transition_sensitivity.png,
        data/test1c_mindrag_transition_sensitivity.csv
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


def _cd_min(run_dir):
    f = run_dir / "flow.force"
    if not f.exists():
        return None
    d = np.loadtxt(f)
    n = len(d)
    seg = d[int(n * 0.5):]
    return float(np.min(seg[:, 2]))


def main():
    rows = []
    for label, path_fn in CONFIGS:
        for a in TRANS_ANGLES:
            rd = path_fn(a)
            rows.append(dict(config=label, alpha=a, cd_min=_cd_min(rd)))

    with open(DATA / "test1c_mindrag_transition_sensitivity.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["config", "alpha", "cd_min"])
        w.writeheader(); w.writerows(rows)
    print("wrote test1c_mindrag_transition_sensitivity.csv")

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for label, _ in CONFIGS:
        sub = [r for r in rows if r["config"] == label and r["cd_min"] is not None]
        sub.sort(key=lambda r: r["alpha"])
        if not sub:
            continue
        a = [r["alpha"] for r in sub]
        ax.plot(a, [r["cd_min"] for r in sub], "o-", color=COLORS[label], label=label, lw=1.8, ms=7)
    ax.axhline(0, color="black", lw=1.0, ls="--", alpha=0.6)
    ax.text(TRANS_ANGLES[0], 0.03, "Cd=0 (thrust below this line -- not reached here)",
            fontsize=7, color="#555555", va="bottom")
    ax.set_xlabel("mean angle of attack (deg)")
    ax.set_ylabel(r"minimum instantaneous $C_d$ over cycle")
    ax.set_title("Does the min-Cd jump's location/size move under grid/domain refinement?")
    ax.set_xticks(TRANS_ANGLES)
    ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle("Test 1c (min-drag companion): is the 34-35deg transition\n"
                 "grid/domain-sensitive in min-Cd, like it is in mean-Cl and Strouhal?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = FIGS / "test_1c_mindrag_transition_sensitivity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out.name}")

    print("\nTest 1c min-Cd (max-thrust) summary:")
    for label, _ in CONFIGS:
        sub = [r for r in rows if r["config"] == label]
        sub.sort(key=lambda r: r["alpha"])
        print(f"  {label}: " + ", ".join(
            f"a{r['alpha']}: min(Cd)={r['cd_min']:.4f}" if r["cd_min"] is not None
            else f"a{r['alpha']}: (missing)" for r in sub))


if __name__ == "__main__":
    main()
