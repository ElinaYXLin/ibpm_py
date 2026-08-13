"""
plot_field_only.py

Field-max-only versions of all 5 figures in this folder. compute_fieldmax.py's
figures always show field-max side by side with the (unreliable) lineout-max
for direct comparison; this strips the lineout bars out entirely, leaving
just the 2-D field-max metric this folder treats as ground truth -- useful
whenever the lineout comparison isn't the point and the field-max numbers
alone are what's being presented (e.g. to the mentor).

Reads the same data/*.csv this folder's compute_fieldmax.py / run_ngrid_sweep.py
already produced -- no new computation, just a different plot.

Usage: python3 plot_field_only.py
Output: figures/field_only/*.png
"""
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures" / "field_only"
FIGS.mkdir(parents=True, exist_ok=True)


def load_csv(name):
    return np.genfromtxt(DATA / name, delimiter=",", names=True, dtype=None, encoding="utf-8")


def bar_fieldmax_only(rows, key, fname, title, group_order=None):
    keys = group_order or sorted(set(rows[key]), key=str)
    fig, axes = plt.subplots(1, 2, figsize=(max(7, 1.6 * len(keys)), 5))
    for ax, region in zip(axes, ("LE", "TE")):
        vals = [float(rows["field_max"][(rows[key] == k) & (rows["region"] == region)
                                         & (rows["impl"] == "py")][0])
                for k in keys]
        x = np.arange(len(keys))
        ax.bar(x, vals, 0.5, color="#c0392b")
        for xi, v in zip(x, vals):
            ax.text(xi, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels([str(k) for k in keys], rotation=30, ha="right", fontsize=8)
        ax.set_title(f"{region} window", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        if ax is axes[0]:
            ax.set_ylabel("2-D field-max |omega|")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(FIGS / fname, dpi=140)
    plt.close(fig)
    print(f"wrote field_only/{fname}")


def to_str_key(arr):
    return np.array([v.decode() if isinstance(v, bytes) else v for v in arr])


def main():
    r1 = load_csv("table1_test0b_angles.csv")
    r1_key = r1["alpha"]
    bar_fieldmax_only({"field_max": r1["field_max"], "region": to_str_key(r1["region"]),
                        "impl": to_str_key(r1["impl"]), "alpha": r1_key},
                       "alpha", "table1_test0b_angles.png",
                       "Test 0b companion (2-D field-max only): alpha=0/9/12",
                       group_order=[0, 9, 12])

    r2 = load_csv("table2_test2a_grid_refinement.csv")
    bar_fieldmax_only({"field_max": r2["field_max"], "region": to_str_key(r2["region"]),
                        "impl": to_str_key(r2["impl"]), "dx": r2["dx"]},
                       "dx", "table2_test2a_grid_refinement.png",
                       "Test 2a companion (2-D field-max only): grid refinement",
                       group_order=[0.02, 0.01, 0.005])

    r3 = load_csv("table3_test3a_spacing.csv")
    bar_fieldmax_only({"field_max": r3["field_max"], "region": to_str_key(r3["region"]),
                        "impl": to_str_key(r3["impl"]), "case": to_str_key(r3["case"])},
                       "case", "table3_test3a_spacing.png",
                       "Test 3a companion (2-D field-max only): boundary-point spacing",
                       group_order=["naca0012_LTEsparse (ds=4dx)", "naca0012_baseline (ds=dx)",
                                    "naca0012_LTEdense (ds=dx/4)"])

    r4 = load_csv("table4_test3b_shape.csv")
    bar_fieldmax_only({"field_max": r4["field_max"], "region": to_str_key(r4["region"]),
                        "impl": to_str_key(r4["impl"]), "case": to_str_key(r4["case"])},
                       "case", "table4_test3b_shape.png",
                       "Test 3b companion (2-D field-max only): shape/thickness family",
                       group_order=["naca0006", "naca0012 (baseline)", "naca0018",
                                    "naca0012_roundTE", "cylinder"])

    r5 = load_csv("table5_ngrid_sweep.csv")
    bar_fieldmax_only({"field_max": r5["field_max"], "region": to_str_key(r5["region"]),
                        "impl": to_str_key(r5["impl"]), "ngrid": r5["ngrid"]},
                       "ngrid", "table5_ngrid_sweep.png",
                       "Table 5 (2-D field-max only): ngrid=1,2,3,4 (dx=0.02, ds=dx, alpha=0 fixed)",
                       group_order=[1, 2, 3, 4])


if __name__ == "__main__":
    main()
