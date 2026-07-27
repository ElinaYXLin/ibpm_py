"""
test5_cylinder_control.py

Test 5: cylinder control. The LE/TE striping is reported to vanish on a
cylinder (gentle curvature everywhere, no sharp nose) -- if upstream
noise also vanishes there, sharp curvature is necessary for the
artifact and this folder's null-region premise ("the true answer is
zero upstream") is empirically validated, not just argued from a
diffusion-length estimate. If the cylinder still shows upstream noise,
this is a generic immersed-boundary artifact, not curvature-specific.
Recomputes Tests 1-4's metrics for the cylinder at dx=0.02, Re=1000,
alpha=0, steady, both implementations, alongside the NACA0012 baseline
for direct comparison. Zero new runs -- reuses
../5-leading_edge/runs/shape_spacing/cylinder(_cpp).

Usage: python3 test5_cylinder_control.py
Output: figures/test5_cylinder_control.png, data/test5_cylinder_control.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c
from test4_spectral import nyquist_fraction_2d

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0
DX = 0.02
NSTEPS = 3000

CASES = {
    "naca0012": dict(geom=c.BASE_GEOM_DX002,
                      py=c.KURT1 / "runs" / "dx0.020" / "steady_py_a00",
                      cpp=c.KURT1 / "runs" / "dx0.020" / "steady_cpp_a00"),
    "cylinder": dict(geom=c.CYLINDER_GEOM_DX002,
                      py=c.KURT5 / "runs" / "shape_spacing" / "cylinder",
                      cpp=c.KURT5 / "runs" / "shape_spacing" / "cylinder_cpp"),
}
SHAPE_COLOR = {"naca0012": "#1f77b4", "cylinder": "#e67e22"}


def main():
    X, Y = c.grid_xy(DX)
    rows = []
    profiles = {}
    for shape, case in CASES.items():
        g = c.load_geom_points(case["geom"])
        x_le = g["x"][g["i_le"]]
        profiles[shape] = {}
        for impl in ("py", "cpp"):
            om = c.load_omega(case[impl], NSTEPS)
            m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            prof = c.upstream_profile(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            mask = c.upstream_mask(X, x_le, BUFFER_DX, DX)
            ys = Y[0, :]
            iy = np.where((ys >= YLIM[0]) & (ys <= YLIM[1]))[0]
            sub = om[np.ix_(np.where(mask)[0], iy)]
            frac = nyquist_fraction_2d(sub)
            profiles[shape][impl] = prof
            rows.append(dict(shape=shape, impl=impl, x_le=x_le, enstrophy=m["enstrophy"],
                              int_abs=m["int_abs"], peak=m["peak"], nyquist_fraction=frac))
            print(f"{shape} {impl}: x_le={x_le:.4f} enstrophy={m['enstrophy']:.5f} "
                  f"int_abs={m['int_abs']:.5f} peak={m['peak']:.3f} nyquist_frac={frac:.3f}")

    with open(c.DATA / "test5_cylinder_control.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["shape", "impl", "x_le", "enstrophy", "int_abs", "peak", "nyquist_fraction"])
        w.writeheader(); w.writerows(rows)
    print("wrote test5_cylinder_control.csv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    for shape in CASES:
        x_le = next(r["x_le"] for r in rows if r["shape"] == shape and r["impl"] == "py")
        for impl, ls in (("py", "-"), ("cpp", "--")):
            prof = profiles[shape][impl]
            order = np.argsort(prof["xs"])
            dist = x_le - prof["xs"][order]
            vals = np.clip(prof["max_abs"][order], 1e-6, None)
            lbl = shape if impl == "py" else None
            ax.plot(dist, vals, ls, color=SHAPE_COLOR[shape], label=lbl,
                    lw=1.6 if impl == "py" else 1.0, alpha=1.0 if impl == "py" else 0.6)
    ax.set_yscale("log")
    ax.set_xlabel("distance upstream of LE (chord)"); ax.set_ylabel("max|omega| over y")
    ax.set_title("Streamwise decay: NACA0012 vs cylinder", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=9)

    ax = axes[1]
    shapes = list(CASES)
    xpos = np.arange(len(shapes))
    py_ens = [next(r["enstrophy"] for r in rows if r["shape"] == s and r["impl"] == "py") for s in shapes]
    cpp_ens = [next(r["enstrophy"] for r in rows if r["shape"] == s and r["impl"] == "cpp") for s in shapes]
    width = 0.35
    ax.bar(xpos - width / 2, py_ens, width, label="py_static", color="#1f77b4")
    ax.bar(xpos + width / 2, cpp_ens, width, label="cpp_static", color="#d62728", alpha=0.7)
    ax.set_xticks(xpos); ax.set_xticklabels(shapes)
    ax.set_yscale("log")
    ax.set_title("Upstream enstrophy: NACA0012 vs cylinder", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both", axis="y")

    ax = axes[2]
    for shape in CASES:
        py_v = np.array([r["enstrophy"] for r in rows if r["shape"] == shape and r["impl"] == "py"])
        cpp_v = np.array([r["enstrophy"] for r in rows if r["shape"] == shape and r["impl"] == "cpp"])
        reldiff = np.abs(py_v - cpp_v) / np.maximum(py_v, 1e-12)
        ax.bar(shape, reldiff[0], color=SHAPE_COLOR[shape])
    ax.set_yscale("log")
    ax.set_ylabel("relative diff, py vs cpp enstrophy")
    ax.set_title("py/cpp agreement check", fontsize=10)
    ax.grid(alpha=0.3, which="both", axis="y")

    fig.suptitle("Test 5: cylinder control -- is sharp curvature necessary for upstream noise?\n"
                 "Re=1000, alpha=0, steady, dx=0.02", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = c.FIGS / "test5_cylinder_control.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
