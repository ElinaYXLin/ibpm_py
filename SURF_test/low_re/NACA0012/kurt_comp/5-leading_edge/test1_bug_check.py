"""
test1_bug_check.py

Group 1 of the LE/TE striping investigation (see README.md): "is it a bug?"
Uses ONLY the already-existing steady, dx=0.02, Re=1000 runs in
../1-paper_based/runs/dx0.020/steady_{py,cpp}_a{00,09,12}/ -- zero new
solver runs.

Test 1a -- py_static vs cpp_static difference field at the LE/TE. Two
independently-written implementations that agree to floating-point
precision at the LE/TE rules out a port-specific bug there; it would have
to live in the shared algorithm/formulation or the shared geometry input
(the repo has already verified src/ is byte-identical to py/'s origin --
see static_test/README.md -- so a shared-algorithm bug would be an
upstream one, i.e. effectively a method limitation, not a coding bug).

Test 1b -- time evolution of the LE-ringing / TE-blob amplitude across the
existing t=0,2.5,...,30 snapshots. A genuine numerical instability grows
unboundedly (or blows up) in time; a fixed discretization artifact
saturates to a bounded, ~steady pattern once the flow itself has settled.

Usage: python3 test1_bug_check.py
Output: figures/test1a_*.png, figures/test1b_*.png, data/test1*.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

ALPHAS = [0, 9, 12]
STEP_FINAL = 3000
STEPS_ALL = list(range(0, 3001, 250))  # t=0,2.5,...,30 (restart cadence = nsteps//12)
DX = 0.02
LE_XLIM, LE_YLIM = (-0.15, 0.35), (-0.25, 0.25)
TE_XLIM, TE_YLIM = (0.65, 1.35), (-0.25, 0.25)


def run_dir(impl, alpha):
    return c.KURT1 / "runs" / "dx0.020" / f"steady_{impl}_a{alpha:02d}"


def window_amplitude(X, Y, om, xlim, ylim):
    """max|omega| - min|omega|... no: use peak-to-trough like test0, but
    over the full 2-D window (not just a 1-D lineout) so it's meaningful
    at every timestep, including early transients that aren't yet
    left-right symmetric about y=0."""
    _, _, sub, _, _ = c.window(X, Y, om, xlim, ylim)
    return float(sub.max() - sub.min()), float(np.abs(sub).max())


def main():
    X, Y, dx = c.grid_xy(DX)

    # ---------------- Test 1a: py vs cpp diff at LE/TE, final (developed) snapshot ----------------
    rows_1a = []
    for alpha in ALPHAS:
        om_py = c.load_omega(run_dir("py", alpha), STEP_FINAL)
        om_cpp = c.load_omega(run_dir("cpp", alpha), STEP_FINAL)
        diff = om_py - om_cpp

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
        for col, (xlim, ylim, name) in enumerate([(LE_XLIM, LE_YLIM, "LE"),
                                                     (TE_XLIM, TE_YLIM, "TE")]):
            Xw, Yw, fw_py, _, _ = c.window(X, Y, om_py, xlim, ylim)
            _, _, fw_cpp, _, _ = c.window(X, Y, om_cpp, xlim, ylim)
            fw_diff = fw_py - fw_cpp
            ax = axes[col]
            vmax = max(np.abs(fw_diff).max(), 1e-300)
            im = ax.pcolormesh(Xw, Yw, fw_diff[:-1, :-1], shading="flat",
                                cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            ax.set_aspect("equal")
            ax.set_title(f"{name}: py-cpp diff, max|diff|={vmax:.2e}\n"
                         f"(max|omega| there = {max(np.abs(fw_py).max(), np.abs(fw_cpp).max()):.3g})",
                         fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.8)
            rows_1a.append(dict(alpha=alpha, region=name,
                                 max_abs_diff=vmax,
                                 max_abs_omega=max(np.abs(fw_py).max(), np.abs(fw_cpp).max()),
                                 relative_diff=vmax / max(np.abs(fw_py).max(), np.abs(fw_cpp).max(), 1e-300)))
        # third panel: whole near-body field diff for context
        Xw, Yw, fw_py, _, _ = c.window(X, Y, om_py, (-1, 4), (-1.2, 1.2))
        _, _, fw_cpp, _, _ = c.window(X, Y, om_cpp, (-1, 4), (-1.2, 1.2))
        fw_diff = fw_py - fw_cpp
        ax = axes[2]
        vmax = max(np.abs(fw_diff).max(), 1e-300)
        im = ax.pcolormesh(Xw, Yw, fw_diff[:-1, :-1], shading="flat", cmap="RdBu_r",
                            vmin=-vmax, vmax=vmax)
        ax.set_aspect("equal")
        ax.set_title(f"whole near-body field, max|diff|={vmax:.2e}", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8)

        fig.suptitle(f"Test 1a: py_static - cpp_static, alpha={alpha}deg steady, "
                     f"Re=1000, dx=0.02, t=30", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.9])
        out = c.FIGS / f"test1a_py_cpp_diff_a{alpha:02d}.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"wrote {out.name}")

    with open(c.DATA / "test1a_py_cpp_diff.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["alpha", "region", "max_abs_diff", "max_abs_omega", "relative_diff"])
        w.writeheader()
        w.writerows(rows_1a)
    print("\nTest 1a summary:")
    for r in rows_1a:
        print(f"  alpha={r['alpha']:>2} {r['region']}: max|py-cpp|={r['max_abs_diff']:.3e}, "
              f"relative to max|omega|={r['relative_diff']:.3e}")

    # ---------------- Test 1b: amplitude vs time ----------------
    rows_1b = []
    for alpha in ALPHAS:
        for impl in ("py", "cpp"):
            rd = run_dir(impl, alpha)
            for step in STEPS_ALL:
                if not (rd / f"flow{step:05d}.bin").exists():
                    continue
                om = c.load_omega(rd, step)
                amp_le, peak_le = window_amplitude(X, Y, om, LE_XLIM, LE_YLIM)
                amp_te, peak_te = window_amplitude(X, Y, om, TE_XLIM, TE_YLIM)
                t = step * 0.01  # dt=0.01 for steady dx=0.02 cases
                rows_1b.append(dict(alpha=alpha, impl=impl, step=step, t=t,
                                     amp_le=amp_le, peak_le=peak_le,
                                     amp_te=amp_te, peak_te=peak_te))

    with open(c.DATA / "test1b_amplitude_vs_time.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["alpha", "impl", "step", "t", "amp_le", "peak_le", "amp_te", "peak_te"])
        w.writeheader()
        w.writerows(rows_1b)
    print(f"\nwrote {c.DATA / 'test1b_amplitude_vs_time.csv'}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    colors = {0: "#2980b9", 9: "#c0392b", 12: "#16a085"}
    for alpha in ALPHAS:
        for impl, ls in (("py", "-"), ("cpp", "--")):
            sub = [r for r in rows_1b if r["alpha"] == alpha and r["impl"] == impl]
            sub.sort(key=lambda r: r["t"])
            tt = [r["t"] for r in sub]
            axes[0].plot(tt, [r["amp_le"] for r in sub], ls, color=colors[alpha],
                         label=f"a{alpha} {impl}", lw=1.3, marker=".")
            axes[1].plot(tt, [r["amp_te"] for r in sub], ls, color=colors[alpha],
                         label=f"a{alpha} {impl}", lw=1.3, marker=".")
    axes[0].set_title("LE window peak-to-trough amplitude vs time"); axes[0].set_xlabel("t")
    axes[1].set_title("TE window peak-to-trough amplitude vs time"); axes[1].set_xlabel("t")
    for ax in axes:
        ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Test 1b: LE/TE artifact amplitude over time -- growing (instability) "
                 "vs. saturating (discretization signature)?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = c.FIGS / "test1b_amplitude_vs_time.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
