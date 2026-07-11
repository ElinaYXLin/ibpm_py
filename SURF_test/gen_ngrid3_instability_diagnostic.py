"""
gen_ngrid3_instability_diagnostic.py

Follow-up to 3-ngrid=3/ and 3-dx0.01/: explains (1) why ngrid=3 diverges
around t~20-21 regardless of dt, and (2) why dx=0.01 shows larger peak
vorticity fluctuations than dx=0.02, using the restart snapshots already
on disk (no new solver runs).

Panel A: domain-RMS vorticity vs time, ngrid=1 baseline (saturates/
plateaus, as expected for a statistically-stationary vortex-shedding wake)
vs. ngrid=3 (grows without bound until blowup) -- shows the ngrid=3 growth
is a real runaway instability, not just noisier output.

Panel B: for ngrid=3, RMS vorticity split into a "core" band (near the
body/wake, |x|<1.5, |y|<0.5) vs. an "edge" band (outermost 3 cells of the
finest grid, i.e. right at the fine/coarse domain interface). The edge
band starts near zero (as expected -- it's far-field) and grows to match
the core by the time the run blows up, showing the runaway energy growth
originates at the multi-domain grid interface, not in the wake itself.

Panel C: max|omega| and RMS vorticity vs dx (0.02 vs 0.01, both ngrid=1,
matched t=30 snapshot) -- shows max|omega| grows substantially with
resolution while RMS stays roughly flat, i.e. finer dx resolves sharper/
more concentrated peaks in the same broadband noise, not more total
energy. Consistent with SD7003/README.md's documented explanation
(under-resolved boundary layer at Re~61100, no subgrid dissipation) rather
than a new artifact introduced by refining the grid.

Usage: python3 SURF_test/gen_ngrid3_instability_diagnostic.py
Output: SURF_test/airfoils/SD7003/3-ngrid=3/instability_diagnostic.png
"""
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
SURF = REPO / "SURF_test"
sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402


def load(path):
    return State(filename=str(path)).omega._data[0].copy()


def rms(f):
    return float(np.sqrt(np.nanmean(f.astype(np.float64) ** 2)))


def main():
    # ---- Panel A data: ngrid=1 baseline vs ngrid=3 RMS(t) ----
    base_dir = SURF / "airfoils" / "SD7003" / "_run_data" / "flowfield"
    base_steps = list(range(0, 3001, 500))
    base_t = [s * 0.01 for s in base_steps]
    base_rms = [rms(load(base_dir / f"flow{s:05d}.bin")) for s in base_steps]

    ng3_dir = SURF / "airfoils" / "SD7003" / "_run_data_ngrid3" / "flowfield"
    ng3_steps = [s for s in range(0, 6001, 1000)
                 if (ng3_dir / f"flow{s:05d}.bin").exists()]
    ng3_t, ng3_rms, ng3_finite = [], [], []
    for s in ng3_steps:
        f = load(ng3_dir / f"flow{s:05d}.bin")
        ng3_t.append(s * 0.005)
        ng3_finite.append(bool(np.isfinite(f).all()))
        ng3_rms.append(rms(f) if np.isfinite(f).all() else np.nan)

    # ---- Panel B data: ngrid=3 core-band vs edge-band RMS(t) (finite steps only) ----
    NX, NY = 300, 152
    dx = 6.0 / NX
    xs = -2.0 + np.arange(1, NX) * dx
    ys = -1.52 + np.arange(1, NY) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    edge = 3
    ii = np.arange(X.shape[0])
    jj = np.arange(X.shape[1])
    mask_edge = ((ii[:, None] < edge) | (ii[:, None] >= X.shape[0] - edge) |
                 (jj[None, :] < edge) | (jj[None, :] >= X.shape[1] - edge))
    mask_core = (np.abs(X) < 1.5) & (np.abs(Y) < 0.5)

    core_t, core_rms, edge_rms = [], [], []
    for s in ng3_steps:
        f = load(ng3_dir / f"flow{s:05d}.bin")
        if not np.isfinite(f).all():
            continue
        core_t.append(s * 0.005)
        core_rms.append(rms(f[mask_core]))
        edge_rms.append(rms(f[mask_edge]))

    # ---- Panel C data: max|omega|, RMS vs dx at matched t=30 ----
    dx02 = load(SURF / "airfoils" / "SD7003" / "_run_data" / "flowfield" / "flow03000.bin")
    dx01 = load(SURF / "airfoils" / "SD7003" / "_run_data_dx001" / "flowfield" / "flow06000.bin")
    dxs = [0.02, 0.01]
    max_vals = [float(np.max(np.abs(dx02))), float(np.max(np.abs(dx01)))]
    rms_vals = [rms(dx02), rms(dx01)]

    # ---- Plot ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    ax.plot(base_t, base_rms, "o-", color="C0", label="ngrid=1 baseline (dt=0.01)")
    ax.plot(ng3_t, ng3_rms, "s-", color="C3", label="ngrid=3 (dt=0.005)")
    for t, isf in zip(ng3_t, ng3_finite):
        if not isf:
            ax.axvline(t, color="C3", ls=":", alpha=0.4)
    ax.set_xlabel("t"); ax.set_ylabel(r"domain-RMS $|\omega|$")
    ax.set_title("A: ngrid=1 saturates, ngrid=3 grows unbounded")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(core_t, core_rms, "o-", color="C2", label="core (near body/wake)")
    ax.plot(core_t, edge_rms, "^-", color="C1", label="edge (fine-grid boundary)")
    ax.set_xlabel("t"); ax.set_ylabel("RMS $|\\omega|$ (banded)")
    ax.set_title("B: ngrid=3 -- edge-band energy catches up to core\n(instability starts at the domain interface)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[2]
    ax2 = ax.twinx()
    ax.plot(dxs, max_vals, "o-", color="C3", label="max$|\\omega|$")
    ax2.plot(dxs, rms_vals, "s-", color="C0", label="RMS $|\\omega|$")
    ax.set_xlabel("dx"); ax.set_ylabel("max$|\\omega|$", color="C3")
    ax2.set_ylabel("RMS $|\\omega|$", color="C0")
    ax.invert_xaxis()
    ax.set_title("C: finer dx -> sharper peaks,\nsimilar domain-averaged energy (t=30)")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

    fig.suptitle("SD7003: diagnosing the ngrid=3 divergence and the dx=0.01 peak-vorticity increase", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    outpath = SURF / "airfoils" / "SD7003" / "3-ngrid=3" / "instability_diagnostic.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")

    print("\nPanel A data:")
    print("  baseline (ngrid=1):", list(zip(base_t, base_rms)))
    print("  ngrid=3:            ", list(zip(ng3_t, ng3_rms)))
    print("\nPanel B data (ngrid=3 core vs edge):")
    print("  core:", list(zip(core_t, core_rms)))
    print("  edge:", list(zip(core_t, edge_rms)))
    print("\nPanel C data (dx=0.02 vs dx=0.01 at t=30):")
    print("  max|omega|:", dict(zip(dxs, max_vals)))
    print("  RMS|omega|:", dict(zip(dxs, rms_vals)))


if __name__ == "__main__":
    main()
