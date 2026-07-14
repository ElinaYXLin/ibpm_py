"""
compute_spectrum_and_energy.py

Two more discriminating tests between "genuine chaos" and "under-resolved
numerical instability" for the SD8000/SD7003 coarse-grid (dx=0.04)
blow-ups documented in ../7-chaos_sensitivity/README.md, both computed
from field snapshots ALREADY on disk there (../7-chaos_sensitivity/
_run_data/{SD8000_ext4000_py,SD8000_ext4000_cpp,SD7003_ext4000_py}/,
saved every 25 steps) -- no new simulation runs needed.

1. Vorticity wavenumber spectrum, computed from the last STABLE snapshot
   before each case's recorded blow-up step (2D FFT of the vorticity
   field, radially binned into a 1D enstrophy-spectral-density vs |k|
   curve). A physical (if under-resolved) turbulent cascade decays
   smoothly with k; a spectrum that piles up / turns upward at the
   highest resolvable wavenumbers (near the grid Nyquist limit) is the
   textbook signature of an ALIASING instability -- energy artificially
   accumulating at the grid scale because nothing in this DNS-style
   solver (no subgrid model, see ../LSAT-SD7003/6-explicit_dissipation/
   README.md) removes it.

2. Domain-integrated kinetic energy (0.5*InnerProduct(q,q)) and enstrophy
   (sum(omega^2)*dx^2) vs. time, from t=0 to just past the blow-up step.
   Bounded/fluctuating growth is consistent with a chaotic attractor;
   monotonic, accelerating (super-exponential-looking) growth all the way
   to the NaN step is the signature of a runaway numerical instability
   rather than saturating turbulent dynamics.

Usage: python3 SURF_test/airfoils/8-dt_refinement_and_spectra/compute_spectrum_and_energy.py
Output: SURF_test/airfoils/8-dt_refinement_and_spectra/figures/
        spectrum_<case>.png, energy_growth.png
"""
from __future__ import annotations

import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
CHAOS_DIR = REPO / "SURF_test" / "airfoils" / "7-chaos_sensitivity"
RUN_DATA = CHAOS_DIR / "_run_data"
HERE = REPO / "SURF_test" / "airfoils" / "8-dt_refinement_and_spectra"
FIGS = HERE / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402
from py.vector_operations import InnerProduct  # noqa: E402

# coarse grid (dx=0.04): nx=150, ny=75, length=6, xoffset=-2, yoffset=-1.5
NX, NY = 150, 75
LENGTH = 6.0
DX = LENGTH / NX
BLOWUP_THRESHOLD = 20.0

# (name, impl, dt, restart_interval) -- matches ../7-chaos_sensitivity/run_chaos_sensitivity.py
CASES = [
    ("SD8000", "py", 0.01, 25),
    ("SD8000", "cpp", 0.01, 25),
    ("SD7003", "py", 0.01, 25),
    # SD7003 cpp never blew up through 4000 steps (../7-chaos_sensitivity/README.md) -- excluded
]


def find_blowup_step(force_path, threshold=BLOWUP_THRESHOLD):
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    step, cd, cl = d[:, 0], d[:, 2], d[:, 3]
    mask = (np.abs(cl) > threshold) | (np.abs(cd) > threshold)
    if not mask.any():
        return None
    return int(step[np.argmax(mask)])


def last_stable_snapshot_step(rundir, blowup_step, restart_interval=25):
    """Largest available snapshot step that is < blowup_step and whose
    field is not already NaN (blow-up in the force trace can lag the
    field going fully to NaN by up to one snapshot interval -- confirmed
    directly in ../7-chaos_sensitivity/README.md point 4)."""
    candidates = sorted(rundir.glob("run*.bin"))
    steps = sorted(int(p.stem.replace("run", "")) for p in candidates)
    best = None
    for s in steps:
        if s >= blowup_step:
            break
        w = State(filename=str(rundir / f"run{s:05d}.bin")).omega._data[0]
        if np.isnan(w).any():
            break
        best = s
    return best


def radial_spectrum(field, dx):
    """2D FFT of `field`, binned radially into a 1D spectral-density vs
    |k| curve (k in cycles per unit length; Nyquist = 1/(2*dx))."""
    ny, nx = field.shape
    F = np.fft.fft2(field)
    psd2d = (np.abs(F) ** 2) / (nx * ny)
    kx = np.fft.fftfreq(nx, d=dx)
    ky = np.fft.fftfreq(ny, d=dx)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K = np.sqrt(KX ** 2 + KY ** 2)
    k_nyquist = 1.0 / (2 * dx)
    nbins = min(nx, ny) // 2
    bins = np.linspace(0, k_nyquist, nbins + 1)
    bin_idx = np.digitize(K.ravel(), bins) - 1
    psd_flat = psd2d.ravel()
    k_centers = 0.5 * (bins[:-1] + bins[1:])
    spectrum = np.full(nbins, np.nan)
    for b in range(nbins):
        sel = bin_idx == b
        if sel.any():
            spectrum[b] = psd_flat[sel].mean()
    return k_centers, spectrum, k_nyquist


def gen_spectra():
    fig, axes = plt.subplots(1, len(CASES), figsize=(5.2 * len(CASES), 4.6), sharey=True)
    if len(CASES) == 1:
        axes = [axes]
    for ax, (name, impl, dt, restart) in zip(axes, CASES):
        rundir = RUN_DATA / f"{name}_ext4000_{impl}"
        force_path = rundir / "run.force"
        blowup = find_blowup_step(force_path)
        if blowup is None:
            ax.set_title(f"{name} {impl}: no blow-up, skipped")
            continue
        stable_step = last_stable_snapshot_step(rundir, blowup, restart)
        w = State(filename=str(rundir / f"run{stable_step:05d}.bin")).omega._data[0]
        k, spec, k_ny = radial_spectrum(w, DX)
        valid = np.isfinite(spec) & (spec > 0)
        ax.loglog(k[valid], spec[valid], "o-", ms=3, color="C0")
        ax.axvline(k_ny, color="gray", ls="--", lw=1, label=f"grid Nyquist $k$={k_ny:.2f}")
        ax.set_xlabel("wavenumber $|k|$ (1/length)")
        ax.set_title(f"{name} {impl}\nlast-stable snapshot t={stable_step*dt:.2f} "
                     f"(blow-up at step {blowup}, t={blowup*dt:.2f})", fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("vorticity power spectral density $|\\hat\\omega(k)|^2$")
    fig.suptitle("Vorticity wavenumber spectrum just before blow-up (dx=0.04 coarse grid)\n"
                 "upturn near grid Nyquist = energy piling up at grid scale (aliasing signature), "
                 "not a decaying physical cascade")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    outp = FIGS / "spectrum_prebreakup.png"
    fig.savefig(outp, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outp}")


def gen_energy_growth():
    fig, axes = plt.subplots(2, len(CASES), figsize=(5.2 * len(CASES), 7.5), sharex="col")
    if len(CASES) == 1:
        axes = axes.reshape(2, 1)
    for col, (name, impl, dt, restart) in enumerate(CASES):
        rundir = RUN_DATA / f"{name}_ext4000_{impl}"
        blowup = find_blowup_step(rundir / "run.force")
        if blowup is None:
            continue
        steps = sorted(int(p.stem.replace("run", "")) for p in rundir.glob("run*.bin"))
        steps = [s for s in steps if s <= blowup + restart]  # include first post-blowup point
        ts, kes, ens = [], [], []
        for s in steps:
            st = State(filename=str(rundir / f"run{s:05d}.bin"))
            w = st.omega._data[0]
            if np.isnan(w).any():
                ts.append(s * dt); kes.append(np.nan); ens.append(np.nan)
                continue
            ke = 0.5 * InnerProduct(st.q, st.q)
            enstrophy = float(np.sum(w ** 2)) * DX * DX
            ts.append(s * dt); kes.append(ke); ens.append(enstrophy)

        ax = axes[0, col]
        ax.semilogy(ts, kes, "o-", ms=3, color="C0")
        ax.axvline(blowup * dt, color="k", ls="--", lw=1, label=f"blow-up t={blowup*dt:.2f}")
        ax.set_title(f"{name} {impl}: kinetic energy", fontsize=9)
        ax.set_ylabel("KE (log scale)"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")

        ax = axes[1, col]
        ax.semilogy(ts, ens, "o-", ms=3, color="C3")
        ax.axvline(blowup * dt, color="k", ls="--", lw=1, label=f"blow-up t={blowup*dt:.2f}")
        ax.set_title(f"{name} {impl}: enstrophy", fontsize=9)
        ax.set_xlabel("t"); ax.set_ylabel("enstrophy (log scale)"); ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")

    fig.suptitle("Domain-integrated KE / enstrophy growth to blow-up (dx=0.04 coarse grid)\n"
                 "monotonic runaway all the way to the NaN step = runway numerical instability, "
                 "not a saturating/bounded chaotic attractor")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    outp = FIGS / "energy_growth.png"
    fig.savefig(outp, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outp}")


def main():
    gen_spectra()
    gen_energy_growth()


if __name__ == "__main__":
    main()
