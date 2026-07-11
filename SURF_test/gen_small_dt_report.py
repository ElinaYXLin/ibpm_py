"""
gen_small_dt_report.py

Report for SURF_test/run_small_dt_validation.py's output: does py/ibpm.py,
now calling the real FFTW3 library with FFTW_EXHAUSTIVE planning (see
py/elliptic_solver_2d.py), still agree with C++ build/ibpm -- both at a
much smaller dt (0.001, vs. this suite's usual 0.01) than used elsewhere?

Two kinds of evidence, both written to SURF_test/{SD7003,SD8000}/3-small_dt/:

  1. Force-trace comparison (small_dt_force_comparison.png) -- Cd/Cl(t) for
     both implementations overlaid, plus their pointwise difference on a log
     scale, the same "bit-identical early, chaotic divergence later" check
     used throughout this test suite (see SD8000/2-c++included/
     port_fidelity_diagnostic.png), now at 10x finer time resolution.

  2. Intermediate FIELD-value interception (small_dt_field_comparison.png +
     small_dt_field_diff_summary.txt) -- unlike the force-trace check (which
     only compares two scalar numbers per step), this loads the actual
     restart snapshots (vorticity + flux fields) both implementations wrote
     every RESTART_EVERY steps and computes the pointwise max/RMS difference
     directly on the full field, at every snapshot -- catching a fidelity
     regression the integrated force alone might average out or miss.

Usage:  python3 SURF_test/gen_small_dt_report.py
"""
import json
import pathlib
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

META = json.loads((REPO / "SURF_test" / "small_dt_run_meta.json").read_text())
DT = META["dt"]
RESTART_EVERY = META["restart_every"]
NSTEPS = META["nsteps"]


def load_force(p):
    d = np.loadtxt(p)
    return d[:, 1], d[:, 2], d[:, 3]  # t, Cd, Cl


def load_omega_q(path):
    s = State(filename=str(path))
    return s.omega._data[0].copy(), s.q._data[0].copy()


def main():
    for name in META["cases"]:
        outdir = REPO / "SURF_test" / "airfoils" / name / "3-small_dt"
        outdir.mkdir(parents=True, exist_ok=True)
        py_dir = REPO / "SURF_test" / "airfoils" / name / "_run_data_smalldt"
        cpp_dir = REPO / "SURF_test" / "airfoils" / name / "_run_data_smalldt_cpp"

        # ---------------- 1. force trace comparison ----------------
        tp, cdp, clp = load_force(py_dir / "run.force")
        tc, cdc, clc = load_force(cpp_dir / "run.force")
        n = min(len(tp), len(tc))
        dcd = np.abs(cdp[:n] - cdc[:n])
        dcl = np.abs(clp[:n] - clc[:n])
        onset_cl = int(np.argmax(dcl > 1e-6)) if np.any(dcl > 1e-6) else n - 1

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
        axes[0].plot(tp[:n], cdp[:n], color="C0", lw=1.1, label="py/ibpm.py (native FFTW3)")
        axes[0].plot(tc[:n], cdc[:n], color="C3", lw=1.0, ls="--", label="C++ build/ibpm")
        axes[0].set_xlabel("time t"); axes[0].set_ylabel("$C_d$")
        axes[0].set_title(f"{name}: $C_d(t)$, dt={DT}"); axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

        axes[1].plot(tp[:n], clp[:n], color="C0", lw=1.1, label="py/ibpm.py (native FFTW3)")
        axes[1].plot(tc[:n], clc[:n], color="C3", lw=1.0, ls="--", label="C++ build/ibpm")
        axes[1].set_xlabel("time t"); axes[1].set_ylabel("$C_l$")
        axes[1].set_title(f"{name}: $C_l(t)$, dt={DT}"); axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3)

        axes[2].semilogy(tp[:n], np.maximum(dcd, 1e-18), color="C1", lw=0.9, label="$|\\Delta C_d|$")
        axes[2].semilogy(tp[:n], np.maximum(dcl, 1e-18), color="#6a1b9a", lw=0.9, label="$|\\Delta C_l|$")
        # only annotate a real threshold crossing; if the diff never gets that
        # large within this run, say so instead (placed within the actual
        # plotted range, not at a hardcoded y that may fall miles outside a
        # near-machine-precision-only y-axis and blow up the saved bbox)
        crossed = bool(np.any(dcl > 1e-6))
        ymax_data = max(dcd.max(), dcl.max(), 1e-18)
        ytext = ymax_data * 3 if ymax_data > 0 else 1e-16
        if crossed:
            axes[2].axvline(tp[onset_cl], color="gray", ls="--", lw=1)
            axes[2].text(tp[onset_cl] + 0.02, ytext,
                         f"|ΔCl| crosses 1e-6\nat t={tp[onset_cl]:.3f} (step {onset_cl})",
                         fontsize=7.5, color="#6a1b9a", va="bottom")
        else:
            axes[2].text(tp[n // 2], ytext,
                         f"|ΔCl| never exceeds 1e-6 anywhere in this run\n"
                         f"(max = {dcl.max():.1e}) -- machine-precision agreement\n"
                         f"throughout, at this dt/step count",
                         fontsize=7.5, color="#6a1b9a", ha="center", va="bottom")
        axes[2].set_xlabel("time t"); axes[2].set_ylabel("|py - cpp|")
        axes[2].set_title("Pointwise difference (log)"); axes[2].legend(fontsize=8)
        axes[2].grid(alpha=0.3, which="both")

        fig.suptitle(f"{name}: py/ibpm.py (native FFTW3, FFTW_EXHAUSTIVE) vs. C++ build/ibpm, "
                     f"dt={DT} (10x smaller than this suite's usual dt=0.01)", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fig.savefig(outdir / "small_dt_force_comparison.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # ---------------- 2. intermediate field-value interception ----------------
        steps = list(range(0, NSTEPS + 1, RESTART_EVERY))
        rows = []
        avail_steps = []
        for s in steps:
            py_f = py_dir / f"run{s:05d}.bin"
            cpp_f = cpp_dir / f"run{s:05d}.bin"
            if not (py_f.exists() and cpp_f.exists()):
                continue
            omega_py, q_py = load_omega_q(py_f)
            omega_cpp, q_cpp = load_omega_q(cpp_f)
            d_omega = np.abs(omega_py - omega_cpp)
            d_q = np.abs(q_py - q_cpp)
            rows.append(dict(step=s, t=s * DT,
                              omega_max=float(np.abs(omega_py).max()),
                              domega_max=float(d_omega.max()), domega_rms=float(np.sqrt((d_omega**2).mean())),
                              dq_max=float(d_q.max()), dq_rms=float(np.sqrt((d_q**2).mean()))))
            avail_steps.append(s)

        if rows:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            t_arr = [r["t"] for r in rows]
            axes[0].semilogy(t_arr, [max(r["domega_max"], 1e-18) for r in rows],
                              "o-", color="C0", label="max |Δω|")
            axes[0].semilogy(t_arr, [max(r["domega_rms"], 1e-18) for r in rows],
                              "s--", color="C1", label="RMS |Δω|")
            axes[0].set_xlabel("time t"); axes[0].set_ylabel("|ω_py - ω_cpp|")
            axes[0].set_title(f"{name}: vorticity field, py vs. cpp\n(from restart snapshots every "
                              f"{RESTART_EVERY} steps)", fontsize=10)
            axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3, which="both")

            axes[1].semilogy(t_arr, [max(r["dq_max"], 1e-18) for r in rows],
                              "o-", color="C0", label="max |Δq|")
            axes[1].semilogy(t_arr, [max(r["dq_rms"], 1e-18) for r in rows],
                              "s--", color="C1", label="RMS |Δq|")
            axes[1].set_xlabel("time t"); axes[1].set_ylabel("|q_py - q_cpp|")
            axes[1].set_title(f"{name}: flux field, py vs. cpp", fontsize=10)
            axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3, which="both")

            fig.suptitle(f"{name}: intermediate-value interception (full vorticity/flux field, "
                         f"not just integrated force), dt={DT}", fontsize=11)
            fig.tight_layout(rect=(0, 0, 1, 0.92))
            fig.savefig(outdir / "small_dt_field_comparison.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

            with open(outdir / "small_dt_field_diff_summary.txt", "w") as f:
                f.write(f"{name}: intermediate field-value comparison, py/ibpm.py (native FFTW3) "
                        f"vs. C++ build/ibpm\n")
                f.write(f"dx={META['dx']} nx={META['nx']} ny={META['ny']} dt={DT} nsteps={NSTEPS} "
                        f"restart_every={RESTART_EVERY}\n\n")
                f.write(f"{'step':>6} {'t':>7} {'max|omega|':>11} {'max|domega|':>12} "
                        f"{'rms|domega|':>12} {'max|dq|':>10} {'rms|dq|':>10}\n")
                for r in rows:
                    f.write(f"{r['step']:6d} {r['t']:7.3f} {r['omega_max']:11.4f} "
                            f"{r['domega_max']:12.3e} {r['domega_rms']:12.3e} "
                            f"{r['dq_max']:10.3e} {r['dq_rms']:10.3e}\n")
            print(f"{name}: field snapshots compared at steps {avail_steps}")
        else:
            print(f"{name}: WARNING no matching restart snapshots found for field comparison")

        print(f"{name}: wrote figures to {outdir}")

    print("done")


if __name__ == "__main__":
    main()
