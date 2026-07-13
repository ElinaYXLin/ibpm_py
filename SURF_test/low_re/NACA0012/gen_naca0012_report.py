"""
gen_naca0012_report.py

Report figures for the NACA0012 low-Reynolds-number validation against
published NON-LSAT CFD-benchmark drag coefficients. Recomputes mean
Cl/Cd directly from each run's .force file (robust to the results-JSON
concurrent-write race), for both py/ibpm.py and C++ build/ibpm.

Reference drag (alpha=0, Cl=0 by NACA0012 symmetry):
  Re=500:  Cd = 0.1762 (Lockard et al.), 0.1759 (Wu et al.), 0.178 (Nita et al. LBM)
  Re=1000: Cd = 0.119  (Di Ilio et al. HLBM & XFOIL), ~0.12 (Kurtulus)

Outputs:
  polar_comparison.png   -- Cl(alpha), Cd(alpha) at Re=500, py vs cpp,
                            with the Re=500 reference Cd(0) band marked
  grid_convergence.png   -- Cd(alpha=0, Re=500) vs dx, py vs cpp, vs the
                            reference band (if the grid-convergence runs exist)
  fidelity_summary.txt   -- py-vs-cpp Cl/Cd agreement table + benchmark comparison

Usage: python3 SURF_test/low_re/NACA0012/gen_naca0012_report.py
"""
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
OUTDIR = REPO / "SURF_test" / "low_re" / "NACA0012"
AVG_FRAC = 0.6

# Reference Cd at alpha=0 (Cl=0 by symmetry), non-LSAT CFD benchmarks
REF_CD = {
    500: {"Lockard et al.": 0.1762, "Wu et al.": 0.1759, "Nita et al. (LBM)": 0.178},
    1000: {"Di Ilio et al. (HLBM)": 0.119, "Di Ilio et al. (XFOIL)": 0.119, "Kurtulus (~)": 0.12},
}
POLAR_ALPHAS = [0, 2, 4, 6, 8, 10]
GRID_DX = [0.04, 0.02, 0.01]


def mean_force(force_path, frac=AVG_FRAC):
    d = np.loadtxt(force_path)
    if d.ndim == 1:
        d = d[None, :]
    seg = d[int(len(d) * (1 - frac)):]
    return float(seg[:, 3].mean()), float(seg[:, 2].mean())  # cl, cd


def load_polar(subdir):
    """Return {(Re,alpha): (cl,cd)} for whichever runs exist under subdir."""
    out = {}
    base = OUTDIR / subdir
    for d in sorted(base.glob("Re*_a*")):
        fp = d / "run.force"
        if not fp.exists():
            continue
        # parse Re<val>_a<±NN>
        name = d.name
        Re = int(name.split("_a")[0][2:])
        alpha = int(name.split("_a")[1])
        out[(Re, alpha)] = mean_force(fp)
    return out


def load_gridconv(subdir):
    """Return {dx: (cl,cd)} for Re=500 alpha=0 grid-convergence runs if present."""
    out = {}
    base = OUTDIR / subdir
    for dx in GRID_DX:
        fp = base / f"gridconv_dx{dx}" / "run.force"
        if fp.exists():
            out[dx] = mean_force(fp)
    return out


def main():
    py = load_polar("_run_data")
    cpp = load_polar("_run_data_cpp")
    if not py or not cpp:
        print("no polar runs found yet")
        return

    # ---------- polar figure (Re=500) ----------
    a = [al for al in POLAR_ALPHAS if (500, al) in py and (500, al) in cpp]
    py_cl = [py[(500, al)][0] for al in a]; py_cd = [py[(500, al)][1] for al in a]
    cpp_cl = [cpp[(500, al)][0] for al in a]; cpp_cd = [cpp[(500, al)][1] for al in a]

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(a, py_cl, "o-", color="C0", label="py/ibpm.py")
    ax[0].plot(a, cpp_cl, "^--", color="C3", label="C++ build/ibpm")
    ax[0].set_xlabel(r"$\alpha$ (deg)"); ax[0].set_ylabel("$C_l$")
    ax[0].set_title("NACA0012 Re=500: lift coefficient"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    ax[1].plot(a, py_cd, "o-", color="C0", label="py/ibpm.py")
    ax[1].plot(a, cpp_cd, "^--", color="C3", label="C++ build/ibpm")
    refs = REF_CD[500]
    lo, hi = min(refs.values()), max(refs.values())
    ax[1].axhspan(lo, hi, color="0.7", alpha=0.5, zorder=0,
                  label=f"ref. Cd($\\alpha$=0): {lo:.4f}-{hi:.4f}\n(Lockard/Wu/Nita)")
    ax[1].set_xlabel(r"$\alpha$ (deg)"); ax[1].set_ylabel("$C_d$")
    ax[1].set_title("NACA0012 Re=500: drag coefficient"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.suptitle("NACA0012 at Re=500: py/ibpm.py vs. C++ build/ibpm vs. published CFD-benchmark drag\n"
                 "(non-LSAT reference; no wind-tunnel data exists at Re in the hundreds -- see README)")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUTDIR / "polar_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 'polar_comparison.png'}")

    # ---------- grid-convergence figure (Re=500, alpha=0) ----------
    gc_py = load_gridconv("_run_data")
    gc_cpp = load_gridconv("_run_data_cpp")
    if gc_py and gc_cpp:
        dxs = sorted(set(gc_py) & set(gc_cpp), reverse=True)
        fig, axg = plt.subplots(figsize=(6.5, 5))
        axg.plot(dxs, [gc_py[d][1] for d in dxs], "o-", color="C0", label="py/ibpm.py")
        axg.plot(dxs, [gc_cpp[d][1] for d in dxs], "^--", color="C3", label="C++ build/ibpm")
        axg.axhspan(lo, hi, color="0.7", alpha=0.5, zorder=0, label="ref. Cd (Lockard/Wu/Nita)")
        axg.invert_xaxis(); axg.set_xlabel("grid spacing dx"); axg.set_ylabel("$C_d$ ($\\alpha$=0)")
        axg.set_title("NACA0012 Re=500, $\\alpha$=0: drag grid convergence")
        axg.legend(fontsize=9); axg.grid(alpha=0.3)
        fig.savefig(OUTDIR / "grid_convergence.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {OUTDIR / 'grid_convergence.png'}")

    # ---------- fidelity + benchmark summary ----------
    with open(OUTDIR / "fidelity_summary.txt", "w") as f:
        f.write("NACA0012 low-Re validation: py/ibpm.py vs. C++ build/ibpm, and vs. CFD benchmark\n\n")
        f.write("Polar (dx=0.02, nx=300, ny=150, t=30, mean over last 60%):\n")
        f.write(f"{'Re':>6} {'alpha':>6} {'Cl_py':>10} {'Cl_cpp':>10} {'Cd_py':>10} {'Cd_cpp':>10} {'|dCd|/Cd_cpp':>14}\n")
        for key in sorted(set(py) & set(cpp)):
            Re, al = key
            clp, cdp = py[key]; clc, cdc = cpp[key]
            rel = abs(cdp - cdc) / max(abs(cdc), 1e-12)
            f.write(f"{Re:6d} {al:6d} {clp:10.4f} {clc:10.4f} {cdp:10.4f} {cdc:10.4f} {rel:14.2e}\n")
        f.write("\nDrag benchmark comparison at alpha=0 (Cl=0 by symmetry):\n")
        for Re in (500, 1000):
            if (Re, 0) in py:
                cd_py = py[(Re, 0)][1]; cd_cpp = cpp[(Re, 0)][1]
                f.write(f"  Re={Re}: Cd_py={cd_py:.4f}  Cd_cpp={cd_cpp:.4f}\n")
                for ref, val in REF_CD[Re].items():
                    f.write(f"      vs {ref:22s} {val:.4f}  ({(cd_py-val)/val:+.1%})\n")
    print(f"wrote {OUTDIR / 'fidelity_summary.txt'}")
    print(open(OUTDIR / "fidelity_summary.txt").read())


if __name__ == "__main__":
    main()
