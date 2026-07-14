"""
gen_sd_gridconv_extended_figs.py

Extends ../LSAT-SD7003/README.md's / ../LSAT-SD8000/README.md's original
grid_convergence.png (py only, dx=0.04/0.02/0.01, from
run_all_airfoils.py's batch_results.json) with the finer dx=0.005 and
dx=0.0025 levels added by run_sd_gridconv_extended.py, and with the C++
(build/ibpm) points from run_all_airfoils_cpp.py / the same extended
sweep overlaid -- following the same convention as
SURF_test/low_re/NACA0012/gen_naca0012_report.py's grid_convergence.png.

Usage: python3 SURF_test/gen_sd_gridconv_extended_figs.py
Output: SURF_test/airfoils/LSAT-<name>/grid_convergence_extended.png
        SURF_test/airfoils/LSAT-<name>/grid_convergence_extended_summary.txt
"""
import json
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
SURF = REPO / "SURF_test"
import sys
sys.path.insert(0, str(SURF))
from parse_uiuc import parse_blocks, nearest_block  # noqa: E402

CASES = {
    "SD7003": dict(Re=61100, conv_alpha=-0.09),
    "SD8000": dict(Re=60800, conv_alpha=-0.81),
}
EXTRA_DX = [0.005, 0.0025]


def load_all_points(name):
    py = json.loads((SURF / "batch_results.json").read_text())
    cpp = json.loads((SURF / "batch_results_cpp.json").read_text())
    pts = {"py": [(c["dx"], c["cd_mean"], c["cd_std"]) for c in py["convergence"][name]],
           "cpp": [(c["dx"], c["cd_mean"], c["cd_std"]) for c in cpp["convergence"][name]]}
    for impl in ("py", "cpp"):
        for dx in EXTRA_DX:
            fp = SURF / f"batch_results_extended_{impl}_{name}.json"
            ext = json.loads(fp.read_text())
            s = ext[str(dx)]
            pts[impl].append((dx, s["cd_mean"], s["cd_std"]))
        pts[impl].sort(key=lambda t: -t[0])
    return pts


def main():
    for name, cfg in CASES.items():
        outdir = REPO / "SURF_test" / "airfoils" / f"LSAT-{name}"
        pts = load_all_points(name)

        drg_blocks = parse_blocks(SURF / "airfoils" / f"LSAT-{name}" / f"{name}.DRG.txt", "drg")
        exp = nearest_block(drg_blocks, cfg["Re"])
        exp_alpha = np.array(exp["alpha"]); exp_cd = np.array(exp["Cd"])
        exp_row = np.argmin(np.abs(exp_alpha - cfg["conv_alpha"]))
        exp_cd_at_alpha = exp_cd[exp_row]

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        for ax, skip_coarse in zip(axes, [False, True]):
            for impl, color, marker in [("py", "C0", "o"), ("cpp", "C3", "^")]:
                rows = [p for p in pts[impl] if not (skip_coarse and p[0] == 0.04)]
                dxs = [p[0] for p in rows]
                cds = [p[1] for p in rows]
                stds = [p[2] for p in rows]
                ax.errorbar(dxs, cds, yerr=stds, fmt=f"{marker}-", color=color, capsize=3,
                            label=f"{'py/ibpm.py' if impl == 'py' else 'C++ build/ibpm'}")
            ax.axhline(exp_cd_at_alpha, color="k", ls="--", lw=1.2,
                       label=f"UIUC LSAT exp. ($\\alpha$={exp_alpha[exp_row]:.2f}°, Re={exp['Re']:.0f})")
            ax.invert_xaxis()
            ax.set_xscale("log")
            ax.set_xlabel("grid spacing dx (log scale)")
            ax.set_ylabel("$C_d$")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3, which="both")
            ax.set_title("dx=0.04 excluded (zoomed; dx=0.04's huge error bars\n"
                         "are the numerical-instability blow-up, see ../8-dt_refinement_and_spectra/)"
                         if skip_coarse else "full range (dx=0.04 included)")
        fig.suptitle(f"{name}: extended grid convergence, $\\alpha$={cfg['conv_alpha']}°, "
                     f"Re={cfg['Re']}\n(dx=0.04 to 0.0025; error bars = ±1 std. dev. of unsteady force trace)")
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        outp = outdir / "grid_convergence_extended.png"
        fig.savefig(outp, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {outp}")

        # ---- summary text: dx, Cd, dCd, %-of-Cd, py-vs-cpp agreement ----
        summary_path = outdir / "grid_convergence_extended_summary.txt"
        with open(summary_path, "w") as f:
            f.write(f"{name}: extended grid convergence, alpha={cfg['conv_alpha']}, Re={cfg['Re']}\n")
            f.write(f"UIUC LSAT experimental Cd at alpha={exp_alpha[exp_row]:.2f}: {exp_cd_at_alpha:.5f}\n\n")
            for impl in ("py", "cpp"):
                f.write(f"[{impl}]\n")
                f.write(f"{'dx':>8} {'Cd':>10} {'dCd':>12} {'|dCd|/Cd':>10}\n")
                prev = None
                for dx, cd, std in pts[impl]:
                    if prev is None:
                        f.write(f"{dx:8.4f} {cd:10.6f} {'':>12} {'':>10}\n")
                    else:
                        dcd = cd - prev
                        f.write(f"{dx:8.4f} {cd:10.6f} {dcd:+12.6f} {abs(dcd)/abs(cd):9.1%}\n")
                    prev = cd
                f.write("\n")
            f.write("py-vs-cpp agreement at each dx:\n")
            for (dx_p, cd_p, _), (dx_c, cd_c, _) in zip(pts["py"], pts["cpp"]):
                assert dx_p == dx_c
                f.write(f"  dx={dx_p:.4f}: Cd_py={cd_p:.6f}  Cd_cpp={cd_c:.6f}  "
                        f"rel.diff={abs(cd_p-cd_c)/abs(cd_c):.2%}\n")
        print(f"wrote {summary_path}")
        print(open(summary_path).read())


if __name__ == "__main__":
    main()
