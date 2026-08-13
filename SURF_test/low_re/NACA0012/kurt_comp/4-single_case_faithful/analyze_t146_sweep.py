"""
analyze_t146_sweep.py

Reduces runs/t146_sweep/ (this folder's t=146 steady angle sweep, see
run_t146_sweep.py) into mean Cl/Cd per angle+impl, and plots it against
both Kurtulus's own digitized Fig 1 steady curve and stage 1's original
t=30 steady sweep -- isolating the effect of averaging-window duration
alone, since domain/dx/dt/angle-set are otherwise identical to stage 1.

Usage: python3 analyze_t146_sweep.py
Output: data/mean_coefficients_t146.csv, figures/fig1_mean_coefficients_t146.png
"""
import io
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
STAGE1 = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "1-paper_based"
STAGE4 = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp" / "4-single_case_faithful"
RUNS = STAGE4 / "runs" / "t146_sweep"
DATA = STAGE4 / "data"
FIGS = STAGE4 / "figures"
DATA.mkdir(exist_ok=True)
FIGS.mkdir(exist_ok=True)

AVG_FRAC = 0.5  # last 50% of the run, same convention as stage 1
IMPL_COLOR = {"py": "#1f77b4", "cpp": "#d62728"}
IMPL_LABEL = {"py": "py_static", "cpp": "cpp_static"}


def read_csv_with_comments(path):
    lines = [l for l in pathlib.Path(path).read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    return np.genfromtxt(io.StringIO("\n".join(lines)), delimiter=",",
                         names=True, dtype=None, encoding="utf-8")


def load_force(run_dir):
    f = run_dir / "flow.force"
    if not f.exists():
        return None
    try:
        d = np.loadtxt(f)
    except Exception:
        return None
    if d.ndim != 2 or len(d) < 10:
        return None
    return d  # cols: step, t, Cd, Cl


def mean_coeffs(d):
    i0 = int(len(d) * (1 - AVG_FRAC))
    seg = d[i0:]
    cd, cl = seg[:, 2], seg[:, 3]
    return dict(cl_mean=float(np.mean(cl)), cl_std=float(np.std(cl)),
                cd_mean=float(np.mean(cd)), cd_std=float(np.std(cd)),
                nan=bool(np.isnan(seg).any()))


def collect():
    rows = []
    for run in sorted(RUNS.iterdir()):
        if not run.is_dir():
            continue
        try:
            _, impl, atag = run.name.split("_")
            angle = int(atag[1:])
        except ValueError:
            continue
        d = load_force(run)
        if d is None:
            continue
        mc = mean_coeffs(d)
        rows.append((impl, angle, mc["cl_mean"], mc["cl_std"],
                     mc["cd_mean"], mc["cd_std"], int(mc["nan"])))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(DATA / "mean_coefficients_t146.csv", "w") as f:
        f.write("impl,alpha_deg,cl_mean,cl_std,cd_mean,cd_std,nan\n")
        for r in rows:
            f.write("%s,%d,%.5f,%.5f,%.5f,%.5f,%d\n" % r)
    print(f"wrote mean_coefficients_t146.csv ({len(rows)} rows)")
    return rows


def sub(mc, impl):
    m = (mc["impl"] == impl) & (mc["nan"] == 0)
    a = mc["alpha_deg"][m]
    o = np.argsort(a)
    return a[o], mc["cl_mean"][m][o], mc["cd_mean"][m][o]


def plot():
    mc_t146 = read_csv_with_comments(DATA / "mean_coefficients_t146.csv")
    mc_t30 = read_csv_with_comments(STAGE1 / "data" / "mean_coefficients_dx0.020.csv")
    k = read_csv_with_comments(STAGE1 / "data" / "kurtulus_fig1_digitized.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, qty, kcol, ylab in [
            (axes[0], "cl", "cl_steady", r"$\overline{C_l}$"),
            (axes[1], "cd", "cd_steady", r"$\overline{C_d}$")]:
        ax.plot(k["alpha_deg"], k[kcol], "k-", marker="^", mfc="white", mec="k",
                ms=5, lw=0.8, label="Kurtulus (2019) Fig 1, steady")
        for impl in ("py", "cpp"):
            m = (mc_t30["motion"] == "steady") & (mc_t30["impl"] == impl) & (mc_t30["nan"] == 0)
            a30 = mc_t30["alpha_deg"][m]; o = np.argsort(a30)
            a30 = a30[o]
            y30 = (mc_t30["cl_mean"] if qty == "cl" else mc_t30["cd_mean"])[m][o]
            ax.plot(a30, y30, "--", color=IMPL_COLOR[impl], lw=1.0, alpha=0.5,
                    marker="x", ms=3, label=f"IBPM t=30 (stage1) {IMPL_LABEL[impl]}")
            a146, cl146, cd146 = sub(mc_t146, impl)
            y146 = cl146 if qty == "cl" else cd146
            ax.plot(a146, y146, "-", color=IMPL_COLOR[impl], lw=1.5,
                    marker=".", ms=6, label=f"IBPM t=146 (this run) {IMPL_LABEL[impl]}")
        ax.set_xlabel(r"$\alpha_0$ [deg]"); ax.set_ylabel(ylab)
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle("Steady mean coefficients: does matching Kurtulus's t=146 "
                 "averaging duration change the sweep's conclusions?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIGS / "fig1_mean_coefficients_t146.png", dpi=140)
    plt.close(fig)
    print("wrote fig1_mean_coefficients_t146.png")


if __name__ == "__main__":
    collect()
    plot()
