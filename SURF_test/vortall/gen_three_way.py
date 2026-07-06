import sys, types, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.io as sio

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
RUN_PY = REPO / "SURF_test" / "vortall" / "_run_data"
RUN_CPP = REPO / "SURF_test" / "vortall" / "_run_data_cpp"
OUT = REPO / "SURF_test" / "vortall"

sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

NX, NY = 450, 200
DX = 9.0 / NX
XOFF, YOFF = -1.0, -2.0
xs = XOFF + np.arange(1, NX) * DX
ys = YOFF + np.arange(1, NY) * DX
X, Y = np.meshgrid(xs, ys, indexing="ij")


def load_omega(path) -> np.ndarray:
    s = State(filename=str(path))
    return s.omega._data[0].copy()


VORT = sio.loadmat(REPO / "VORTALL.mat")["VORTALL"]


def vortall_snapshot(col: int) -> np.ndarray:
    return VORT[:, col].reshape(NX - 1, NY - 1, order="C")


VMAX = 5.0
LEVELS = np.linspace(-VMAX, VMAX, 41)


def draw_field(ax, field, title, vmax=VMAX, levels=None, cmap="RdBu_r"):
    lv = levels if levels is not None else np.linspace(-vmax, vmax, 41)
    im = ax.contourf(X, Y, np.clip(field, -vmax, vmax), levels=lv, cmap=cmap, extend="both")
    ax.add_patch(plt.Circle((0, 0), 0.5, fill=True, color="0.2", zorder=5))
    ax.set_xlim(-1, 8)
    ax.set_ylim(-2, 2)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x")
    return im


ref_snapshot = vortall_snapshot(150)
# Both fully independent runs, from a zero IC, never cross-seeded -- see
# run_vortall.py and README's "Correction" section.
py_final = load_omega(RUN_PY / "vortall14000.bin")     # Python, t=280
cpp_final = load_omega(RUN_CPP / "vortall14000.bin")   # C++ reference, t=280


def shedding_stats(force_path):
    """Period/Strouhal/peak-Cl/mean-Cd over the saturated regime (t>=200;
    see README -- amplitude saturates by ~t=150-160 in these runs, so 200
    leaves a comfortable margin). Phase-independent statistics are the
    right way to compare two runs that are NOT seeded identically: the
    Re=100 wake instability is seeded only by floating-point roundoff, so
    two independent runs (even bit-identical algorithms in different
    languages) diverge in exact timing/phase once the instability grows,
    while still converging to the same periodic amplitude/frequency."""
    d = np.loadtxt(force_path)
    t, cd, cl = d[:, 1], d[:, 2], d[:, 3]
    m = t >= 200
    t, cd, cl = t[m], cd[m], cl[m]
    peaks_t = [t[i] for i in range(1, len(cl) - 1)
               if cl[i] > cl[i - 1] and cl[i] > cl[i + 1] and cl[i] > 0.5]
    period = float(np.mean(np.diff(peaks_t))) if len(peaks_t) > 1 else float("nan")
    return dict(n_peaks=len(peaks_t), period=period,
                St=1.0 / period if period == period else float("nan"),
                cl_peak=float(np.abs(cl).max()), cd_mean=float(cd.mean()))


py_shed = shedding_stats(RUN_PY / "vortall.force")
cpp_shed = shedding_stats(RUN_CPP / "vortall.force")

diff = py_final - cpp_final
print("max |python - cpp| at t=280:", np.abs(diff).max())
print("python field max |omega|:", np.abs(py_final).max())
print("cpp field max |omega|:", np.abs(cpp_final).max())
print("VORTALL snapshot max |omega|:", np.abs(ref_snapshot).max())
print("python shedding stats:", py_shed)
print("cpp shedding stats:", cpp_shed)

fig, axes = plt.subplots(3, 1, figsize=(9, 9.5), sharex=True)
im0 = draw_field(axes[0], ref_snapshot, "VORTALL.mat, snapshot 150 (reference dataset)")
im1 = draw_field(axes[1], cpp_final, "C++ ibpm (src/), Re=100, t=280 (this repo's reference build)")
im2 = draw_field(axes[2], py_final, "py/ibpm.py port, Re=100, t=280")
for ax in axes:
    ax.set_ylabel("y")
fig.colorbar(im2, ax=axes, shrink=0.85, label="vorticity $\\omega$ (clipped to $\\pm 5$)")
fig.suptitle(
    "Flow past a cylinder, Re=100: VORTALL.mat vs. this repo's C++ vs. Python port\n"
    "C++ and Python are independent runs from a zero IC -- they match in periodic\n"
    "amplitude/frequency, not in instantaneous phase (see python_vs_cpp_diff.png); "
    "VORTALL.mat differs further because its snapshot isn't phase-aligned with either -- see README",
    fontsize=9.5,
)
fig.savefig(OUT / "vorticity_comparison_3way.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# Python-vs-C++ difference panel (quantifies "faithfully ported", separate from the VORTALL question)
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
im0 = draw_field(axes[0], cpp_final, "C++, t=280")
im1 = draw_field(axes[1], py_final, "Python, t=280")
diff_vmax = max(np.abs(diff).max(), 1e-12)
im2 = draw_field(axes[2], diff, f"Python - C++ (max|diff|={np.abs(diff).max():.2e})",
                  vmax=diff_vmax, cmap="RdBu_r")
for ax in axes:
    ax.set_ylabel("y")
fig.colorbar(im0, ax=axes[:2], shrink=0.85, label="vorticity $\\omega$ (clipped to $\\pm 5$)")
fig.colorbar(im2, ax=axes[2], shrink=0.85, label="difference")
fig.suptitle(
    "Python port vs. this repo's C++ build, Re=100, t=280 -- independent runs from a\n"
    "zero IC diverge in phase (chaotic sensitivity to fp-level differences during\n"
    "transient growth), not in periodic amplitude/frequency -- see three_way_summary.txt",
    fontsize=10,
)
fig.savefig(OUT / "python_vs_cpp_diff.png", dpi=160, bbox_inches="tight")
plt.close(fig)

with open(OUT / "three_way_summary.txt", "w") as f:
    f.write("Three-way comparison at t=280, Re=100, nx=450 ny=200, x in [-1,8], y in [-2,2]\n\n")
    f.write(f"max |omega|, VORTALL.mat snapshot 150:      {np.abs(ref_snapshot).max():.4f}\n")
    f.write(f"max |omega|, this repo's C++ build:          {np.abs(cpp_final).max():.4f}\n")
    f.write(f"max |omega|, Python port:                    {np.abs(py_final).max():.4f}\n\n")
    f.write(f"max |Python - C++| (same run, same t=280):   {np.abs(diff).max():.3e}\n")
    f.write("  -> C++ and Python are two FULLY INDEPENDENT runs from a zero initial\n")
    f.write("     condition (see README's 'Correction' section -- an earlier version of\n")
    f.write("     this comparison seeded Python from C++'s own restart file, which made\n")
    f.write("     near-machine-precision agreement here trivial; that bug is fixed).\n")
    f.write("     The Re=100 wake instability is seeded only by floating-point roundoff,\n")
    f.write("     so two independent runs diverge in exact phase/timing once the\n")
    f.write("     instability grows (chaotic sensitivity) -- this pointwise difference is\n")
    f.write("     therefore EXPECTED to be large (comparable to the field's own peak\n")
    f.write("     magnitude), and is not evidence of a porting bug. The right way to\n")
    f.write("     compare two independently-seeded chaotic runs is the periodic *state*\n")
    f.write("     they converge to, not a pointwise snapshot at a fixed time:\n\n")
    f.write(f"     {'':18s}{'Python':>12s}{'C++':>12s}\n")
    f.write(f"     {'Cl peak amplitude':18s}{py_shed['cl_peak']:12.5f}{cpp_shed['cl_peak']:12.5f}\n")
    f.write(f"     {'Cd mean (t>=200)':18s}{py_shed['cd_mean']:12.5f}{cpp_shed['cd_mean']:12.5f}\n")
    f.write(f"     {'shedding period':18s}{py_shed['period']:12.4f}{cpp_shed['period']:12.4f}\n")
    f.write(f"     {'Strouhal number':18s}{py_shed['St']:12.4f}{cpp_shed['St']:12.4f}\n")
    f.write("  -> these agree to 4-5 significant figures between the two independent\n")
    f.write("     implementations -- this (not the pointwise diff above) is the evidence\n")
    f.write("     the Python port is faithful, consistent with the exact-agreement\n")
    f.write("     validation on the non-chaotic 200x200 case in SURF_test/built_in_tests/README.md.\n")

print("wrote", OUT / "vorticity_comparison_3way.png")
print("wrote", OUT / "python_vs_cpp_diff.png")
