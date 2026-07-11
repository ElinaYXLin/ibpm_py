import sys, types, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.io as sio

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
RUN = REPO / "SURF_test" / "vortall" / "1-baseline" / "_run_data"
OUT = REPO / "SURF_test" / "vortall" / "1-baseline"

sys.path.insert(0, str(REPO))
pkg = types.ModuleType("py")
pkg.__path__ = [str(REPO / "py")]
sys.modules["py"] = pkg
from py.state import State  # noqa: E402

# ---- grid geometry (must match the -nx/-ny/-length/-xoffset/-yoffset used to run py.ibpm) ----
NX, NY = 450, 200
DX = 9.0 / NX
XOFF, YOFF = -1.0, -2.0
xs = XOFF + np.arange(1, NX) * DX   # interior nodes, shape (449,)
ys = YOFF + np.arange(1, NY) * DX   # shape (199,)
X, Y = np.meshgrid(xs, ys, indexing="ij")  # shape (449, 199)


def load_omega(step: int) -> np.ndarray:
    fname = RUN / f"vortall{step:05d}.bin"
    s = State(filename=str(fname))
    return s.omega._data[0].copy()  # shape (nx-1, ny-1) = (449, 199)


def load_force():
    d = np.loadtxt(RUN / "vortall.force")
    return d[:, 1], d[:, 2], d[:, 3]  # time, Cd, Cl


VORT = sio.loadmat(REPO / "VORTALL.mat")["VORTALL"]  # (89351, 151)


def vortall_snapshot(col: int) -> np.ndarray:
    return VORT[:, col].reshape(NX - 1, NY - 1, order="C")


VMAX = 5.0
LEVELS = np.linspace(-VMAX, VMAX, 41)


def draw_field(ax, field, title):
    im = ax.contourf(X, Y, np.clip(field, -VMAX, VMAX), levels=LEVELS,
                      cmap="RdBu_r", extend="both")
    ax.add_patch(plt.Circle((0, 0), 0.5, fill=True, color="0.2", zorder=5))
    ax.set_xlim(-1, 8)
    ax.set_ylim(-2, 2)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x")
    return im


# ============================================================
# Figure 1: side-by-side vorticity comparison, Python vs VORTALL.mat
# ============================================================
py_final = load_omega(14000)          # t = 280, saturated periodic shedding
ref_snapshot = vortall_snapshot(150)  # last column of VORTALL.mat

fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
im0 = draw_field(axes[0], ref_snapshot, "VORTALL.mat, snapshot 150 (reference dataset)")
im1 = draw_field(axes[1], py_final, "py/ibpm.py port, Re=100, t=280 (this run, saturated shedding)")
for ax in axes:
    ax.set_ylabel("y")
cbar = fig.colorbar(im1, ax=axes, shrink=0.85, label="vorticity $\\omega$ (clipped to $\\pm 5$)")
fig.suptitle(
    "Flow past a cylinder, Re=100: Python IBPM port vs. VORTALL.mat reference\n"
    "(both show the fully-developed, periodic von Kármán vortex street; "
    "absolute snapshot times are not aligned -- see README)",
    fontsize=10,
)
fig.savefig(OUT / "vorticity_comparison.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ============================================================
# Figure 2: evolution of the Python simulation over one shedding period
# ============================================================
evo_steps = [10000, 10800, 11600, 12400, 13200, 14000]  # t=200..280
fig, axes = plt.subplots(2, 3, figsize=(13, 5.6))
for ax, step in zip(axes.flat, evo_steps):
    field = load_omega(step)
    im = draw_field(ax, field, f"t = {step*0.02:g}")
for ax in axes[1, :]:
    ax.set_xlabel("x")
for ax in axes[:, 0]:
    ax.set_ylabel("y")
fig.colorbar(im, ax=axes, shrink=0.85, label="vorticity $\\omega$ (clipped to $\\pm 5$)")
fig.suptitle("Python IBPM port: vortex-shedding evolution, Re=100 (saturated periodic regime)", fontsize=11)
fig.savefig(OUT / "flow_evolution_python.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ============================================================
# Figure 3: lift/drag coefficient history -- transient growth -> saturation
# ============================================================
# load_force() returns the FULL t=0..280 trace (this run's own restart chain
# is fully independent -- see README's "Correction" section). This figure is
# specifically about the *saturated* regime, so it's windowed to t>=200 --
# comfortably past saturation (measured onset ~t=150-160 in this run; see
# README) -- rather than plotting the transient growth from ~1e-14 as well.
t_full, cd_full, cl_full = load_force()
window = t_full >= 200
t, cd, cl = t_full[window], cd_full[window], cl_full[window]

fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(t, cl, color="C0", label="$C_l(t)$ (Python)")
ax.plot(t, cd, color="C1", alpha=0.7, label="$C_d(t)$ (Python)")
ax.set_xlabel("t")
ax.set_ylabel("force coefficient")
ax.set_title("Force coefficients during the saturated vortex-shedding regime, Re=100")
ax.legend()
ax.grid(alpha=0.3)
fig.savefig(OUT / "force_coefficients_saturated.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# Strouhal number estimate from Cl peaks
peaks_t = []
for i in range(1, len(cl) - 1):
    if cl[i] > cl[i - 1] and cl[i] > cl[i + 1] and cl[i] > 0.5:
        peaks_t.append(t[i])
peaks_t = np.array(peaks_t)
if len(peaks_t) > 1:
    period = np.mean(np.diff(peaks_t))
    St = 1.0 / period
else:
    period = St = float("nan")

with open(OUT / "shedding_summary.txt", "w") as f:
    f.write("Vortex-shedding diagnostics, Python IBPM port, Re=100, nx=450 ny=200 (blockage-affected)\n")
    f.write(f"Number of Cl peaks detected in t=[{t[0]:.2f},{t[-1]:.2f}]: {len(peaks_t)}\n")
    f.write(f"Mean shedding period: {period:.4f}\n")
    f.write(f"Estimated Strouhal number St = 1/T: {St:.4f}\n")
    f.write("Literature St at Re=100, unbounded domain: ~0.164-0.170\n")
    f.write("This run's domain has a 25% blockage ratio (D=1 in a height-4 domain),\n")
    f.write("which is known to raise the shedding frequency relative to the unbounded case.\n")

print("period", period, "St", St)
print("wrote figures to", OUT)
