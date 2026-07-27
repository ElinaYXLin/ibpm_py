"""
test4_spectral.py

Test 4: spectral character of the upstream noise. FFTs the upstream
window along x (per-row, then averaged) and in 2-D, and reports what
fraction of total power/enstrophy sits at the 2-cell (Nyquist)
wavelength -- `5-leading_edge`'s Test 0b already saw 1-2 cell
wavelengths at the LE itself; this checks whether the same odd-even
grid-decoupling signature extends into the upstream region, which would
name a specific, well-known numerical mechanism rather than leaving it
as generic "error". A cheap proxy (sign changes per unit length along
x, at y=0) is reported alongside since it requires no FFT machinery and
should tell the same story. NACA0012, alpha=0, steady, Re=1000,
dx=0.02/0.01/0.005, both implementations. Zero new runs.

Usage: python3 test4_spectral.py
Output: figures/test4_spectral.png, data/test4_spectral.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0

CASES = [
    dict(dx=0.02, py=c.KURT1 / "runs" / "dx0.020" / "steady_py_a00",
         cpp=c.KURT1 / "runs" / "dx0.020" / "steady_cpp_a00", nsteps=3000),
    dict(dx=0.01, py=c.KURT5 / "runs" / "grid_refine" / "dx0.0100",
         cpp=c.KURT5 / "runs" / "grid_refine" / "dx0.0100_cpp", nsteps=6000),
    dict(dx=0.005, py=c.KURT5 / "runs" / "grid_refine" / "dx0.0050",
         cpp=c.KURT5 / "runs" / "grid_refine" / "dx0.0050_cpp", nsteps=12000),
]
DX_COLOR = {0.02: "#1f77b4", 0.01: "#e67e22", 0.005: "#27ae60"}


def nyquist_fraction_2d(sub):
    """2-D FFT of the upstream window; fraction of total power in the
    outermost quarter-ring of frequency space (wavelengths < 4 cells,
    which brackets the 2-cell Nyquist mode this checks for)."""
    if sub.shape[0] < 4 or sub.shape[1] < 4:
        return float("nan")
    F = np.fft.fft2(sub - sub.mean())
    P = np.abs(F) ** 2
    nx, ny = sub.shape
    kx = np.fft.fftfreq(nx) * 2  # cycles per 2 cells (Nyquist = 1.0)
    ky = np.fft.fftfreq(ny) * 2
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    Kmag = np.sqrt(KX ** 2 + KY ** 2) / np.sqrt(2)  # normalize so Nyquist corner = 1
    high = Kmag > 0.75
    return float(P[high].sum() / P.sum())


def sign_changes_per_length(xs, vals_row, dx):
    s = np.sign(vals_row)
    s = s[s != 0]
    if s.size < 2:
        return 0.0
    n_changes = int(np.sum(np.abs(np.diff(s)) > 1))
    length = (s.size - 1) * dx
    return n_changes / length if length > 0 else 0.0


def main():
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    rows = []
    row_om = {}  # dx -> impl -> (xs, y0 row values) for the sign-change plot
    for case in CASES:
        dx = case["dx"]
        X, Y = c.grid_xy(dx)
        row_om[dx] = {}
        for impl in ("py", "cpp"):
            om = c.load_omega(case[impl], case["nsteps"])
            m = c.upstream_mask(X, x_le, BUFFER_DX, dx)
            ys = Y[0, :]
            iy = np.where((ys >= YLIM[0]) & (ys <= YLIM[1]))[0]
            sub = om[np.ix_(np.where(m)[0], iy)]
            frac = nyquist_fraction_2d(sub)

            iy0 = c.nearest_index(ys, 0.0)
            xs_row = X[m, 0]
            vals_row = om[np.ix_(np.where(m)[0], [iy0])].flatten()
            order = np.argsort(xs_row)
            sc = sign_changes_per_length(xs_row[order], vals_row[order], dx)
            row_om[dx][impl] = (xs_row[order], vals_row[order])

            rows.append(dict(dx=dx, impl=impl, nyquist_fraction=frac, sign_changes_per_chord=sc))
            print(f"dx={dx} {impl}: Nyquist-quarter power fraction={frac:.3f}, "
                  f"sign changes/chord at y=0 = {sc:.1f} (2/dx={2/dx:.0f} would be alternating every cell)")

    with open(c.DATA / "test4_spectral.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dx", "impl", "nyquist_fraction", "sign_changes_per_chord"])
        w.writeheader(); w.writerows(rows)
    print("wrote test4_spectral.csv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    dxs = sorted(row_om)
    for impl, marker in (("py", "o"), ("cpp", "x")):
        vals = [next(r["nyquist_fraction"] for r in rows if r["dx"] == dx and r["impl"] == impl) for dx in dxs]
        ax.plot(dxs, vals, marker + "-", label=impl, lw=1.8 if impl == "py" else 1.0,
                alpha=1.0 if impl == "py" else 0.6, ms=8)
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xlabel("dx (log, refining -->)"); ax.set_ylabel("fraction of upstream power at ~2-cell wavelength")
    ax.set_title("Nyquist-band power fraction", fontsize=10)
    ax.set_ylim(0, 1); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    ax = axes[1]
    for impl, marker in (("py", "o"), ("cpp", "x")):
        vals = [next(r["sign_changes_per_chord"] for r in rows if r["dx"] == dx and r["impl"] == impl) for dx in dxs]
        ax.plot(dxs, vals, marker + "-", label=impl, lw=1.8 if impl == "py" else 1.0,
                alpha=1.0 if impl == "py" else 0.6, ms=8)
    max_possible = [1.0 / dx for dx in dxs]  # every-cell alternation
    ax.plot(dxs, max_possible, "k--", lw=1, label="every-cell alternation (1/dx)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("dx (log, refining -->)"); ax.set_ylabel("sign changes per chord, y=0 row")
    ax.set_title("Cheap proxy: sign-change rate vs every-cell alternation", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    ax = axes[2]
    dx_show = 0.02
    for impl, ls in (("py", "-"), ("cpp", "--")):
        xs_row, vals_row = row_om[dx_show][impl]
        dist = x_le - xs_row
        ax.plot(dist, vals_row, ls, color=DX_COLOR[dx_show], label=impl,
                lw=1.6 if impl == "py" else 1.0, alpha=1.0 if impl == "py" else 0.6)
    ax.set_xlim(0, 0.3)
    ax.set_xlabel("distance upstream of LE (chord)"); ax.set_ylabel("omega at y=0")
    ax.set_title(f"Raw y=0 row, dx={dx_show} (visual check for cell-to-cell sign flips)", fontsize=10)
    ax.axhline(0, color="gray", lw=0.5)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)

    fig.suptitle("Test 4: spectral character of upstream noise -- is it odd-even grid decoupling?", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = c.FIGS / "test4_spectral.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
