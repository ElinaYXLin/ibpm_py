"""
test4_eigenmode_projection.py

Follow-up test #4 from ../7-stripe_investigation/README.md's Proposals
section: Test 8 (in that folder) showed the noise is seeded by the
projection step; Test 6/this folder's Test 2 showed severity tracks
cond(M) -- but nothing yet shows the noise specifically lives in M's
ill-conditioned (small-eigenvalue) modes.

Judgment call worth flagging: the OBVIOUS first check -- "how much of
the converged boundary force f's energy sits in M's small-eigenvalue
modes" -- turns out not to be diagnostic on its own. Since the
projection step solves f = M^-1 * rhs, and M^-1's eigenvalues are 1/lambda,
f is algebraically GUARANTEED to load disproportionately onto small
eigenvalues for any reasonably-spread rhs, entirely independent of
whether those modes are "causing a problem" -- this folder confirmed
that directly (99.8% of f's energy sits in the smallest 10% of
eigenvalues) and it would come out looking exactly the same even in a
perfectly healthy solve. So this script instead runs the sharper, more
direct test: reconstruct the boundary force using ONLY the largest-K
eigenvalue (well-conditioned) modes, spread that truncated force with
model.B() (the same operator Test 8 used), and check whether the
upstream footprint specifically shrinks as the ill-conditioned tail is
dropped, while the near-body force pattern is preserved. This doubles
as a direct proof-of-concept for the "truncate M's small eigenvalues"
mitigation proposed in ../7-stripe_investigation/README.md.

Usage: python3 test4_eigenmode_projection.py
Output: figures/test4_eigenmode_projection.png,
        data/test4_eigenmode_projection.csv
"""
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

DX = 0.02
BUFFER_DX = 2.0
YLIM = (-0.5, 0.5)
NSTEPS = 3000
RUN_DIR = c.KURT1 / "runs" / "dx0.020" / "steady_py_a00"
KEEP_FRACTIONS = [1.0, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01]  # fraction of LARGEST-eigenvalue modes kept


def main():
    grid = c.Grid(300, 150, 1, c.DOMAIN["length"], c.DOMAIN["xoffset"], c.DOMAIN["yoffset"])
    geom = c.Geometry(str(c.BASE_GEOM_DX002))
    M, numPoints = c.build_projection_matrix(grid, geom, c.RE)
    eigvals, eigvecs = np.linalg.eigh(M)  # ascending: eigvals[0] smallest (worst-conditioned)
    n = len(eigvals)
    print(f"M is {M.shape}, numPoints={numPoints}, "
          f"eigenvalue range [{eigvals.min():.4e}, {eigvals.max():.4e}], "
          f"cond={eigvals.max()/eigvals.min():.3e}")

    state = c.load_state(RUN_DIR, NSTEPS)
    f_flat = state.f._data.copy()
    f_norm2 = float(np.sum(f_flat ** 2))
    coeffs = eigvecs.T @ f_flat
    energy_frac = coeffs ** 2 / f_norm2
    print(f"(context, not diagnostic on its own -- see docstring): "
          f"bottom 10% of eigenvalues carry {energy_frac[:n//10].sum()*100:.1f}% of f's energy")

    model = c.NavierStokesModel(grid, geom, c.RE)
    model.init()
    X, Y = c.grid_xy(DX)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    def spread(vec_flat):
        bv = c.BoundaryVector(numPoints)
        bv._data[:] = vec_flat
        out = c.Scalar(grid)
        model.B(bv, out)
        return out._data[0].copy()

    full_footprint = spread(f_flat)

    # note (superseding an earlier version of this script): "upstream
    # enstrophy" of the raw spread footprint is trivially 0 at EVERY
    # truncation level -- consistent with, not contradicting, Test 8
    # Probe A (../7-stripe_investigation), which already showed the raw
    # B()-spread footprint has zero support beyond ~1-2 cells from the LE
    # regardless of the boundary force's content. That's the wrong region
    # to look in here. The informative comparison is WITHIN the compact
    # support near the LE (x_le to x_le+0.3c): does the oscillatory
    # (checkerboard-like) content there specifically shrink under
    # truncation, while the bulk/smooth force magnitude survives?
    ys = Y[0, :]
    iy0 = c.nearest_index(ys, 0.0)
    nb_mask = (X[:, 0] >= x_le) & (X[:, 0] <= x_le + 0.3)

    rows = []
    footprints = {}
    for keep_frac in KEEP_FRACTIONS:
        k = max(1, int(round(n * keep_frac)))
        # keep the k LARGEST-eigenvalue modes (indices n-k .. n-1), zero the rest
        mask = np.zeros(n, dtype=bool)
        mask[n - k:] = True
        coeffs_trunc = np.where(mask, coeffs, 0.0)
        f_trunc = eigvecs @ coeffs_trunc
        footprint = spread(f_trunc)
        footprints[keep_frac] = footprint

        row_vals = footprint[np.ix_(np.where(nb_mask)[0], [iy0])].flatten()
        s = np.sign(row_vals)
        s = s[s != 0]
        autocorr = float(np.mean(s[:-1] * s[1:])) if len(s) > 1 else float("nan")
        nb_peak = float(np.abs(footprint[nb_mask, :]).max())
        nb_bulk = float(np.sum(np.abs(footprint[nb_mask, :])) * DX * DX)  # smooth/bulk magnitude, less sensitive to sign-flip noise
        energy_kept = float(coeffs_trunc[mask] @ coeffs_trunc[mask] / f_norm2)
        rows.append(dict(keep_fraction=keep_frac, n_modes_kept=k, energy_kept_fraction=energy_kept,
                          near_le_lag1_autocorr=autocorr, near_body_peak=nb_peak, near_body_bulk=nb_bulk))
        print(f"keep top {keep_frac*100:.0f}% of modes ({k}/{n}): energy kept={energy_kept*100:.1f}%, "
              f"near-LE y=0 lag-1 autocorr={autocorr:.3f}, near-body peak={nb_peak:.3f}, "
              f"near-body bulk (sum|f|*dA)={nb_bulk:.4f}")

    with open(c.DATA / "test4_eigenmode_projection.csv", "w", newline="") as fcsv:
        w = csv.DictWriter(fcsv, fieldnames=["keep_fraction", "n_modes_kept", "energy_kept_fraction",
                                               "near_le_lag1_autocorr", "near_body_peak", "near_body_bulk"])
        w.writeheader(); w.writerows(rows)
    print("wrote test4_eigenmode_projection.csv")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    ax = axes[0, 0]
    kf = [r["keep_fraction"] * 100 for r in rows]
    ax.plot(kf, [r["near_le_lag1_autocorr"] for r in rows], "o-", color="#c0392b", lw=1.8, ms=7)
    ax.axhline(-1, color="black", ls=":", lw=1, alpha=0.6, label="-1 = perfect checkerboard")
    ax.set_xscale("log")
    ax.set_xlabel("% of modes kept (largest-eigenvalue first)")
    ax.set_ylabel("lag-1 sign autocorrelation, near-LE row")
    ax.set_title("Does dropping the ill-conditioned tail\nremove the near-LE checkerboard?", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(kf, [r["near_body_bulk"] for r in rows], "o-", color="#2980b9", lw=1.8, ms=7,
             label="near-body bulk magnitude (sum|f|*dA)")
    ax.axhline(rows[0]["near_body_bulk"], color="gray", ls="--", lw=1,
               label="full (untruncated) value")
    ax.set_xscale("log")
    ax.set_xlabel("% of modes kept (largest-eigenvalue first)")
    ax.set_ylabel("near-body bulk |footprint| magnitude")
    ax.set_title("...while the bulk near-body force survives?", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    for ax, keep_frac, title in [(axes[1, 0], 1.0, "Full (untruncated) force, spread"),
                                    (axes[1, 1], 0.1, "Top 10% of modes only, spread")]:
        Xw, Yw, fw, _, _ = c.window(X, Y, footprints[keep_frac], (-0.5, 1.35), (-0.35, 0.35))
        V = np.abs(fw).max() or 1.0
        ax.pcolormesh(Xw, Yw, np.clip(fw[:-1, :-1], -V, V), shading="flat", cmap="jet", vmin=-V, vmax=V)
        ax.set_aspect("equal"); ax.set_title(title, fontsize=10)

    fig.suptitle("Test 4: truncating M's ill-conditioned modes from the boundary force --\n"
                 "does the upstream footprint go away while the real force survives?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = c.FIGS / "test4_eigenmode_projection.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
