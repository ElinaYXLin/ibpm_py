"""
test3_checkerboard_quantified.py

Follow-up test #3 from ../7-stripe_investigation/README.md's Proposals
section: Test 4's bulk Nyquist-fraction metric came out weak (~5%) only
because the alternating pattern rides on top of a smooth decaying
envelope, which spreads the signal's spectral power across many bands.
This divides out the envelope first (a local moving-RMS normalization
along each upstream row), then recomputes the Nyquist-band power
fraction on the ENVELOPE-NORMALIZED signal, and separately computes the
lag-1 sign autocorrelation (a simpler, complementary proxy: a pure
checkerboard has lag-1 sign autocorrelation of exactly -1). Turns
"visually looks like a checkerboard" into an actual number. Zero new
runs -- reuses the same dx=0.02/0.01/0.005 fields as
../7-stripe_investigation Tests 1/4.

Usage: python3 test3_checkerboard_quantified.py
Output: figures/test3_checkerboard_quantified.png,
        data/test3_checkerboard_quantified.csv
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
    dict(dx=0.02, run=c.KURT1 / "runs" / "dx0.020" / "steady_py_a00", nsteps=3000),
    dict(dx=0.01, run=c.KURT5 / "runs" / "grid_refine" / "dx0.0100", nsteps=6000),
    dict(dx=0.005, run=c.KURT5 / "runs" / "grid_refine" / "dx0.0050", nsteps=12000),
]


def envelope_normalize(vals, window=5):
    """Divide the signal by a local moving RMS (window cells), to remove
    the smooth decaying envelope and isolate the oscillatory part. Cells
    where the local RMS is ~0 are left as 0 (avoid divide-by-~0 blowup)."""
    n = len(vals)
    out = np.zeros(n)
    half = window // 2
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        local_rms = np.sqrt(np.mean(vals[lo:hi] ** 2))
        out[i] = vals[i] / local_rms if local_rms > 1e-8 else 0.0
    return out


def nyquist_fraction_1d(vals):
    if len(vals) < 4:
        return float("nan")
    F = np.fft.rfft(vals - vals.mean())
    P = np.abs(F) ** 2
    # top quarter of the frequency band = near-Nyquist (short wavelength)
    n_hi = max(1, len(P) // 4)
    return float(P[-n_hi:].sum() / P.sum()) if P.sum() > 0 else float("nan")


def lag1_sign_autocorr(vals):
    s = np.sign(vals)
    s = s[s != 0]
    if len(s) < 2:
        return float("nan")
    return float(np.mean(s[:-1] * s[1:]))


def main():
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    rows = []
    example_raw = example_norm = example_dx = None
    for case in CASES:
        dx = case["dx"]
        X, Y = c.grid_xy(dx)
        om = c.load_omega(case["run"], case["nsteps"])
        ys = Y[0, :]
        iy0 = c.nearest_index(ys, 0.0)
        mask = c.upstream_mask(X, x_le, BUFFER_DX, dx)
        xs = X[mask, 0]
        order = np.argsort(xs)
        row_vals = om[np.ix_(np.where(mask)[0], [iy0])].flatten()[order]

        raw_nyq = nyquist_fraction_1d(row_vals)
        raw_autocorr = lag1_sign_autocorr(row_vals)
        norm_vals = envelope_normalize(row_vals)
        norm_nyq = nyquist_fraction_1d(norm_vals)
        norm_autocorr = lag1_sign_autocorr(norm_vals)

        rows.append(dict(dx=dx, raw_nyquist_fraction=raw_nyq, raw_lag1_autocorr=raw_autocorr,
                          envelope_norm_nyquist_fraction=norm_nyq, envelope_norm_lag1_autocorr=norm_autocorr))
        print(f"dx={dx}: raw Nyquist frac={raw_nyq:.3f}, raw lag-1 autocorr={raw_autocorr:.3f} | "
              f"envelope-normalized Nyquist frac={norm_nyq:.3f}, "
              f"envelope-normalized lag-1 autocorr={norm_autocorr:.3f}")

        if dx == 0.01:
            # dx=0.01 is the most informative example: raw autocorr (~0.06)
            # looks like NOT a checkerboard, envelope-normalized (-1.0)
            # proves it IS one -- dx=0.02 doesn't show this contrast since
            # its raw signal was already clean.
            example_raw, example_norm, example_dx = row_vals, norm_vals, dx

    with open(c.DATA / "test3_checkerboard_quantified.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dx", "raw_nyquist_fraction", "raw_lag1_autocorr",
                                            "envelope_norm_nyquist_fraction", "envelope_norm_lag1_autocorr"])
        w.writeheader(); w.writerows(rows)
    print("wrote test3_checkerboard_quantified.csv")

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    ax = axes[0]
    dxs = [r["dx"] for r in rows]
    ax.plot(dxs, [r["raw_nyquist_fraction"] for r in rows], "o-", color="gray", label="raw", lw=1.6, ms=7)
    ax.plot(dxs, [r["envelope_norm_nyquist_fraction"] for r in rows], "o-", color="#c0392b",
             label="envelope-normalized", lw=1.8, ms=7)
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xlabel("dx (log, refining -->)"); ax.set_ylabel("Nyquist-band power fraction")
    ax.set_title("Nyquist fraction: raw vs. envelope-normalized", fontsize=10)
    ax.set_ylim(0, 1); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(dxs, [r["raw_lag1_autocorr"] for r in rows], "o-", color="gray", label="raw", lw=1.6, ms=7)
    ax.plot(dxs, [r["envelope_norm_lag1_autocorr"] for r in rows], "o-", color="#2980b9",
             label="envelope-normalized", lw=1.8, ms=7)
    ax.axhline(-1, color="black", ls=":", lw=1, alpha=0.6, label="-1 = perfect checkerboard")
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xlabel("dx (log, refining -->)"); ax.set_ylabel("lag-1 sign autocorrelation")
    ax.set_title("Lag-1 sign autocorrelation (-1 = alternates every cell)", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    ax = axes[2]
    # row arrays are ordered by increasing x (index 0 = far upstream domain
    # edge, last index = nearest the LE, since the mask spans domain_min to
    # x_le-buffer) -- take the LAST n_show entries and reverse them so index
    # 0 = nearest the LE, matching the axis label
    n_show = min(60, len(example_norm))
    raw_near_le = example_raw[-n_show:][::-1]
    norm_near_le = example_norm[-n_show:][::-1]
    ax.plot(np.arange(n_show), raw_near_le / (np.abs(raw_near_le).max() or 1),
             "-", color="gray", lw=1.2, label="raw (normalized to its own peak)")
    ax.plot(np.arange(n_show), norm_near_le / (np.abs(norm_near_le).max() or 1),
             "-", color="#c0392b", lw=1.4, label="envelope-normalized")
    ax.set_xlabel("cell index (from LE outward)"); ax.set_ylabel("normalized value")
    ax.set_title(f"Example row, dx={example_dx}\n(envelope-normalization isolates the oscillation)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)

    fig.suptitle("Test 3: is the upstream noise a checkerboard pattern? (quantified)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = c.FIGS / "test3_checkerboard_quantified.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
