"""
analyze_followup.py

Follow-up analysis to 1-paper_based/README.md's "Anomalies" section. Two
kinds of test, per the mentor's ask (understand ibpm's limitations, not
propose a new method):

ZERO NEW RUNS -- reprocess the existing 1-paper_based/runs/dx0.020/ force
traces with a different (still standard) analysis choice, to check whether
an anomaly is a property of the ANALYSIS (averaging window, peak-picking,
a summary statistic like min) rather than of ibpm itself:
  A. Thrust window (f4hz): recompute using a phase-average over the
     developed cycles and the 5th-percentile, instead of the raw minimum,
     which is a single-outlier-sensitive statistic.
  B. Post-stall mean Cl/Cd jaggedness (steady): recompute the time-average
     window locked to an integer number of measured shedding periods,
     instead of a fixed last-50%-of-run window that can straddle partial
     cycles.
  C. Shedding Strouhal plateaus: recompute with a finer frequency
     resolution (zero-padded FFT + parabolic peak interpolation) instead
     of the raw FFT bin spacing used in 1-paper_based/analyze_kurt.py.

CHEAP STEADY RUNS -- reprocess 2-follow_up/runs/ (run_followup.py's output)
against the 1-paper_based/runs/dx0.020/steady_py_a00..05 baseline:
  D. ngrid sweep (far-field/blockage via the multi-domain method)
  E. domain-size sweep (far-field/blockage via a bigger uniform domain)
  F. grid-alignment (half-cell domain shift) at alpha=0
  G. dx refinement (dx=0.01) at alpha=0

Usage: python3 SURF_test/low_re/NACA0012/kurt_comp/2-follow_up/analyze_followup.py
Output: 2-follow_up/data/*.csv, *.txt
"""
import pathlib
import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp"
PAPER_RUNS = KURT / "1-paper_based" / "runs" / "dx0.020"
FOLLOWUP = KURT / "2-follow_up"
FUP_RUNS = FOLLOWUP / "runs"
DATA = FOLLOWUP / "data"
DATA.mkdir(exist_ok=True)

FREQ = {"f1hz": 0.684931506849315, "f4hz": 2.73972602739726}


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


# ============================================================
# A. Thrust window: phase-average + 5th-percentile vs raw min
# ============================================================
def test_a_thrust_window():
    out = ["impl,alpha_deg,cd_min,cd_p05,cd_phase_avg_min"]
    for impl in ("py", "cpp"):
        for a in list(range(0, 41)) + [50, 60]:
            run = PAPER_RUNS / f"f4hz_{impl}_a{a:02d}"
            d = load_force(run)
            if d is None:
                continue
            t, cd = d[:, 1], d[:, 2]
            n = len(d)
            seg = d[int(n * 0.4):]  # last 60%, same window as 1-paper_based
            cd_seg = seg[:, 2]
            cd_min = float(np.min(cd_seg))
            cd_p05 = float(np.percentile(cd_seg, 5))
            # phase-average: bin instantaneous Cd by phase within the pitch
            # cycle (developed portion only, last 2 periods), then take the
            # min of the PHASE-AVERAGED curve -- averages out chaotic
            # single-step spikes that dominate the raw min/percentile.
            f = FREQ["f4hz"]
            period = 1.0 / f
            tmax = t[-1]
            dev = t > (tmax - 4 * period)  # a few periods for a stable phase-avg
            if dev.sum() < 20:
                cd_phase_min = float("nan")
            else:
                phase = np.mod(2 * np.pi * f * t[dev], 2 * np.pi)
                nbins = 24
                bin_idx = np.clip((phase / (2 * np.pi) * nbins).astype(int), 0, nbins - 1)
                cd_dev = cd[dev]
                phase_avg = np.array([cd_dev[bin_idx == k].mean() if (bin_idx == k).any()
                                       else np.nan for k in range(nbins)])
                cd_phase_min = float(np.nanmin(phase_avg))
            out.append(f"{impl},{a},{cd_min:.5f},{cd_p05:.5f},{cd_phase_min:.5f}")
    (DATA / "thrust_window_reanalysis.csv").write_text("\n".join(out) + "\n")
    print("wrote thrust_window_reanalysis.csv")


# ============================================================
# B. Period-locked averaging for steady mean Cl/Cd
# ============================================================
def _shedding_freq(d):
    """Same peak-picking as 1-paper_based/analyze_kurt.py's
    shedding_strouhal, but returns the frequency (not just a prominence
    check) so it can be used to size an integer-period averaging window."""
    n = len(d)
    seg = d[int(n * 0.5):]
    t, cl = seg[:, 1], seg[:, 3]
    if np.isnan(cl).any():
        return None
    cl = cl - np.mean(cl)
    if np.std(cl) < 5e-3:
        return None  # steady/attached, no shedding
    dt = np.mean(np.diff(t))
    win = np.hanning(len(cl))
    amp = np.abs(np.fft.rfft(cl * win))
    freqs = np.fft.rfftfreq(len(cl), d=dt)
    band = (freqs > 0.1) & (freqs < 3.0)
    if not band.any():
        return None
    pk_amp = np.max(amp[band])
    if pk_amp < 5.0 * np.mean(amp[band]):
        return None
    return float(freqs[band][np.argmax(amp[band])])


def test_b_period_locked_mean():
    out = ["impl,alpha_deg,cl_mean_last50pct,cl_mean_periodlocked,cd_mean_last50pct,"
           "cd_mean_periodlocked,n_periods_used"]
    for impl in ("py", "cpp"):
        for a in list(range(0, 41)) + [50, 60]:
            run = PAPER_RUNS / f"steady_{impl}_a{a:02d}"
            d = load_force(run)
            if d is None:
                continue
            t, cd, cl = d[:, 1], d[:, 2], d[:, 3]
            n = len(d)
            seg50 = d[int(n * 0.5):]
            cl_50, cd_50 = float(np.mean(seg50[:, 3])), float(np.mean(seg50[:, 2]))

            fst = _shedding_freq(d)
            if fst is None or fst <= 0:
                # no clean shedding tone -- period-locking isn't meaningful,
                # fall back to the same last-50% window (n_periods=0 flags this)
                out.append(f"{impl},{a},{cl_50:.5f},{cl_50:.5f},{cd_50:.5f},{cd_50:.5f},0")
                continue
            period = 1.0 / fst
            tmax = t[-1]
            # use as many WHOLE periods as fit in the last 50% of the run
            avail = tmax * 0.5
            n_periods = max(1, int(avail / period))
            t0 = tmax - n_periods * period
            mask = t >= t0
            cl_pl, cd_pl = float(np.mean(cl[mask])), float(np.mean(cd[mask]))
            out.append(f"{impl},{a},{cl_50:.5f},{cl_pl:.5f},{cd_50:.5f},{cd_pl:.5f},{n_periods}")
    (DATA / "period_locked_mean_reanalysis.csv").write_text("\n".join(out) + "\n")
    print("wrote period_locked_mean_reanalysis.csv")


# ============================================================
# C. Finer-resolution Strouhal (zero-padding + parabolic peak interp)
# ============================================================
def test_c_fine_strouhal():
    out = ["impl,alpha_deg,st_raw_bin,st_zeropad_interp"]
    for impl in ("py", "cpp"):
        for a in list(range(0, 41)) + [50, 60]:
            run = PAPER_RUNS / f"steady_{impl}_a{a:02d}"
            d = load_force(run)
            if d is None:
                continue
            n = len(d)
            seg = d[int(n * 0.5):]
            t, cl = seg[:, 1], seg[:, 3]
            if np.isnan(cl).any():
                continue
            cl = cl - np.mean(cl)
            dt = np.mean(np.diff(t))
            if np.std(cl) < 5e-3:
                out.append(f"{impl},{a},0.0000,0.0000")
                continue
            m = len(cl)
            win = np.hanning(m)

            # raw (same as 1-paper_based): FFT bin resolution = 1/(m*dt)
            amp_raw = np.abs(np.fft.rfft(cl * win))
            freqs_raw = np.fft.rfftfreq(m, d=dt)
            band_raw = (freqs_raw > 0.1) & (freqs_raw < 3.0)
            if not band_raw.any() or np.max(amp_raw[band_raw]) < 5.0 * np.mean(amp_raw[band_raw]):
                out.append(f"{impl},{a},0.0000,0.0000")
                continue
            st_raw = float(freqs_raw[band_raw][np.argmax(amp_raw[band_raw])])

            # zero-padded to 8x length -> 8x finer bin spacing, then
            # parabolic (quadratic) interpolation across the 3 bins around
            # the peak for sub-bin accuracy
            nfft = 8 * m
            amp = np.abs(np.fft.rfft(cl * win, n=nfft))
            freqs = np.fft.rfftfreq(nfft, d=dt)
            band = (freqs > 0.1) & (freqs < 3.0)
            k = np.argmax(amp[band])
            idx = np.flatnonzero(band)[k]
            if 0 < idx < len(amp) - 1:
                y0, y1, y2 = amp[idx - 1], amp[idx], amp[idx + 1]
                denom = (y0 - 2 * y1 + y2)
                delta = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-12 else 0.0
            else:
                delta = 0.0
            st_fine = float(freqs[idx] + delta * (freqs[1] - freqs[0]))
            out.append(f"{impl},{a},{st_raw:.4f},{st_fine:.4f}")
    (DATA / "strouhal_fine_reanalysis.csv").write_text("\n".join(out) + "\n")
    print("wrote strouhal_fine_reanalysis.csv")


# ============================================================
# D/E/F/G: cheap steady runs (2-follow_up/runs/)
# ============================================================
def _lift_slope(pts):
    """pts: list of (alpha_deg, cl_mean). Fit dCl/dalpha (per rad) through
    the origin-anchored small-angle points, same method as 1-paper_based."""
    if len(pts) < 3:
        return None
    aa = np.radians([p[0] for p in pts])
    cc = np.array([p[1] for p in pts])
    return float(np.polyfit(aa, cc, 1)[0])


def _mean_cl(run_dir):
    d = load_force(run_dir)
    if d is None:
        return None
    n = len(d)
    seg = d[int(n * 0.5):]
    return float(np.mean(seg[:, 3])), float(np.mean(seg[:, 2]))


def test_d_ngrid_sweep():
    out = ["ngrid,alpha_deg,cl_mean,cd_mean"]
    pts_by_ngrid = {}
    # ngrid=1 baseline reused from 1-paper_based (py only)
    pts1 = []
    for a in [0, 1, 2, 3, 4, 5]:
        run = PAPER_RUNS / f"steady_py_a{a:02d}"
        r = _mean_cl(run)
        if r:
            out.append(f"1,{a},{r[0]:.5f},{r[1]:.5f}")
            pts1.append((a, r[0]))
    pts_by_ngrid[1] = pts1
    for ngrid in (2, 3):
        pts = []
        for a in [0, 1, 2, 3, 4, 5]:
            run = FUP_RUNS / "ngrid_sweep" / f"ngrid{ngrid}_a{a:02d}"
            r = _mean_cl(run)
            if r:
                out.append(f"{ngrid},{a},{r[0]:.5f},{r[1]:.5f}")
                pts.append((a, r[0]))
        pts_by_ngrid[ngrid] = pts
    (DATA / "ngrid_sweep_reanalysis.csv").write_text("\n".join(out) + "\n")
    print("wrote ngrid_sweep_reanalysis.csv")

    slope_lines = ["Lift-curve slope dCl/dalpha vs ngrid (far-field/blockage treatment),",
                   "steady NACA0012 Re=1000, dx=0.02, py_static only",
                   "(paper: ~pi = %.4f per rad; 1-paper_based ngrid=1 baseline: 3.5720)\n" % np.pi]
    for ngrid in (1, 2, 3):
        slope = _lift_slope(pts_by_ngrid[ngrid])
        if slope is not None:
            slope_lines.append(f"ngrid={ngrid}: dCl/dalpha = {slope:.4f} per rad  ({slope/np.pi:.3f} x pi)")
    (DATA / "ngrid_liftslope.txt").write_text("\n".join(slope_lines) + "\n")
    print("wrote ngrid_liftslope.txt")


def test_e_domain_sweep():
    out = ["domain,alpha_deg,cl_mean,cd_mean"]
    pts_by_domain = {}
    pts_base = []
    for a in [0, 1, 2, 3, 4, 5]:
        run = PAPER_RUNS / f"steady_py_a{a:02d}"
        r = _mean_cl(run)
        if r:
            out.append(f"baseline_L6,{a},{r[0]:.5f},{r[1]:.5f}")
            pts_base.append((a, r[0]))
    pts_by_domain["baseline_L6"] = pts_base
    pts_big = []
    for a in [0, 1, 2, 3, 4, 5]:
        run = FUP_RUNS / "domain_sweep" / f"large_a{a:02d}"
        r = _mean_cl(run)
        if r:
            out.append(f"large_L10,{a},{r[0]:.5f},{r[1]:.5f}")
            pts_big.append((a, r[0]))
    pts_by_domain["large_L10"] = pts_big
    (DATA / "domain_sweep_reanalysis.csv").write_text("\n".join(out) + "\n")
    print("wrote domain_sweep_reanalysis.csv")

    slope_lines = ["Lift-curve slope dCl/dalpha vs domain size (blockage), ngrid=1,",
                   "steady NACA0012 Re=1000, dx=0.02, py_static only",
                   "(paper: ~pi = %.4f per rad)\n" % np.pi]
    for label in ("baseline_L6", "large_L10"):
        slope = _lift_slope(pts_by_domain[label])
        if slope is not None:
            slope_lines.append(f"{label}: dCl/dalpha = {slope:.4f} per rad  ({slope/np.pi:.3f} x pi)")
    (DATA / "domain_liftslope.txt").write_text("\n".join(slope_lines) + "\n")
    print("wrote domain_liftslope.txt")


def test_f_grid_alignment():
    base = PAPER_RUNS / "steady_py_a00"
    shifted = FUP_RUNS / "grid_alignment" / "yshift_a00"
    rb, rs = _mean_cl(base), _mean_cl(shifted)
    lines = ["Grid-alignment check: Cl(alpha=0) with domain shifted by dx/2=0.01 in y,",
             "steady NACA0012 Re=1000, dx=0.02, py_static only\n"]
    if rb:
        lines.append(f"baseline   (yoffset=-1.50): Cl(0) = {rb[0]:.5f}")
    if rs:
        lines.append(f"half-cell-shifted (yoffset=-1.49): Cl(0) = {rs[0]:.5f}")
    if rb and rs:
        lines.append(f"\ndifference: {rs[0]-rb[0]:.5f}  (sign flip: {np.sign(rb[0])!=np.sign(rs[0])})")
    (DATA / "grid_alignment_check.txt").write_text("\n".join(lines) + "\n")
    print("wrote grid_alignment_check.txt")


def test_g_dx_refine():
    base = PAPER_RUNS / "steady_py_a00"
    fine = FUP_RUNS / "dx_refine" / "dx0.010_a00"
    rb, rf = _mean_cl(base), _mean_cl(fine)
    lines = ["dx refinement check: Cl(alpha=0), steady NACA0012 Re=1000, py_static only\n"]
    if rb:
        lines.append(f"dx=0.02 (nx=300,ny=150): Cl(0) = {rb[0]:.5f}")
    if rf:
        lines.append(f"dx=0.01 (nx=600,ny=300): Cl(0) = {rf[0]:.5f}")
    if rb and rf:
        shrink = (1 - abs(rf[0]) / abs(rb[0])) * 100 if rb[0] != 0 else float("nan")
        lines.append(f"\n|Cl(0)| shrink from dx=0.02 to dx=0.01: {shrink:.1f}%")
    (DATA / "dx_refine_check.txt").write_text("\n".join(lines) + "\n")
    print("wrote dx_refine_check.txt")


if __name__ == "__main__":
    print("=== zero-new-runs reanalysis ===")
    test_a_thrust_window()
    test_b_period_locked_mean()
    test_c_fine_strouhal()
    print("\n=== cheap-steady-runs reanalysis ===")
    test_d_ngrid_sweep()
    test_e_domain_sweep()
    test_f_grid_alignment()
    test_g_dx_refine()
