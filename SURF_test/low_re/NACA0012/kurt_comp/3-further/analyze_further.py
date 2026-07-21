"""
analyze_further.py

The 7 further tests proposed as follow-up to 2-follow_up/README.md's "Open
questions". Zero-new-run tests (1a, 1b, 2a, 3a) work immediately from
1-paper_based/ and 2-follow_up/ data. Tests 2b, 2c, 3b need run_further.py's
new simulations first.

Usage: python3 analyze_further.py [zero|new|all]
Output: 3-further/data/*.csv, *.txt
"""
import sys
import pathlib
import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp"
PAPER_RUNS = KURT / "1-paper_based" / "runs" / "dx0.020"
FOLLOWUP_DATA = KURT / "2-follow_up" / "data"
FURTHER = KURT / "3-further"
FUP_RUNS = FURTHER / "runs"
DATA = FURTHER / "data"
DATA.mkdir(exist_ok=True)

FREQ = {"f1hz": 0.684931506849315, "f4hz": 2.73972602739726}
REP_ANGLES = [15, 20, 30, 40]
JAGGED_ANGLES = [25, 28, 30]
NEIGHBORS = {25: (24, 26), 28: (27, 29), 30: (29, 31)}
PERTURBATIONS = [-1e-2, -1e-3, 1e-2]


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
    return d


# ================================================================
# 1a. Full-record phase-averaged thrust window (many more cycles)
# ================================================================
def test_1a_full_phase_average():
    out = ["impl,alpha_deg,cd_phase_avg_min_4period,cd_phase_avg_min_fullrecord,n_bins_populated"]
    for impl in ("py",):
        for a in range(28, 43):
            run = PAPER_RUNS / f"f4hz_{impl}_a{a:02d}"
            d = load_force(run)
            if d is None:
                continue
            t, cd = d[:, 1], d[:, 2]
            f = FREQ["f4hz"]
            period = 1.0 / f
            tmax = t[-1]

            def phase_avg_min(mask):
                if mask.sum() < 20:
                    return float("nan"), 0
                phase = np.mod(2 * np.pi * f * t[mask], 2 * np.pi)
                nbins = 24
                bin_idx = np.clip((phase / (2 * np.pi) * nbins).astype(int), 0, nbins - 1)
                cd_m = cd[mask]
                populated = 0
                vals = []
                for k in range(nbins):
                    sel = bin_idx == k
                    if sel.any():
                        vals.append(cd_m[sel].mean())
                        populated += 1
                return (float(np.min(vals)) if vals else float("nan")), populated

            dev4 = t > (tmax - 4 * period)
            m4, _ = phase_avg_min(dev4)
            devfull = t > (tmax * 0.1)  # skip only the initial transient
            mfull, npop = phase_avg_min(devfull)
            out.append(f"{impl},{a},{m4:.5f},{mfull:.5f},{npop}")
    (DATA / "test1a_full_phase_average.csv").write_text("\n".join(out) + "\n")
    print("wrote test1a_full_phase_average.csv")


# ================================================================
# 1b. Does the dip band line up with a shedding-mode transition?
# ================================================================
def test_1b_dip_vs_transitions():
    # 1a's full-record phase-avg min
    d1a = np.genfromtxt(DATA / "test1a_full_phase_average.csv", delimiter=",",
                         names=True, dtype=None, encoding="utf-8")
    # 2-follow_up's fine Strouhal (steady) and 1-paper_based's mean Cl
    dst = np.genfromtxt(FOLLOWUP_DATA / "strouhal_fine_reanalysis.csv", delimiter=",",
                         names=True, dtype=None, encoding="utf-8")
    dmc = np.genfromtxt(KURT / "1-paper_based" / "data" / "mean_coefficients_dx0.020.csv",
                         delimiter=",", names=True, dtype=None, encoding="utf-8")

    out = ["alpha_deg,thrust_dip_cd,steady_strouhal_fine,steady_cl_mean"]
    for a in range(28, 43):
        row_dip = d1a[d1a["alpha_deg"] == a]
        dip = float(row_dip["cd_phase_avg_min_fullrecord"][0]) if len(row_dip) else float("nan")
        row_st = dst[(dst["impl"] == "py") & (dst["alpha_deg"] == a)]
        st = float(row_st["st_zeropad_interp"][0]) if len(row_st) else float("nan")
        row_cl = dmc[(dmc["motion"] == "steady") & (dmc["impl"] == "py") & (dmc["alpha_deg"] == a)]
        cl = float(row_cl["cl_mean"][0]) if len(row_cl) else float("nan")
        out.append(f"{a},{dip:.5f},{st:.4f},{cl:.5f}")
    (DATA / "test1b_dip_vs_transitions.csv").write_text("\n".join(out) + "\n")
    print("wrote test1b_dip_vs_transitions.csv")


# ================================================================
# 2a. Full spectral map (spectrogram): angle x frequency amplitude
# ================================================================
def test_2a_spectrogram():
    angles = list(range(0, 41)) + [50, 60]
    freq_grid = np.linspace(0.0, 1.0, 201)  # St from 0 to 1
    amp_map = np.full((len(angles), len(freq_grid)), np.nan)
    for i, a in enumerate(angles):
        run = PAPER_RUNS / f"steady_py_a{a:02d}"
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
        win = np.hanning(len(cl))
        amp = np.abs(np.fft.rfft(cl * win))
        freqs = np.fft.rfftfreq(len(cl), d=dt)
        if amp.max() > 0:
            amp_norm = amp / amp.max()
        else:
            amp_norm = amp
        amp_map[i] = np.interp(freq_grid, freqs, amp_norm, left=0, right=0)

    np.savez(DATA / "test2a_spectrogram.npz", angles=np.array(angles),
             freq_grid=freq_grid, amp_map=amp_map)
    # also a flat CSV for easy inspection
    lines = ["alpha_deg," + ",".join(f"{f:.3f}" for f in freq_grid)]
    for a, row in zip(angles, amp_map):
        lines.append(f"{a}," + ",".join(f"{v:.4f}" if not np.isnan(v) else "" for v in row))
    (DATA / "test2a_spectrogram.csv").write_text("\n".join(lines) + "\n")
    print("wrote test2a_spectrogram.csv/.npz")


# ================================================================
# 3a. Running-mean convergence: has t=30 actually converged?
# ================================================================
def test_3a_running_mean_convergence():
    out = ["impl,alpha_deg,cl_mean_at_t15,cl_mean_at_t22p5,cl_mean_at_t30,"
           "pct_change_last_third_vs_first_two_thirds"]
    for a in range(15, 41):
        run = PAPER_RUNS / f"steady_py_a{a:02d}"
        d = load_force(run)
        if d is None:
            continue
        t, cl = d[:, 1], d[:, 3]
        tmax = t[-1]

        def mean_up_to(frac):
            mask = t <= tmax * frac
            return float(np.mean(cl[mask])) if mask.sum() > 5 else float("nan")

        m_half = mean_up_to(0.5)     # ~t=15, using only the FIRST half (impulsive transient + a bit)
        m_75 = mean_up_to(0.75)      # ~t=22.5, first 75%
        m_full = float(np.mean(cl))  # full record mean (t=0-30)
        # convergence check: mean of last third vs mean of first two-thirds
        mask_last3rd = t > tmax * (2 / 3)
        mask_first2 = t <= tmax * (2 / 3)
        m_last3rd = float(np.mean(cl[mask_last3rd])) if mask_last3rd.sum() > 5 else float("nan")
        m_first2 = float(np.mean(cl[mask_first2])) if mask_first2.sum() > 5 else float("nan")
        pct = 100 * (m_last3rd - m_first2) / abs(m_first2) if m_first2 not in (0, float("nan")) else float("nan")
        out.append(f"py,{a},{m_half:.5f},{m_75:.5f},{m_full:.5f},{pct:.2f}")
    (DATA / "test3a_running_mean_convergence.csv").write_text("\n".join(out) + "\n")
    print("wrote test3a_running_mean_convergence.csv")


# ================================================================
# 2b. dx refinement of shedding frequency at representative angles
# ================================================================
def _steady_st(run_dir):
    d = load_force(run_dir)
    if d is None:
        return None
    n = len(d)
    seg = d[int(n * 0.5):]
    t, cl = seg[:, 1], seg[:, 3]
    if np.isnan(cl).any():
        return None
    cl = cl - np.mean(cl)
    if np.std(cl) < 5e-3:
        return 0.0
    dt = np.mean(np.diff(t))
    m = len(cl)
    win = np.hanning(m)
    nfft = 8 * m
    amp = np.abs(np.fft.rfft(cl * win, n=nfft))
    freqs = np.fft.rfftfreq(nfft, d=dt)
    band = (freqs > 0.05) & (freqs < 3.0)
    if not band.any() or np.max(amp[band]) < 5.0 * np.mean(amp[band]):
        return 0.0
    idx = np.flatnonzero(band)[np.argmax(amp[band])]
    if 0 < idx < len(amp) - 1:
        y0, y1, y2 = amp[idx - 1], amp[idx], amp[idx + 1]
        denom = y0 - 2 * y1 + y2
        delta = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-12 else 0.0
    else:
        delta = 0.0
    return max(0.0, float(freqs[idx] + delta * (freqs[1] - freqs[0])))


def test_2b_dx_refine():
    out = ["alpha_deg,st_dx0.02,st_dx0.01,st_dx0.005"]
    st005 = {}
    r = _steady_st(FUP_RUNS / "dx_refine" / "dx0.005_a20")
    if r is not None:
        st005[20] = r
    for a in REP_ANGLES:
        st02 = _steady_st(PAPER_RUNS / f"steady_py_a{a:02d}")
        st01 = _steady_st(FUP_RUNS / "dx_refine" / f"dx0.010_a{a:02d}")
        s02 = f"{st02:.4f}" if st02 is not None else ""
        s01 = f"{st01:.4f}" if st01 is not None else ""
        s005 = f"{st005[a]:.4f}" if a in st005 else ""
        out.append(f"{a},{s02},{s01},{s005}")
    (DATA / "test2b_dx_refine_strouhal.csv").write_text("\n".join(out) + "\n")
    print("wrote test2b_dx_refine_strouhal.csv")


# ================================================================
# 2c. ngrid sweep of shedding frequency at representative angles
# ================================================================
def test_2c_ngrid_strouhal():
    out = ["alpha_deg,st_ngrid1,st_ngrid2,st_ngrid3"]
    for a in REP_ANGLES:
        st1 = _steady_st(PAPER_RUNS / f"steady_py_a{a:02d}")
        st2 = _steady_st(FUP_RUNS / "ngrid_sweep" / f"ngrid2_a{a:02d}")
        st3 = _steady_st(FUP_RUNS / "ngrid_sweep" / f"ngrid3_a{a:02d}")
        s1 = f"{st1:.4f}" if st1 is not None else ""
        s2 = f"{st2:.4f}" if st2 is not None else ""
        s3 = f"{st3:.4f}" if st3 is not None else ""
        out.append(f"{a},{s1},{s2},{s3}")
    (DATA / "test2c_ngrid_strouhal.csv").write_text("\n".join(out) + "\n")
    print("wrote test2c_ngrid_strouhal.csv")


# ================================================================
# 3b. IC ensemble: multistability / hysteresis check
# ================================================================
def test_3b_ic_ensemble():
    out = ["alpha_deg,ic_type,cl_mean,cd_mean"]
    for a in JAGGED_ANGLES:
        run = PAPER_RUNS / f"steady_py_a{a:02d}"
        d = load_force(run)
        if d is not None:
            n = len(d)
            seg = d[int(n * 0.5):]
            out.append(f"{a},impulsive,{np.mean(seg[:,3]):.5f},{np.mean(seg[:,2]):.5f}")

        for direction in ("from_below", "from_above"):
            run = FUP_RUNS / "ic_ensemble" / "approach" / f"a{a:02d}_{direction}"
            d = load_force(run)
            if d is not None:
                n = len(d)
                seg = d[int(n * 0.5):]
                out.append(f"{a},{direction},{np.mean(seg[:,3]):.5f},{np.mean(seg[:,2]):.5f}")

        for rel in PERTURBATIONS:
            tag = f"rel{rel:+.0e}".replace("+", "p").replace("-", "m")
            run = FUP_RUNS / "ic_ensemble" / "perturb" / f"a{a:02d}_{tag}"
            d = load_force(run)
            if d is not None:
                n = len(d)
                seg = d[int(n * 0.5):]
                out.append(f"{a},perturb_{tag},{np.mean(seg[:,3]):.5f},{np.mean(seg[:,2]):.5f}")
    (DATA / "test3b_ic_ensemble.csv").write_text("\n".join(out) + "\n")
    print("wrote test3b_ic_ensemble.csv")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("zero", "all"):
        print("=== zero-new-run tests ===")
        test_1a_full_phase_average()
        test_1b_dip_vs_transitions()
        test_2a_spectrogram()
        test_3a_running_mean_convergence()
    if which in ("new", "all"):
        print("=== new-run-dependent tests ===")
        test_2b_dx_refine()
        test_2c_ngrid_strouhal()
        test_3b_ic_ensemble()
