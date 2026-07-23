"""
analyze_followup_hm.py

Tests H-M: why is ibpm's Cd oscillation amplitude (fig13_14_hysteresis.png,
f=4Hz, alpha0=0) +123% vs. the paper's, much bigger than Cl's +31%? (see
../1-paper_based/README.md's "Open question" and this folder's README
"H-M" section). H and M need no new runs; I/J/K/L need
run_followup_hm.py's new runs first.

Usage: python3 analyze_followup_hm.py
Output: data/test_H_*.csv, data/test_I_*.csv, ..., data/test_M_*.csv
"""
import io
import pathlib

import numpy as np

REPO = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main")
KURT = REPO / "SURF_test" / "low_re" / "NACA0012" / "kurt_comp"
PAPER = KURT / "1-paper_based"
HERE = KURT / "2-follow_up"
RUNS_HM = HERE / "runs" / "hm"
BASE_F4HZ_A00 = PAPER / "runs" / "dx0.020" / "f4hz_py_a00"
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)

FREQ_F4HZ = 2.73972602739726


def read_csv_with_comments(path):
    lines = [l for l in pathlib.Path(path).read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    return np.genfromtxt(io.StringIO("\n".join(lines)), delimiter=",",
                          names=True, dtype=None, encoding="utf-8")


def load_force(run_dir, min_t=29.0):
    """Returns None unless the run has actually reached ~t=30 -- a
    still-running background job's partial flow.force otherwise gets
    silently picked up and produces garbage "last 2 periods" metrics from
    the very start of the transient."""
    f = pathlib.Path(run_dir) / "flow.force"
    if not f.exists():
        return None
    d = np.loadtxt(f)
    if d.ndim != 2 or d[-1, 1] < min_t:
        return None
    return d


def hysteresis_metrics(d, kt):
    """Peak-to-peak ratio (ibpm/paper) and RMS error (branch-matched,
    interpolated), same method as ../1-paper_based/gen_kurt_figs.py's
    _hysteresis_error_metrics -- duplicated here (not imported) since this
    folder deliberately doesn't import from 1-paper_based's script."""
    t, cd, cl = d[:, 1], d[:, 2], d[:, 3]
    period = 1.0 / FREQ_F4HZ
    m = t > (t[-1] - 2 * period)
    alpha = np.sin(2 * np.pi * FREQ_F4HZ * t[m])
    dadt = np.cos(2 * np.pi * FREQ_F4HZ * t[m])
    cd_m, cl_m = cd[m], cl[m]
    out = {}
    for name, paper_col, arr in [("Cl", "cl", cl_m), ("Cd", "cd", cd_m)]:
        errs = []
        for i in range(len(kt)):
            a0 = kt["alpha_deg"][i]
            br = kt["branch"][i]
            mask = (dadt > 0) if br == "up" else (dadt < 0)
            if mask.sum() < 2:
                continue
            idx = np.argsort(alpha[mask])
            val = np.interp(a0, alpha[mask][idx], arr[mask][idx])
            errs.append(val - kt[paper_col][i])
        errs = np.array(errs)
        paper_range = float(kt[paper_col].max() - kt[paper_col].min())
        ibpm_range = float(arr.max() - arr.min())
        out[name] = dict(ibpm_pk2pk=ibpm_range, paper_pk2pk=paper_range,
                          pk2pk_ratio=ibpm_range / paper_range,
                          rms_abs_error=float(np.sqrt(np.mean(errs ** 2))))
    return out


def main():
    kt = read_csv_with_comments(PAPER / "data" / "kurtulus_fig13_14_table.csv")
    fig1 = read_csv_with_comments(PAPER / "data" / "kurtulus_fig1_digitized.csv")
    mc = read_csv_with_comments(PAPER / "data" / "mean_coefficients_dx0.020.csv")

    baseline = load_force(BASE_F4HZ_A00)
    baseline_metrics = hysteresis_metrics(baseline, kt) if baseline is not None else None

    # ---------------- Test H: steady vs. dynamic baseline ----------------
    rows_h = []
    for a in range(0, 6):
        i_fig1 = int(np.argmin(np.abs(fig1["alpha_deg"] - a)))
        paper_steady = float(fig1["cd_steady"][i_fig1])
        paper_f4hz = float(fig1["cd_f4hz"][i_fig1])
        m_steady = (mc["motion"] == "steady") & (mc["impl"] == "py") & (mc["alpha_deg"] == a)
        m_f4hz = (mc["motion"] == "f4hz") & (mc["impl"] == "py") & (mc["alpha_deg"] == a)
        ibpm_steady = float(mc["cd_mean"][m_steady][0]) if m_steady.any() else float("nan")
        ibpm_f4hz = float(mc["cd_mean"][m_f4hz][0]) if m_f4hz.any() else float("nan")
        rows_h.append(dict(
            alpha_deg=a, ibpm_cd_steady=ibpm_steady, ibpm_cd_f4hz=ibpm_f4hz,
            paper_cd_steady=paper_steady, paper_cd_f4hz=paper_f4hz,
            pct_gap_steady=100 * (ibpm_steady - paper_steady) / paper_steady if paper_steady else float("nan"),
            pct_gap_f4hz=100 * (ibpm_f4hz - paper_f4hz) / paper_f4hz if paper_f4hz else float("nan"),
        ))
    import csv
    with open(DATA / "test_H_steady_vs_dynamic.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_h[0].keys()))
        w.writeheader(); w.writerows(rows_h)
    print("wrote test_H_steady_vs_dynamic.csv")
    print("\nTest H summary (digitization uncertainty on paper's Fig 1 is stated as +/-0.05):")
    for r in rows_h:
        print(f"  a={r['alpha_deg']}: ibpm steady={r['ibpm_cd_steady']:.4f} "
              f"(paper {r['paper_cd_steady']:.3f}, {r['pct_gap_steady']:+.0f}%), "
              f"ibpm f4hz_mean={r['ibpm_cd_f4hz']:.4f} (paper {r['paper_cd_f4hz']:.3f}, "
              f"{r['pct_gap_f4hz']:+.0f}%)")

    # ---------------- Tests I/J/K/L: each vs. the existing baseline ----------------
    variants = [
        ("I_ngrid2", "ngrid=2 (domain)"),
        ("I_ngrid3", "ngrid=3 (domain)"),
        ("J_LTEdense", "LE+TE boundary refined (ds=dx/4)"),
        ("K_dt0.0025", "dt refined (0.005->0.0025)"),
        ("L_Re1010", "Re +1% (1000->1010)"),
    ]
    rows_ijkl = []
    if baseline_metrics is not None:
        rows_ijkl.append(dict(case="baseline (dx=0.02,ngrid=1,dt=0.005,Re=1000)", impl="py",
                               cl_pk2pk_ratio=baseline_metrics["Cl"]["pk2pk_ratio"],
                               cd_pk2pk_ratio=baseline_metrics["Cd"]["pk2pk_ratio"],
                               cl_rms_err=baseline_metrics["Cl"]["rms_abs_error"],
                               cd_rms_err=baseline_metrics["Cd"]["rms_abs_error"]))
    # baseline cpp, for the same py-vs-cpp fidelity check done everywhere else
    d_base_cpp = load_force(PAPER / "runs" / "dx0.020" / "f4hz_cpp_a00")
    if d_base_cpp is not None:
        m_base_cpp = hysteresis_metrics(d_base_cpp, kt)
        rows_ijkl.append(dict(case="baseline (dx=0.02,ngrid=1,dt=0.005,Re=1000)", impl="cpp",
                               cl_pk2pk_ratio=m_base_cpp["Cl"]["pk2pk_ratio"],
                               cd_pk2pk_ratio=m_base_cpp["Cd"]["pk2pk_ratio"],
                               cl_rms_err=m_base_cpp["Cl"]["rms_abs_error"],
                               cd_rms_err=m_base_cpp["Cd"]["rms_abs_error"]))
    for name, label in variants:
        for impl, dirname in (("py", name), ("cpp", f"{name}_cpp")):
            d = load_force(RUNS_HM / dirname)
            if d is None:
                print(f"{dirname}: no data yet (run run_followup_hm.py first)")
                continue
            m = hysteresis_metrics(d, kt)
            rows_ijkl.append(dict(case=label, impl=impl,
                                   cl_pk2pk_ratio=m["Cl"]["pk2pk_ratio"], cd_pk2pk_ratio=m["Cd"]["pk2pk_ratio"],
                                   cl_rms_err=m["Cl"]["rms_abs_error"], cd_rms_err=m["Cd"]["rms_abs_error"]))
    if len(rows_ijkl) > 1:
        with open(DATA / "test_IJKL_knob_sensitivity.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_ijkl[0].keys()))
            w.writeheader(); w.writerows(rows_ijkl)
        print("wrote test_IJKL_knob_sensitivity.csv")
        print("\nTests I/J/K/L summary (pk2pk_ratio = ibpm/paper, 1.0 = perfect match):")
        for r in rows_ijkl:
            print(f"  {r['case']} ({r['impl']}): Cl ratio={r['cl_pk2pk_ratio']:.3f}, Cd ratio={r['cd_pk2pk_ratio']:.3f}")

    # ---------------- Test M: paper's own internal self-consistency ----------------
    # Fig 13/14's table rows are NOT evenly spaced in time (see the raw CSV --
    # gaps range from 0.01s to 0.03s), so a naive np.mean() over the listed
    # values double-counts the densely-sampled stretches; a time-weighted
    # (trapezoidal) average over the table's own t_s column is the fair
    # comparison to Fig 1's genuinely time-averaged mean coefficients.
    order = np.argsort(kt["t_s"])
    t_ord, cd_ord, cl_ord = kt["t_s"][order], kt["cd"][order], kt["cl"][order]
    cd_table_mean = float(np.trapz(cd_ord, t_ord) / (t_ord[-1] - t_ord[0]))
    cl_table_mean = float(np.trapz(cl_ord, t_ord) / (t_ord[-1] - t_ord[0]))
    i0 = int(np.argmin(np.abs(fig1["alpha_deg"] - 0)))
    fig1_cd_f4hz_a0 = float(fig1["cd_f4hz"][i0])
    fig1_cl_f4hz_a0 = float(fig1["cl_f4hz"][i0])
    rows_m = [dict(quantity="Cd", fig13_14_table_mean=cd_table_mean, fig1_digitized_value=fig1_cd_f4hz_a0,
                    abs_diff=cd_table_mean - fig1_cd_f4hz_a0,
                    within_digitization_uncertainty=abs(cd_table_mean - fig1_cd_f4hz_a0) <= 0.05),
              dict(quantity="Cl", fig13_14_table_mean=cl_table_mean, fig1_digitized_value=fig1_cl_f4hz_a0,
                    abs_diff=cl_table_mean - fig1_cl_f4hz_a0,
                    within_digitization_uncertainty=abs(cl_table_mean - fig1_cl_f4hz_a0) <= 0.05)]
    with open(DATA / "test_M_paper_selfcheck.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_m[0].keys()))
        w.writeheader(); w.writerows(rows_m)
    print("wrote test_M_paper_selfcheck.csv")
    print("\nTest M summary (paper Fig 13/14 table's own mean vs. paper's Fig 1 digitized point, "
          "both at alpha0=0, f=4Hz; Fig 1 stated digitization uncertainty is +/-0.05):")
    for r in rows_m:
        print(f"  {r['quantity']}: Fig13/14 table mean={r['fig13_14_table_mean']:.4f}, "
              f"Fig1 digitized={r['fig1_digitized_value']:.4f}, diff={r['abs_diff']:+.4f}, "
              f"within +/-0.05 digitization band: {r['within_digitization_uncertainty']}")


if __name__ == "__main__":
    main()
