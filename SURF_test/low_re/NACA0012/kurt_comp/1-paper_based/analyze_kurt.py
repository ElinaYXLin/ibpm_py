"""
analyze_kurt.py

Reduces the raw runs from run_kurt_suite.py into comparison tables vs Kurtulus
(2019). Writes CSV/txt into kurt_comp/data/. Run after the sweep completes (or
partway -- it uses whatever runs are present).

Outputs:
  mean_coefficients_<grid>.csv   mean Cl, Cd (+/-std) per (motion, impl, angle)
  shedding_strouhal_<grid>.csv   steady vortex-shedding Strouhal per (impl, angle)
  thrust_check_<grid>.csv        min instantaneous Cd for f4hz (paper: <0 for 3-37 deg)
  liftslope_<grid>.txt           dCl/dalpha near alpha=0 (paper: ~pi)
  fig13_14_comparison.csv        instantaneous Cl,Cd at alpha0=0 f4hz vs paper table

Usage: python3 SURF_test/low_re/NACA0012/kurt_comp/analyze_kurt.py
"""
import io
import pathlib
import numpy as np

KURT = pathlib.Path("/Users/elina/Desktop/SURF2026/ibpm_py-main/SURF_test/low_re/NACA0012/kurt_comp/1-paper_based")
RUNS = KURT / "runs"
DATA = KURT / "data"
DATA.mkdir(exist_ok=True)


def read_csv_with_comments(path):
    """np.genfromtxt(..., names=True, comments='#') mis-parses files that
    have '#'-prefixed description lines ABOVE a real (uncommented) header
    row -- strip full-line comments ourselves first, then hand genfromtxt a
    clean names=True table. Same fix as gen_kurt_figs.py's helper of the
    same name."""
    lines = [l for l in pathlib.Path(path).read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    return np.genfromtxt(io.StringIO("\n".join(lines)), delimiter=",",
                         names=True, dtype=None, encoding="utf-8")

FREQ = {"f1hz": 0.684931506849315, "f4hz": 2.73972602739726}
AVG_FRAC = 0.5  # average over the last 50% of the run


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


def tail(d):
    i0 = int(len(d) * (1 - AVG_FRAC))
    return d[i0:]


def mean_coeffs(d):
    seg = tail(d)
    cd, cl = seg[:, 2], seg[:, 3]
    return dict(cd_mean=float(np.mean(cd)), cl_mean=float(np.mean(cl)),
                cd_std=float(np.std(cd)), cl_std=float(np.std(cl)),
                cd_min=float(np.min(seg[:, 2])), nan=bool(np.isnan(seg).any()))


def shedding_strouhal(d):
    """Non-dim shedding freq (= Strouhal, since c=U=1) from FFT of Cl over the
    developed segment. Returns 0 if no clear peak (steady/attached)."""
    seg = tail(d)
    t, cl = seg[:, 1], seg[:, 3]
    if np.isnan(cl).any():
        return float("nan")
    cl = cl - np.mean(cl)
    n = len(cl)
    dt = np.mean(np.diff(t))
    if n < 32 or dt <= 0:
        return 0.0
    # require a developed oscillation at all (below shedding onset the wake is
    # steady/attached and Cl is essentially flat)
    if np.std(cl) < 5e-3:
        return 0.0
    win = np.hanning(n)
    amp = np.abs(np.fft.rfft(cl * win))
    freqs = np.fft.rfftfreq(n, d=dt)
    # search the plausible shedding band and require the peak to be prominent
    # (dominates the band mean) -- otherwise there is no real shedding tone
    band = (freqs > 0.1) & (freqs < 3.0)
    if not band.any():
        return 0.0
    pk_amp = np.max(amp[band])
    if pk_amp < 5.0 * np.mean(amp[band]):
        return 0.0
    return float(freqs[band][np.argmax(amp[band])])


def collect(grid):
    gdir = RUNS / grid
    if not gdir.exists():
        return
    rows_mc, rows_st, rows_th = [], [], []
    entries = {}
    for run in sorted(gdir.iterdir()):
        if not run.is_dir():
            continue
        try:
            motion, impl, atag = run.name.split("_")
            angle = int(atag[1:])
        except ValueError:
            continue
        d = load_force(run)
        if d is None:
            continue
        mc = mean_coeffs(d)
        entries[(motion, impl, angle)] = (d, mc)
        rows_mc.append((motion, impl, angle, mc["cl_mean"], mc["cl_std"],
                        mc["cd_mean"], mc["cd_std"], int(mc["nan"])))
        if motion == "steady":
            st = shedding_strouhal(d)
            rows_st.append((impl, angle, st))
        if motion == "f4hz":
            rows_th.append((impl, angle, mc["cd_min"]))

    if rows_mc:
        rows_mc.sort(key=lambda r: (r[0], r[1], r[2]))
        with open(DATA / f"mean_coefficients_{grid}.csv", "w") as f:
            f.write("motion,impl,alpha_deg,cl_mean,cl_std,cd_mean,cd_std,nan\n")
            for r in rows_mc:
                f.write("%s,%s,%d,%.5f,%.5f,%.5f,%.5f,%d\n" % r)
        print(f"wrote mean_coefficients_{grid}.csv ({len(rows_mc)} rows)")
    if rows_st:
        rows_st.sort(key=lambda r: (r[0], r[1]))
        with open(DATA / f"shedding_strouhal_{grid}.csv", "w") as f:
            f.write("impl,alpha_deg,strouhal_nondim\n")
            for r in rows_st:
                f.write("%s,%d,%.4f\n" % r)
        print(f"wrote shedding_strouhal_{grid}.csv ({len(rows_st)} rows)")
    if rows_th:
        rows_th.sort(key=lambda r: (r[0], r[1]))
        with open(DATA / f"thrust_check_{grid}.csv", "w") as f:
            f.write("impl,alpha_deg,cd_min\n")
            for r in rows_th:
                f.write("%s,%d,%.5f\n" % r)
        print(f"wrote thrust_check_{grid}.csv ({len(rows_th)} rows)")

    # lift-curve slope near alpha=0 (steady), per impl
    with open(DATA / f"liftslope_{grid}.txt", "w") as f:
        f.write("Lift-curve slope dCl/dalpha near alpha=0, steady NACA0012 Re=1000\n")
        f.write("(paper: approximately pi = %.4f per rad)\n\n" % np.pi)
        for impl in ("py", "cpp"):
            pts = [(a, entries[("steady", impl, a)][1]["cl_mean"])
                   for a in range(0, 6) if ("steady", impl, a) in entries]
            if len(pts) >= 3:
                aa = np.radians([p[0] for p in pts])
                cc = np.array([p[1] for p in pts])
                slope = np.polyfit(aa, cc, 1)[0]
                f.write("%s: dCl/dalpha = %.4f per rad  (%.3f x pi)\n"
                        % (impl, slope, slope / np.pi))
    print(f"wrote liftslope_{grid}.txt")

    # Fig 13/14 instantaneous comparison (alpha0=0, f4hz), dx0.020 only
    if grid == "dx0.020":
        fig1314_compare(entries)


def fig1314_compare(entries):
    kt = read_csv_with_comments(DATA / "kurtulus_fig13_14_table.csv")
    out_lines = ["branch,alpha_deg,kurt_cl,py_cl,cpp_cl,kurt_cd,py_cd,cpp_cd"]
    for impl in ("py", "cpp"):
        key = ("f4hz", impl, 0)
        if key not in entries:
            print(f"fig13_14: missing {key}, skipping {impl}")
            return
    for impl_data in [None]:  # single pass building combined table
        pass

    def branch_interp(impl, alpha_query, want_down):
        d, _ = entries[("f4hz", impl, 0)]
        t, cd, cl = d[:, 1], d[:, 2], d[:, 3]
        f = FREQ["f4hz"]
        # my convention: alpha(t) = A*sin(2 pi f t), A=1 deg; pitch-down = dalpha/dt<0
        phase = 2 * np.pi * f * t
        alpha_deg = np.sin(phase)  # amplitude 1 deg
        ddt = np.cos(phase)  # sign of dalpha/dt
        # use only the last two full periods (developed state)
        tmax = t[-1]
        period = 1.0 / f
        mask_dev = t > (tmax - 2 * period)
        down = mask_dev & (ddt < 0)
        up = mask_dev & (ddt > 0)
        sel = down if want_down else up
        if sel.sum() < 3:
            return np.nan, np.nan
        order = np.argsort(alpha_deg[sel])
        a_s = alpha_deg[sel][order]
        cl_i = np.interp(alpha_query, a_s, cl[sel][order])
        cd_i = np.interp(alpha_query, a_s, cd[sel][order])
        return cl_i, cd_i

    for row in kt:
        want_down = (row["branch"] == "down")
        py_cl, py_cd = branch_interp("py", row["alpha_deg"], want_down)
        cpp_cl, cpp_cd = branch_interp("cpp", row["alpha_deg"], want_down)
        out_lines.append("%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f" % (
            row["branch"], row["alpha_deg"], row["cl"], py_cl, cpp_cl,
            row["cd"], py_cd, cpp_cd))
    (DATA / "fig13_14_comparison.csv").write_text("\n".join(out_lines) + "\n")
    print("wrote fig13_14_comparison.csv")


if __name__ == "__main__":
    for grid in ("dx0.020", "dx0.010"):
        collect(grid)
