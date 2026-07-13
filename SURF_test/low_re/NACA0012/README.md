# NACA0012 at low Reynolds number vs. published CFD-benchmark drag (non-LSAT)

**Note on this folder's location:** this validation originally lived at
`airfoils/Lockard-NACA0012/`, named per the `<dataset>-<airfoil>`
convention used for the LSAT wind-tunnel cases (see below). It has since
been merged into `low_re/NACA0012/` (alongside the qualitative
`flow_evolution.png` already here, from the same Re=500 case) so that all
low-Re NACA0012 results -- qualitative flow-field and quantitative
Cl/Cd -- live in one place instead of being split across two directories.
Nothing about the data or analysis changed, only its location; see
`../README.md` for how this fits alongside `../SD7003/`.

This is the suite's only validation case against a **non-LSAT reference
dataset**: published low-Reynolds-number *computational* benchmark drag
coefficients for the NACA0012 airfoil, at Reynolds numbers genuinely in
the hundreds. Every airfoil in `../../airfoils/` (`LSAT-SD7003`,
`LSAT-SD8000`, `LSAT-ClarkY`, `LSAT-GM15`) validates against the UIUC LSAT
wind-tunnel dataset; this one cannot, because **no wind-tunnel data exists
at Re in the hundreds** (UIUC LSAT's practical floor is Re~40,000-60,000
-- see `../../airfoils/README.md` and `../README.md`). Below that floor,
the only published Cl/Cd data is CFD/DNS-computed. That is exactly the
dataset used here.

## Folder-naming convention (for the LSAT cases)

To make the reference-dataset provenance explicit and consistent, every
airfoil folder under `../../airfoils/` is named `<dataset>-<airfoil>`:

- `LSAT-SD7003`, `LSAT-SD8000`, `LSAT-ClarkY`, `LSAT-GM15` -- validated
  against the UIUC **LSAT** (Low-Speed Airfoil Tests) wind-tunnel dataset.
- This folder (formerly `Lockard-NACA0012`) -- validated against the
  **Lockard et al.** low-Re CFD drag benchmark (and corroborating
  references below), now living under `low_re/` instead since it isn't an
  LSAT case.

## The reference dataset

At alpha=0, a symmetric airfoil produces zero lift, so the benchmark is a
pure **drag** comparison there. The published low-Re NACA0012 drag
coefficients (all computational; gathered from the modern low-Re CFD
literature) are:

| Re | Cd (alpha=0) | Source |
|---|---|---|
| 500 | **0.1762** | Lockard, Luo, Milder & Singer |
| 500 | 0.1759 | Wu et al. |
| 500 | 0.178 | Nita et al. (two-relaxation-time LBM, arXiv:1901.08766) |
| 1000 | **0.119** | Di Ilio et al. (hybrid LBM), arXiv:2006.10487 |
| 1000 | 0.119 | Di Ilio et al. (XFOIL, same paper) |
| 1000 | ~0.12 | Kurtulus 2015 (Int. J. Micro Air Vehicles) |

The Re=500 references agree to within 0.2% of each other, and the Re=1000
references to within ~1% -- a tight, well-established benchmark.

For the lift side, this folder runs a small angle-of-attack polar at
Re=500 (alpha = 0, 2, 4, 6, 8, 10 deg) so both a Cl(alpha) and a
Cd(alpha) curve are produced. Published Re=500 *lift* data at nonzero
alpha exists only in figures (Kurtulus, Di Ilio), not machine-tabulated,
so the quantitative anchor is the alpha=0 drag; the lift curve is a
physical-plausibility / py-vs-cpp check.

## What was run

Both `py/ibpm.py` and C++ `build/ibpm`, identical setup to the rest of
the suite (dx=0.02, nx=300, ny=150, domain length=6, `ngrid=1`, dt=0.01,
nsteps=3000 = t=30, Cl/Cd time-averaged over the last 60%). Driver:
[`../../run_naca0012_polar.py`](../../run_naca0012_polar.py); figures:
[`gen_naca0012_report.py`](gen_naca0012_report.py).

- **Re=500 polar**: alpha = 0, 2, 4, 6, 8, 10 deg -> `polar_comparison.png`
- **Re=1000, alpha=0**: a second independent drag anchor
- **Re=500, alpha=0 grid convergence**: dx = 0.04, 0.02, 0.01, 0.005 ->
  `grid_convergence.png` (does the immersed-boundary drag converge toward
  the benchmark as the grid refines?). Each level halves dx (and,
  starting at dx=0.02, halves dt to keep the run CFL-stable and rescales
  nsteps to hold t=30 fixed); driven by
  [`run_gridconv.py`](run_gridconv.py), which regenerates the matching
  resampled-boundary-point geometry for each new dx via
  [`../../make_airfoil_raw.py`](../../make_airfoil_raw.py) and skips any
  dx already present on disk.

## Results

**Python vs. C++: exact fidelity.** At Re=500 the flow is steady/laminar
and fully deterministic, so py/ibpm.py and C++ build/ibpm agree to machine
precision at every alpha (see `fidelity_summary.txt`) -- e.g. at alpha=0,
both give Cd=0.1891 to four+ significant figures, and Cl=-0.006 (~0, the
expected symmetry result). This is the cleanest kind of fidelity check:
no chaos to amplify differences (unlike the high-Re cases in
`../../airfoils/LSAT-SD7003/4-Re_sweep/`), so agreement is exact, not just
statistical.

**Vs. the CFD benchmark.** At dx=0.02, the immersed-boundary drag at
Re=500, alpha=0 is Cd~0.189, about 6-7% above the benchmark band
(0.176-0.178). That offset is the expected direction for an
immersed-boundary method at this resolution: the regularized (smeared)
boundary adds a small spurious drag that shrinks as dx -> 0. See
`grid_convergence.png`/`fidelity_summary.txt` for the exact per-dx
numbers; both implementations converge identically (Cd_py=Cd_cpp at
every dx tested, to 6 decimal places -- Re=500 is steady/deterministic,
so there's no chaos to cause disagreement here, same as the polar above).

**The convergence is not monotonic in step size, and that's worth
showing rather than hiding.** Cd(dx): 0.191858 (dx=0.04) -> 0.189095
(dx=0.02) -> 0.183971 (dx=0.01) -> 0.183510 (dx=0.005). The successive
differences (`|dCd|` per halving) are 0.002763, then 0.005124, then
0.000461 -- the step size nearly *doubled* between the first two
halvings before dropping to under a tenth of its previous size at the
third. A 3-point sequence (dx=0.04/0.02/0.01 only, as this study
originally had) can look like it's "converging" just because Cd is
monotonically decreasing, while actually still growing in *how much* it
moves each halving -- which is not yet asymptotic convergence, just a
still-changing trend that happens to be pointed the right direction.
Extending to dx=0.005 resolves this: the large dx=0.02->0.01 jump
turned out to be a pre-asymptotic transient (this resolution range is
where the immersed boundary's regularized delta-function support first
becomes narrow enough, relative to the body's curvature, to resolve
some feature it previously smeared over -- plausible but not
independently isolated here), and by dx=0.005 the sequence has entered
a genuinely flattening regime: the last step is small relative to Cd's
own magnitude, and the offset from the benchmark band has closed from
+6.2% (dx=0.02) to +3.1% (dx=0.005) -- essentially halved. This is
consistent with (not conclusive proof of) approaching a converged
Cd asymptotically as dx -> 0, which is what "resolution effect, not a
modeling error" requires as evidence.

**Why the sweep stops at dx=0.005, not dx=0.0025.** Wall-clock time per
grid level was measured directly (not estimated): dx=0.01 took ~6.7 min
per implementation; dx=0.005 (4x the cells, 2x the timesteps to hold
t=30 fixed at the smaller CFL-stable dt) took ~53 min -- an ~8x
increase, matching the 4x(cells) * 2x(steps) expectation almost exactly.
Extrapolating that same, empirically-consistent 8x-per-halving scaling,
dx=0.0025 (nx=2400, ny=1200) would cost roughly 7 hours per
implementation, ~14 hours total for both -- and by dx=0.005 the
turnover in step size described above is already large and unambiguous
(an 11x drop, not a marginal one), so a 5th point was judged not to
justify roughly half a day of additional compute. `run_gridconv.py`
already lists 0.0025 in its `GRID_DX`, so `python3
SURF_test/low_re/NACA0012/run_gridconv.py 0.0025` reproduces this
decision point exactly and extends the sequence further if wanted.

**Lift curve.** Cl(alpha) at Re=500 is smooth, near-linear at small alpha
with a reduced (sub-2*pi) slope characteristic of this low-Re regime, and
identical between the two implementations -- physically consistent with
the published low-Re NACA0012 lift behavior (Kurtulus, Di Ilio), even
though those references are figure-only for the quantitative values.

## Provenance

Coordinates: `naca0012.dat.txt`, converted from
`https://m-selig.ae.illinois.edu/ads/coord/n0012.dat` (Lednicer format ->
Selig closed loop, same conversion as `../../airfoils/LSAT-ClarkY/`). Reference Cd
values are transcribed from the papers cited in the table above (the
Re=500 trio via Nita et al.'s comparison table, arXiv:1901.08766; the
Re=1000 values via Di Ilio et al., arXiv:2006.10487, and Kurtulus 2015).
