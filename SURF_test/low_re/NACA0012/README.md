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
- **Re=500, alpha=0 grid convergence**: dx = 0.04, 0.02, 0.01 ->
  `grid_convergence.png` (does the immersed-boundary drag converge toward
  the benchmark as the grid refines?)

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
Re=500, alpha=0 is Cd~0.189, about 7% above the benchmark band
(0.176-0.178). That offset is the expected direction and magnitude for an
immersed-boundary method at this resolution: the regularized (smeared)
boundary adds a small spurious drag that shrinks as dx -> 0. See
`grid_convergence.png` -- refining dx moves Cd toward the reference band,
confirming the offset is a resolution effect, not a modeling error, and
that both implementations converge identically. (See the `fidelity_summary.txt`
for the exact per-dx numbers.)

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
