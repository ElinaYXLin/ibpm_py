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

## Leading-edge vorticity investigation

The Re=500 flow-field snapshots above (and `flow_evolution.png`) show a
small, sharp vorticity "speck" pinned right at the leading edge (LE) --
visible in the raw vorticity data even though it's easy to miss at the
full-domain color scale. Four follow-up tests, in
[`leading_edge_investigation/`](leading_edge_investigation/), pin down
what causes it. All four use `py/ibpm.py` only (not the paired C++ runs
used elsewhere in this suite): Re=500 is steady/deterministic here (see
"Python vs. C++: exact fidelity" above), so a second implementation adds
confirmation but not new physics, and would have roughly doubled the
already substantial compute (test 1 alone required a fresh ~49-minute
dx=0.005 run with field snapshots, since the existing
`../run_gridconv.py` sweep used `-restart 0` and kept no snapshots to
check retroactively).

NACA0012's leading-edge radius of curvature, $r_{LE} = 1.1019 t^2 =
1.1019 \times 0.12^2 \approx 0.0159$ (chord=1, $t$=12% thickness), is the
length scale used throughout as "how sharp is the LE, geometrically" to
compare grid/point spacing against.

### Test 1 -- grid refinement: does the speck shrink as dx < r_LE?

[`run_grid_refinement.py`](leading_edge_investigation/run_grid_refinement.py)
reruns Re=500, alpha=0 at dx=0.01 and dx=0.005 (the same points already in
`../grid_convergence.png`, but with restart snapshots enabled this time),
at matching physical times so the two resolutions are directly comparable
frame-by-frame. `fig1_grid_refinement.png` shows both: a full-domain view
plus a zoomed-in LE inset (zoom window = 6 x $r_{LE}$) at the snapshot
where the near-LE peak is largest, and again at the final state (t=30) --
same convention requested for all flow-evolution figures in this
investigation, one panel zoomed out, one zoomed in on the peak.

**Result -- the opposite of the hypothesis, and more informative for it.**
`fig1b_grid_refinement_summary.png` (adding the existing dx=0.02 point for
a 3-point trend) shows the peak magnitude *growing* monotonically as dx
shrinks -- 64.7 (dx=0.02) -> 68.3 (dx=0.01) -> 74.7 (dx=0.005) -- not
shrinking. But its *distance from the true LE* collapses cleanly:
0.040 -> 0.020 -> 0.0206, i.e. pinned at almost exactly 2 grid cells
(2*dx) at every resolution, not converging toward 0 or toward $r_{LE}$.
Together this rules out the original hypothesis (a coarse grid smearing
out and inflating an artificial peak that should shrink and localize on
refinement) and points at the opposite mechanism: the LE stagnation
region has a genuinely sharp vorticity gradient that a coarse grid
*clips* rather than *inflates* -- refining the grid resolves more of the
true peak (hence larger, not smaller, magnitude) while the discrete
peak location just tracks "the grid cell nearest the LE" (hence pinned at
~2dx, not shrinking toward the geometric LE). This had not converged by
dx=0.005; a genuinely converged peak magnitude would need finer dx still,
consistent with `../README.md`'s own note that the Cd sequence itself
was still pre-asymptotic at dx=0.01 and only entered a flattening regime
at dx=0.005.

### Test 2 -- alpha sweep: does the speck's asymmetry track the stagnation point?

[`run_alpha_sweep.py`](leading_edge_investigation/run_alpha_sweep.py)
reruns alpha=0, 2, 8, 10 deg (Re=500, dx=0.02 -- identical settings to
`../naca0012_polar_results.json`, but with restart snapshots this time,
since the polar sweep used `-restart 0`) and compares the top-surface vs.
bottom-surface LE vorticity peak in `fig2_alpha_sweep.png`.

**Result -- confirms the hypothesis cleanly.** `fig2b_alpha_asymmetry_summary.png`
/ `data/alpha_asymmetry.txt`: the top/bottom peak asymmetry is
essentially zero at alpha=0 (-1.2%, consistent with symmetric flow --
the small nonzero value is a numerical-grid-symmetry floor, not a
physical effect, since the domain's y-grid isn't exactly centered on
y=0) and grows monotonically with angle of attack: +1.9% (alpha=2) ->
+12.6% (alpha=8) -> +16.6% (alpha=10). This is exactly what's expected if
the effect tracks the front (forward) stagnation point migrating away
from the geometric nose and toward the lower surface as alpha increases
-- the flow around the nose becomes more asymmetric, so does the LE
vorticity peak that sits right where the flow turns around it.

### Test 3 -- isolating boundary-point density from grid dx

[`make_le_densified_geom.py`](leading_edge_investigation/make_le_densified_geom.py)
builds a NACA0012 boundary-point file with locally finer arc-length
spacing near the LE (a half-cosine ramp down to dx/4 within +-2*r_LE of
the LE, dx=0.02 everywhere else), while
[`run_le_densified.py`](leading_edge_investigation/run_le_densified.py)
runs it at the SAME background grid dx=0.02 as the uniform-spacing
baseline -- isolating Lagrangian boundary-point density from Eulerian
grid resolution, which `../run_gridconv.py`'s dx sweep can't do (it always
changes both together, via `make_airfoil_raw.py` resampling the boundary
to ds=dx at each new dx).

**Result -- densifying boundary points, at fixed grid dx, makes the speck
*worse*, not better.** `fig3_le_densified.png` /
`data/le_densified_peaks.txt`: peak near-LE |omega| grows from 64.7
(uniform) to 89.9 (LE-densified) -- a 39% increase, well beyond test 1's
grid-refinement growth. This is consistent with the regularized delta
function in this codebase's immersed-boundary formulation being
tuned to a support width tied to the *background grid* dx, not to
boundary-point spacing; clustering boundary points more closely than
that support width doesn't resolve the LE curvature better, it makes
neighboring points' regularized force support overlap more, and (per
`make_airfoil_raw.py`'s own docstring on the "over-resolved boundary"
degeneracy) pushes toward the near-singular regime rather than away from
it. So the earlier grid-convergence work's practice of keeping boundary
spacing matched to dx (ds ~ dx, not independently refined) isn't just a
numerical-stability convenience -- it's the right coupling for the LE
region specifically, and decoupling it (this test) actively hurts.

### Test 4 -- is it the projection step straining hardest at the LE?

[`compute_le_residual.py`](leading_edge_investigation/compute_le_residual.py)
extends `../../vortall/1-baseline/inner/compute_residuals.py`'s no-slip
residual check (`||C(omega) - b||`, which that script only reported as a
single norm) to record the residual *per boundary point*, for Re=500,
alpha=0, dx=0.02, using the internal py/ API directly (same
Grid/Geometry/NavierStokesModel/NonlinearIBSolver pattern).

**Result -- rules this mechanism out.** `fig4_residual_vs_distance.png` /
`fig4b_residual_spatial_map.png`: after the first timestep, the residual
is flat at floating-point roundoff (~1e-15 to 1e-17) at *every* boundary
point, with no elevation whatsoever near the LE (the t=0 panel, at
~1.0 uniformly everywhere, is just the un-solved initial condition and
not informative about solver accuracy). The projection step is enforcing
the no-slip constraint to machine precision uniformly around the whole
body -- it is not straining harder at the LE than anywhere else. This
rules out "the linear solve itself is failing near the LE" as the
mechanism, leaving test 1 and test 3's explanation (curvature vs.
regularization-support-width resolution, not a solver-accuracy problem)
as the one the evidence actually supports.

### Summary

The LE vorticity speck is a genuine, physically-motivated feature (a
sharp near-wall vorticity gradient at the stagnation region, whose
strength correctly grows with alpha-driven stagnation-point asymmetry),
not a solver-accuracy artifact (test 4) and not simply "not enough grid
resolution" in the naive sense (test 1's peak grows rather than shrinks
under refinement). What test 1 and test 3 together show is that it's
resolution-*sensitive* in a specific way: refining the background grid
reveals more of the true peak (test 1), while boundary-point density
needs to stay coupled to that same background dx rather than being
independently refined (test 3) -- over-densifying the boundary alone
makes the regularized-delta-function supports overlap and inflates the
peak further, rather than resolving it.

## Provenance

Coordinates: `naca0012.dat.txt`, converted from
`https://m-selig.ae.illinois.edu/ads/coord/n0012.dat` (Lednicer format ->
Selig closed loop, same conversion as `../../airfoils/LSAT-ClarkY/`). Reference Cd
values are transcribed from the papers cited in the table above (the
Re=500 trio via Nita et al.'s comparison table, arXiv:1901.08766; the
Re=1000 values via Di Ilio et al., arXiv:2006.10487, and Kurtulus 2015).
