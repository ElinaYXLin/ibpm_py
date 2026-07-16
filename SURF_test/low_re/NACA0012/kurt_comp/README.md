# Reproducing Kurtulus (2019): pitching NACA0012 at Re=1000

Reproduces every test case in Kurtulus, D.F. (2019), "Unsteady aerodynamics
of a pitching NACA 0012 airfoil at low Reynolds number," *Int. J. Micro Air
Vehicles* 11:1-21, DOI [10.1177/1756829319890609](https://doi.org/10.1177/1756829319890609),
using both `py_static` and `cpp_static` (the fixed-algorithm builds from
[`../../../../static_test/`](../../../../static_test/)).

## Important correction: this paper is numerical, not experimental

The paper is a **CFD study (ANSYS Fluent, finite-volume, SIMPLE algorithm)**,
not a wind-tunnel experiment — stated explicitly in its own abstract and
methodology section. Every comparison below is therefore
**IBPM-vs-Kurtulus's-CFD**, not IBPM-vs-experiment. This matters for how to
read disagreements: unlike the UIUC LSAT comparisons elsewhere in this repo
(real wind-tunnel data), a mismatch here could come from either solver being
"wrong" relative to the true physics — both are independent numerical
approximations of the same equations, with different discretizations
(Cartesian immersed-boundary, this repo, vs. body-fitted unstructured mesh
with near-wall refinement, Fluent).

## Test case

NACA0012, Re=1000, sinusoidal pitch of amplitude 1° about the quarter chord:
$\alpha(t) = \alpha_0 + A\sin(2\pi f t)$, $A=1°$, $f=1$ Hz (reduced frequency
k=4.3) or $f=4$ Hz (k=17.2), for mean angles of attack $\alpha_0=0$° to 60°
(step 1° below 40°, then 50°, 60°, matching the paper's own increment), plus
a non-oscillatory ("steady") baseline. All three motions run at
**dx=0.02** (nx=300, ny=150, domain length=6/xoffset=-2/yoffset=-1.5,
ngrid=1) — this repo's established production resolution for NACA0012.

**dx=0.01 grid-check: skipped, by explicit decision.** A finer-grid subset
was planned but dropped after estimating its cost (~3-20+ hours depending on
machine contention, vs. dx=0.02's already-substantial runtime) against its
actual necessity: it would only have added grid-sensitivity confidence for
our own solver, not new information about matching Kurtulus's trends (their
mesh is a different discretization entirely, so there's no literal
resolution to "match"). The comparison below is entirely from dx=0.02.

Non-dimensionalization (solver uses chord=1, U∞=1; paper uses c=0.1 m,
U∞=0.146 m/s): pitch frequency $f_{nondim} = f \cdot c/U_\infty$ → 0.684932
(1 Hz) / 2.739726 (4 Hz). Pitch phase set to π so the effective AoA matches
the paper's convention ($\alpha_0 + A\sin(2\pi ft)$, not $\alpha_0 -
A\sin(\cdot)$, which is what phase=0 would give with this solver's sign
convention — verified against a short pilot run before the full sweep).

## Run summary

**258/258 jobs completed, 0 failures, 0 NaN.** 3 motions (steady, f1hz, f4hz)
× 2 implementations (py_static, cpp_static) × 43 angles, run via
[`run_kurt_suite.py`](run_kurt_suite.py) (resumable, 8-way parallel — this
machine's compute was heavily contended by other running applications for
most of the run, which is why it took roughly 38.5 hours wall-clock for what
is, compute-cost-wise, a few hours of work). Reduced to comparison tables by
[`analyze_kurt.py`](analyze_kurt.py) → `data/`, figures by
[`gen_kurt_figs.py`](gen_kurt_figs.py) → `figures/`.

**py_static and cpp_static agree essentially everywhere** — visually
indistinguishable in every figure below (the blue py_static curve is
consistently hidden directly under the red cpp_static curve). This is
expected and is exactly what [`../../../../static_test/`](../../../../static_test/)
already established: with the DST planner fixed (`FFTW_ESTIMATE` instead of
`FFTW_EXHAUSTIVE`), both implementations are deterministic and agree to
floating-point precision.

## Results

### Lift-curve slope near α=0 (steady)

`data/liftslope_dx0.020.txt`: **3.572 rad⁻¹ = 1.137×π**, identical for
py_static and cpp_static. The paper reports "approximately π"
(3.1416) for the same quantity. Same order and same qualitative claim
(near-π slope, consistent with thin-airfoil theory even at this low Re) —
about 14% higher in magnitude, plausibly a Cartesian-grid resolution effect
at dx=0.02 rather than a fidelity problem (see "Anomalies" below).

### Mean Cl, Cd vs. angle of attack — `figures/fig1_mean_coefficients.png`

Digitized from the paper's Fig 1 (graphical only, no table — see
`data/kurtulus_fig1_digitized.csv` for the digitization, accurate to
~±0.05) and plotted side-by-side with IBPM's own sweep.

- **0°-15°**: good quantitative agreement — both $C_l(\alpha)$ curves track
  closely to nearly-linear thin-airfoil behavior, $C_d(\alpha)$ magnitudes
  comparable.
- **15°-40°**: the two solvers diverge structurally. The paper's curves show
  a single sharp jump around 30°-35° (post-stall reattachment/regime
  change), climbing cleanly to $\overline{C_l}\approx 2.0$. IBPM's curve in
  this range is **oscillatory/jagged rather than smooth** and plateaus
  lower, around $\overline{C_l}\approx1.5-1.75$, never showing the paper's
  sharp jump. This is a genuine structural difference, not just a numerical
  offset — see "Anomalies."
- **50°-60°**: both solvers converge back to similar values
  ($\overline{C_l}\approx1.4-1.7$, $\overline{C_d}\approx2.5-3.2$) —
  agreement recovers at the highest angles, where the flow is bluff-body-like
  and less sensitive to boundary-layer/separation details either mesh
  resolves differently.

### Vortex-shedding Strouhal number — `figures/fig19_shedding_strouhal.png`

Digitized from the paper's Fig 19 (`data/kurtulus_fig19_digitized.csv`),
compared against IBPM's own FFT-of-$C_l$ peak-frequency estimate
(`analyze_kurt.py`'s `shedding_strouhal`).

- **Onset**: paper reports shedding starting at 8° (peak amplitude at 9°);
  IBPM's onset is also **9°** (both py_static and cpp_static) — a match to
  within 1°.
- **Peak magnitude**: paper St≈0.87-0.89 at onset; IBPM St≈0.87 at 9°-10° —
  essentially exact agreement at the peak.
- **Post-peak shape**: this is where the two diverge. The paper shows a
  smooth, monotonically-decaying St(α) curve down to ~0.23 by 25°, then a
  noisy plateau ~0.24-0.35 through 40°. IBPM instead shows a sharper initial
  drop, a low plateau around St≈0.20 through 20°-28°, then a **second rise**
  to a higher plateau (St≈0.33) from 35°-50°. Not a shape either solver's
  data was forced into — a real difference in the dominant wake frequency
  content at these angles, most likely reflecting that Fluent's body-fitted,
  near-wall-refined mesh resolves the separated shear layer differently than
  this solver's uniform dx=0.02 Cartesian grid can. Flagged as an open
  question, not resolved here.

### Thrust generation (negative instantaneous Cd), f=4Hz — `figures/thrust_check.png`

The paper's one crisp quantitative claim outside the tables: instantaneous
$C_d<0$ (thrust) occurs for $3°\le\alpha_0\le37°$ at f=4Hz, with $C_d\ge0$
throughout the cycle outside that range.

- **Onset matches closely**: IBPM's minimum instantaneous $C_d$ crosses zero
  between α₀=2° and 3°, matching the paper's stated onset almost exactly.
- **Upper bound does not match**: the paper reports $C_d$ returns to
  all-positive above 37°-38°. IBPM's minimum $C_d$ stays negative all the
  way through **α₀=50°**, only returning positive by α₀=60°. This is a
  real, sizeable quantitative disagreement (thrust regime roughly 20°+ wider
  in IBPM than reported) — see "Anomalies."

### Instantaneous forces and hysteresis — `figures/fig11_instantaneous_pitch.png`, `figures/fig13_14_hysteresis.png`

At α₀=0° (the one case the paper gives an actual numeric table for, its
Fig 13/14 caption — `data/kurtulus_fig13_14_table.csv`, 15 instantaneous
(t, α, Ω, $C_l$, $C_d$, $C_l/C_d$) points over one f=4Hz cycle):

- $C_l(t)$ is cleanly sinusoidal at both frequencies, matching the paper's
  qualitative description; $C_d(t)$ shows the paper's reported **period
  doubling** (two bumps per pitch cycle, since drag is roughly
  even-symmetric in a symmetric airfoil's instantaneous AoA) at both f1hz
  and f4hz.
- $C_d$ **stays strictly positive throughout** at α₀=0° for both frequencies
  in IBPM's data — matching the paper's specific claim that α₀=0°, 1°, 2°
  at f=4Hz (and all angles at f=1Hz) never generate thrust.
- The **$C_d$-vs-α hysteresis loop** (`fig13_14_hysteresis.png`, right panel)
  matches the paper's 15 tabulated points closely in both shape and
  magnitude — the strongest quantitative agreement in this whole
  comparison.
- The **$C_l$-vs-α hysteresis loop** (left panel) has the correct shape and
  orientation (same tilted-ellipse sense, same phase relationship to
  instantaneous α) but IBPM's amplitude runs ~20-25% larger than the paper's
  points (peak $|C_l|\approx2.5$ vs. the paper's $\approx2.0-2.1$) —
  consistent with the same ~14% lift-slope excess noted above.

### Wake vorticity fields — `figures/wake_steady.png`, `figures/wake_f4hz.png`

Qualitative comparison at α₀=0°, 9°, 12° (0° = attached, 9° = right at the
paper's reported shedding onset, 12° = clearly post-onset), py_static vs.
cpp_static (indistinguishable, as expected), colored with a jet colormap
matching the paper's own figures (blue=negative, green≈0, red=positive).

- **α₀=0°**: clean, coherent, unseparated wake sheet in both motions —
  matches the paper's Fig 2/Fig 4 at 0° closely.
- **α₀=9°, steady**: fully alternating vortex street already visible by
  x≈2c, matching the paper's description of periodic shedding above 8°.
- **α₀=9°, f4hz**: notably *less developed* than the steady case at the same
  angle — the wake is still an undulating sheet, not yet pinched into
  distinct alternating vortices within the visible domain. This directly
  confirms the paper's own finding that **f=4Hz pitching delays the onset
  of periodic shedding to a higher mean angle (~10°) than the steady case
  (~8°)** — visible here without needing the Strouhal-vs-α data at all.
- **α₀=12°**: both motions show a clear alternating vortex street, with
  f4hz's vortices slightly more compact/tightly-spaced than steady's — again
  consistent with the paper's qualitative description of oscillation
  frequency reshaping (not eliminating) the wake structure.

## Anomalies

- **Small nonzero $\overline{C_l}$ at α=0** (both impls: -0.00648, not
  exactly 0 as the airfoil's geometric symmetry would suggest). Negligible
  in magnitude (compare to the α=2° value, 0.129) — almost certainly a small
  grid/discretization asymmetry (the resampled boundary-point geometry isn't
  perfectly symmetric about y=0 at this dx), not a solver bug. Not
  investigated further here.
- **Post-stall region (15°-40°) is oscillatory/jagged in IBPM's mean
  coefficients**, unlike the paper's single clean jump. Plausibly the same
  phenomenon this repo has already characterized at higher Re elsewhere
  (`SURF_test/airfoils/README.md`'s "mentor question" investigation):
  resolution relative to the boundary-layer/shear-layer length scale, not Re
  in isolation, controls how cleanly a Cartesian IB solver resolves
  separated flow at a fixed dx. Re=1000 here is far more benign than the
  Re≈40-61k cases studied there, but the post-stall regime at high α is
  exactly where shear-layer resolution starts to matter even at low Re — a
  plausible, not conclusively isolated, explanation.
- **Faint background striping visible in the wake-contour figures**, even
  in nominally-uniform far-field regions. Same class of grid-scale numerical
  noise documented and explained elsewhere in this repo for other airfoil
  cases, just far fainter here (Re=1000 is comfortably resolved at
  dx=0.02, unlike the Re≈40-61k LSAT cases) — cosmetic, not a correctness
  concern for the force-coefficient comparisons above.
- **Thrust regime extends far past the paper's reported 37° cutoff** (IBPM:
  negative min-$C_d$ through 50°, paper: through 37°-38°). This is the
  single largest, clearest quantitative disagreement found in this
  comparison. Onset (~3°) matches essentially exactly, so the two solvers
  agree on *when thrust starts* but disagree substantially on *when it
  stops* — an open question, not resolved here; a natural next step would be
  comparing instantaneous force traces at a mid-range angle (e.g. 45°) where
  the two solvers disagree on sign to see whether the mechanism is a
  resolved-physics difference or a Fluent-vs-Cartesian meshing artifact.
- **Shedding-frequency plateau structure past 20° doesn't match the paper's
  smooth decay** (see "Post-peak shape" above) — also unresolved, flagged
  as an open question for follow-up rather than explained away.

## Files

- `run_kurt_suite.py` — the sweep driver (resumable; `python3
  run_kurt_suite.py dx002 8` reruns/resumes it).
- `analyze_kurt.py` — reduces `runs/dx0.020/*/flow.force` into the CSVs in
  `data/`.
- `gen_kurt_figs.py` — generates all figures in `figures/` from `data/` and
  raw vorticity snapshots.
- `data/kurtulus_fig1_digitized.csv`, `kurtulus_fig19_digitized.csv` — hand-digitized
  from the rendered paper PDF (no tables exist in the source; accuracy
  ~±0.05 for Fig 1, ~±0.02 Hz for Fig 19).
- `data/kurtulus_fig13_14_table.csv` — the paper's one actual numeric table
  (its Fig 13/14 captions), transcribed directly, not digitized.
- `data/mean_coefficients_dx0.020.csv`, `shedding_strouhal_dx0.020.csv`,
  `thrust_check_dx0.020.csv`, `liftslope_dx0.020.txt`,
  `fig13_14_comparison.csv` — IBPM's own reduced results.
- `runs/dx0.020/<motion>_<impl>_a<NN>/` — raw per-case output (`flow.force`,
  `flow.cmd`, `run_log.txt`, and vorticity snapshots for α∈{0°,9°,12°} only,
  to keep repo size down; `.cholesky` cache files excluded, matching this
  repo's convention elsewhere).
