# Reproducing Kurtulus (2019): pitching NACA0012 at Re=1000

Reproduces every test case in Kurtulus, D.F. (2019), "Unsteady aerodynamics
of a pitching NACA 0012 airfoil at low Reynolds number," *Int. J. Micro Air
Vehicles* 11:1-21, DOI [10.1177/1756829319890609](https://doi.org/10.1177/1756829319890609),
using both `py_static` and `cpp_static` (the fixed-algorithm builds from
[`../../../../../static_test/`](../../../../../static_test/)).

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

**Non-dimensionalization, explained term by term.** $c$ is the airfoil's
chord length (the straight-line distance from leading to trailing edge —
the paper's physical airfoil has $c=0.1$ m) and $U_\infty$ is the
free-stream speed the air approaches at ($U_\infty=0.146$ m/s in the
paper). Neither this solver nor most CFD codes work in physical units
directly — instead they solve in *non-dimensional* units where the chord is
rescaled to exactly 1 and the free-stream speed to exactly 1 ("solver uses
chord=1, U∞=1"). This is standard practice (it's what makes Reynolds number
comparisons meaningful in the first place) but it means any frequency given
in physical Hz has to be converted before it means anything to the solver:
the *non-dimensional* pitch frequency is $f_{nondim} = f \cdot c/U_\infty$
— literally "how many chord-lengths of flow pass by per pitch cycle."
Plugging in the paper's $c$ and $U_\infty$: 1 Hz → 0.684932, 4 Hz →
2.739726 (these are the frequency values actually passed to ibpm). This
same non-dimensional frequency, written as $k = 2\pi f c/U_\infty$ instead
(same idea, different normalization convention), is what the paper calls
the "reduced frequency" — $k=4.3$ and $k=17.2$ for 1 Hz and 4 Hz, quoted in
the "Test case" paragraph above as a cross-check that the conversion here
matches the paper's own numbers.

Separately, **pitch phase** controls *when in the cycle* the pitching
starts, and it's not just a bookkeeping detail — it can flip the sign of
the whole motion. The paper defines instantaneous angle of attack as
$\alpha(t) = \alpha_0 + A\sin(2\pi f t)$ (starts at the mean angle $\alpha_0$
and increases first). This solver's own sign convention for the pitching
motion, with phase left at its default of 0, would instead produce
$\alpha_0 - A\sin(2\pi f t)$ — the mirror image, decreasing first instead
of increasing. Setting phase $=\pi$ (180°) flips that back to match the
paper's convention (since $-\sin(x) = \sin(x+\pi)$). This was confirmed by
running one short pilot case and checking the sign of the resulting motion
directly, before committing to the full 258-run sweep — getting this wrong
wouldn't crash anything or produce obviously-bad numbers, it would just
silently compare against the paper's motion running backwards.

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
expected and is exactly what [`../../../../../static_test/`](../../../../../static_test/)
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
- **15°-40°**: the two solvers diverge structurally. **"Stall"** is the
  point where lift stops increasing smoothly with angle of attack because
  the flow can no longer stay attached to the airfoil's upper surface and
  separates instead (this is the same separation event that starts the
  vortex shedding described below); **"post-stall"** just means the angle
  range past that point, where the airfoil is operating in this
  separated-flow regime rather than the smooth, thin-airfoil-theory regime
  below stall. The paper's curves show a single sharp jump around 30°-35°
  (post-stall reattachment/regime change — lift recovering as the flow
  restructures itself further past stall), climbing cleanly to
  $\overline{C_l}\approx 2.0$. IBPM's curve in
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

**What "vortex shedding" and "Strouhal number" mean, for anyone new to
this**: past a certain angle of attack, the flow separates from the top of
the airfoil and the wake stops being a smooth, steady sheet — instead,
swirling vortices peel off alternately (first from one side, then the
other) at a regular rate, like a flag flapping. This alternating shedding
makes the lift and drag forces oscillate in time even when the airfoil
itself is held perfectly still ("steady" in this file's terminology means
the *airfoil* isn't moving — the *flow* around it can still be, and
usually is, unsteady). The **Strouhal number** ($St = f_{shed}\cdot
c/U_\infty$) is just that shedding rate, non-dimensionalized the same way
as the pitch frequency above — but here it's a property the *flow itself*
settles into, not something imposed from outside. It's measured, not set:
by watching how the lift coefficient oscillates over time and finding the
dominant frequency in that signal (via an FFT — see "Post-stall region"
below and `2-follow_up/README.md` for more on how, and where that
frequency-finding step itself can be a source of error).

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

**What this section is about**: everywhere else in this file, $C_l$/$C_d$
mean *time-averaged* lift/drag — one number per run, summarizing the whole
oscillating cycle. "Instantaneous" here means the opposite: the raw,
moment-by-moment force values as the airfoil actually pitches back and
forth, not yet averaged away. Plotting instantaneous $C_l$ or $C_d$ against
the *instantaneous* angle of attack $\alpha(t)$ (rather than against time)
produces a **hysteresis loop**: because the flow has "memory" (it responds
to how the airfoil got to its current angle, not just the angle itself),
the force on the way up ($\alpha$ increasing) differs from the force on the
way down ($\alpha$ decreasing) at that *same* instantaneous angle — so
instead of the up-stroke and down-stroke tracing the same line back and
forth, they trace two different paths that together form a loop. The shape
and width of that loop is itself physically meaningful (a wider loop means
a bigger difference between the up- and down-stroke response, i.e. a more
history-dependent, less quasi-steady flow), which is why the paper reports
it and why it's compared here.

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
- **Correction:** an earlier version of this section claimed the opposite
  of what the numbers below show — that $C_d$ was the strongest
  quantitative match and $C_l$ ran 20-25% high. Actually measuring it
  (`hysteresis_error_metrics.csv`, computed by interpolating the IBPM
  curve onto the paper's own tabulated α values, matched by up/down branch
  so the double-valued loop isn't averaged across both halves) shows the
  reverse: $C_d$ is the *worse* match by a wide margin. This was caught
  because the figure itself visibly shows $C_d$'s "tent" rising well above
  the paper's points at its peak, which is inconsistent with the old
  "closely matches" claim — the plot was right, the prose describing it
  was not.
- The **$C_d$-vs-α hysteresis loop** (`fig13_14_hysteresis.png`, right
  panel, now annotated with these numbers directly): IBPM's peak-to-peak
  $C_d$ oscillation is **2.23x the paper's (+123%)** — more than double —
  and the RMS error between the IBPM curve and the paper's 15 points is
  **42% of the paper's own $C_d$ range**. The shape (a symmetric "tent"
  peaking near α=0) matches qualitatively, but the amplitude does not.
- The **$C_l$-vs-α hysteresis loop** (left panel) has the correct shape and
  orientation (same tilted-ellipse sense, same phase relationship to
  instantaneous α), and is the *better*-matching of the two: IBPM's
  peak-to-peak $C_l$ oscillation is **1.31x the paper's (+31%)**, with RMS
  error **11% of the paper's own $C_l$ range** — consistent with (if
  somewhat larger than) the ~14% lift-slope excess noted above for the
  mean coefficients.

### Wake vorticity fields — `figures/wake_steady.png`, `figures/wake_f4hz.png`

Qualitative comparison at α₀=0°, 9°, 12° (0° = attached, 9° = right at the
paper's reported shedding onset, 12° = clearly post-onset), py_static vs.
cpp_static (indistinguishable, as expected), colored with a jet colormap
matching the paper's own figures (blue=negative, green≈0, red=positive).

**Mentor note, resolved: the airfoil-horizontal/flow-tilted look is a
plotting reference frame, not a solver difference.** `wake_steady.png` and
`wake_f4hz.png` plot the raw solver grid, in which the airfoil is horizontal
and the wake bends upward by ≈α₀. Kurtulus's own Fig. 2/Fig. 6 instead hold
the free-stream horizontal and pitch the airfoil. Tracing why: the solver
imposes α₀ by rotating the **free-stream**, not the body — `BaseFlow`
constructs `Flux.UniformFlow(grid, mag, alpha)` at angle α₀ relative to the
fixed grid axes, and `ibpm.py`'s force decomposition
(`drag = Fx·cos(a)+Fy·sin(a)`, `lift = -Fx·sin(a)+Fy·cos(a)`) confirms
`(cos α, sin α)` is the free-stream direction in grid coordinates. The
`.geom` file itself is never rotated per α₀ (`run_kurt_suite.py`'s
`ensure_geom` reuses the same raw NACA0012 geometry at every angle — this is
also why the airfoil silhouette in `wake_*.png` is horizontal in every row).
Kurtulus's figures use the opposite, more visually intuitive convention. Both
describe the exact same flow field, just related by a rigid rotation of the
whole picture by α₀ — not a physical discrepancy between the two ways of
imposing angle of attack.

To confirm this is *only* a plotting convention (and not, say, the solver
secretly producing an asymmetric or wrongly-rotated flow field),
`figures/wake_steady_paperframe.png` and `figures/wake_f4hz_paperframe.png`
re-plot the identical snapshots with `R(-α₀)` applied to the grid coordinates
and the airfoil outline before drawing (vorticity itself is a scalar and
doesn't transform). If the mismatch were a solver-level issue, rotating the
plot coordinates couldn't fix it — the wake would still fail to line up with
Kurtulus's layout. Instead, the rotated figures land on the paper's frame
almost exactly: free-stream horizontal, airfoil pitched nose-up by α₀, wake
undulating about the horizontal centerline rather than deflecting away from
it. This confirms the discrepancy is entirely the plotting frame; no solver
or physics difference is implicated. (The `*_paperframe.png` panels show a
white triangular gap at two corners — that's just the rotated square domain
no longer filling a rectangular axes bounding box, not missing data.)

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
- **Post-stall region (15°-40°, see "post-stall" defined above) is
  oscillatory/jagged in IBPM's mean coefficients**, unlike the paper's
  single clean jump. Plausibly the same phenomenon this repo has already
  characterized at higher Re elsewhere (`SURF_test/airfoils/README.md`'s
  "mentor question" investigation) — worth spelling out the connection,
  since it's not obvious at a glance: that investigation found that a
  Cartesian immersed-boundary solver like this one only resolves separated
  flow *cleanly* (smooth, non-chaotic) when the grid is fine enough,
  *relative to the boundary-layer thickness*, to actually capture the thin
  shear layer at the point of separation — not resolution or Reynolds
  number in isolation, but the two together (boundary-layer thickness
  itself shrinks as Re grows, so a fixed grid resolves it worse and worse
  as Re increases). That investigation worked at Re≈40,000-61,000, where
  under-resolution was severe enough to produce visibly broadband,
  grid-scale-speckled noise across the whole flow field. Re=1000 here is
  far more benign — nowhere near that speckled regime — but stall and
  post-stall angles are specifically where the flow separates and the
  boundary layer becomes a thin shear layer requiring finer resolution to
  capture cleanly than the smooth, attached flow at low angles does. So the
  same underlying mechanism (grid resolution relative to a shrinking
  length scale) could plausibly still be operating here in a much milder
  form — jagged, angle-to-angle-sensitive averages instead of the
  higher-Re case's outright visible speckle — even though Re=1000 is
  otherwise well-resolved. This is a plausible mechanism, not a
  conclusively isolated one — see `2-follow_up/README.md` for a direct test
  of whether this jaggedness is real or just a measurement artifact.
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
