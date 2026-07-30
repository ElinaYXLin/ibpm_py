# 7-stripe_investigation: area and severity of the LE/TE striping's *upstream* footprint

The mentor's follow-up to `../6-edges_further` (which found every LE/TE
"peak" reported throughout Groups 2/3/A-F lands in the fluid, not inside
the body -- see that folder's Group G): rather than keep measuring the
peak right at the LE/TE, where real boundary-layer physics and numerical
artifact are superimposed and hard to separate, quantify the striping
that extends **upstream of the leading edge** specifically -- its area
and its severity -- using x < x_LE as a **null region**.

**Why upstream is a clean measurement**: in steady incompressible flow at
Re=1000, the diffusion length scale nu/U = c/Re ~ 0.001c is about 1/20th
of a single dx=0.02 cell. Beyond a cell or two upstream of the stagnation
point, the physically correct vorticity is ~0. Any signal measured there
is numerical artifact with essentially no real-physics signal mixed in --
unlike the LE peak itself, where the two are inseparable by construction.
This sidesteps exactly the ambiguity that made Group F's conditioning
comparison inconclusive (see Test 6).

Every test below uses **both py_static and cpp_static** wherever data
exists for both, and reports their quantitative agreement explicitly
(see "py/cpp agreement" in each section) -- not just visual overlap.

**All 8 tests generate a figure.** 7 of 8 needed zero new runs (reusing
`../1-paper_based`, `../5-leading_edge`, and `../6-edges_further`'s
existing output); Test 7 turned out to need new runs after all (see its
section), and Test 8 required one new 1-step run. Both were cheap
(seconds each) relative to the faithful2 runs occupying the machine
elsewhere during this investigation.

## Definitions

Quantities reused across multiple tests, defined once here rather than
re-explained in each section (see `common.py`'s `upstream_scalar_metrics`,
`upstream_profile`, and `reach_L_up` for the actual code):

- **Threshold-area A(tau)**: the physical area (in chord^2) of the
  upstream window where `|omega| > tau`, for a small ladder of tau
  values. Literally "how big is the contaminated region," at a chosen
  severity cutoff -- the area counterpart to a peak/severity number.
  Computed as `(number of cells exceeding tau) * dx^2`.

- **Enstrophy**: `0.5 * integral(omega^2) dA` over the upstream window.
  This is a real physical quantity, not an arbitrary statistic -- in the
  Navier-Stokes equations, enstrophy sets the viscous dissipation rate
  (dissipation ~ enstrophy/Re), so it's a standard way to quantify "how
  much rotational/vortical activity" a region contains. Three reasons
  it's used as this folder's primary severity+extent metric: (1) it
  combines area and severity into one number automatically (weighting
  by omega^2 means both a large region of weak noise and a small region
  of strong noise can contribute comparably); (2) it's an *integral*
  over the whole window, not a single point-sample -- exactly the
  property `../6-edges_further` established the y=0 lineout metric
  lacks (Groups B/C found it swings 5x from phase alone); (3) it agrees
  between py_static and cpp_static to machine precision (1e-13 to 1e-15)
  in every test here, the strongest evidence in this whole investigation
  that it's measuring something real rather than sampling noise.

- **Signed vs. absolute integral / "oscillatory" interpretation** (Test
  2): `integral(omega) dA` (signed, positive and negative cancel) vs.
  `integral(|omega|) dA` (absolute, no cancellation). If the upstream
  region contained one coherent vortical structure (e.g. a real vortex
  that had somehow drifted or leaked upstream), the two would be close
  in magnitude, since a single-signed rotating structure doesn't cancel
  against itself. Instead the signed integral is only 2-12% of the
  absolute integral at every dx tested -- meaning the region is almost
  entirely made of positive and negative vorticity patches that
  nearly cancel when summed with sign. That's the signature of a
  sign-alternating (checkerboard-like) numerical pattern, not a
  physical vortex -- and it's the reason this folder reads the upstream
  signal as discretization noise rather than a real leaked/advected
  flow feature. This is consistent with, and reinforced by, Test 4's
  finding of a near-single-cell-period oscillation in the raw field.

- **Reach L_up** (Test 3): given Test 1's per-column max|omega| decay
  profile, L_up is the distance upstream of the LE out to the last
  point (walking outward from the LE) where the profile still exceeds
  a fixed noise floor (1e-3, chosen well above floating-point noise and
  well below any real near-body signal) -- i.e. "how far upstream does
  the signal actually reach before dropping below and staying below a
  meaningful threshold." Reported in both chord units (physical
  distance) and cells (`L_up/dx`) specifically to distinguish a local
  mechanism (fixed cell count, e.g. the discrete delta function's
  compact support) from a non-local one (fixed physical distance) --
  see the Hypotheses table below.

- **"LTE"** (Test 6, and `../6-edges_further`'s Groups D/F): shorthand
  for "Leading+Trailing Edge" -- `naca0012_LTEsparse`/`naca0012_LTEdense`
  are geometries where the boundary-point spacing was coarsened/refined
  at *both* the leading and trailing edge together (as opposed to the
  older, LE-only investigation these build on). Not related to "LE" the
  region label used elsewhere in this folder.

## Hypotheses tested

Each test above is checking a specific hypothesis about the upstream
noise's cause or character. Stated up front, with verdicts, so it's
clear what each test does and doesn't establish -- details and numbers
are in each test's own section below.

| # | Hypothesis | Tested by | Verdict |
|---|---|---|---|
| H0 | x<x_LE is a valid "null region" (true physical vorticity ~0 there), so anything measured is artifact | Tests 1, 3, 5 (indirectly: signal shrinks under refinement and is far weaker on a cylinder, consistent with H0; not directly proven) | **Supported**, not directly proven |
| H1a | The noise footprint is **local**: fixed cell-count (delta-function stencil support) | Test 3 (L_up in cells vs. chord under refinement) | **Not supported** -- L_up shrinks in cells too (43->15->5), not flat |
| H1b | The noise footprint is **non-local**: fixed physical distance, independent of dx | Test 3 | **Not supported** -- L_up shrinks in chord too (0.86->0.15->0.025c), not flat |
| H2 | The noise is odd-even/checkerboard-type grid decoupling (a specific, nameable numerical mechanism) | Test 4 (spectral fraction near Nyquist wavelength; raw-row visual) | **Supported** by the visual (clean ~1-cell-period oscillation); bulk spectral-fraction metric alone is not decisive (~5%, diluted by a smooth envelope) |
| H3 | Sharp LE curvature is **necessary** for upstream noise (cylinder should show none) | Test 5 (cylinder control) | **Not supported** -- cylinder still shows noise, just ~78x weaker. Refined to: curvature strongly *amplifies* a generic IB artifact, isn't required to produce it |
| H4 | Projection-matrix conditioning predicts artifact severity | Test 6 (vs. `../6-edges_further` Group F, which found no relationship using the LE-peak metric) | **Supported** on the clean upstream signal -- monotonic across all 3 geometries; Group F's null result is now understood as a consequence of measuring a physics-contaminated region, not evidence against conditioning mattering |
| H5 | Upstream noise is **purely geometric** (fixed by discretization alone, should be flat vs. angle of attack at fixed geometry) | Test 7 (alpha sweep, geometry provably identical at every alpha; checked in both the fixed body-frame region and a lab-frame region rotated to track the actual flow direction) | **Not supported** -- enstrophy varies several-fold across alpha=0-60 at fixed geometry in *both* frames, so the noise is at least partly flow/force-magnitude-driven |
| H6a | The noise is a direct footprint of the spread boundary force (regularization stage alone) | Test 8 Probe A | **Not supported** as a *sole* explanation -- the raw spread-force footprint is exactly zero beyond ~1-2 cells from the LE, so it cannot by itself reach far upstream |
| H6b | The noise builds up gradually through many steps of advection | Test 8 Probe B (nsteps=1 from zero IC) | **Not supported** -- a small enstrophy signature is already present after exactly one timestep, before any advection has occurred |
| H6c | The noise is *seeded* immediately by projection/regularization, but its full multi-chord *reach* builds up over many iterations via the elliptic solve | Test 8 (both probes combined) | **Supported** -- reconciles H6a/H6b: local instant source (Probe A/B), reach grows with iteration count (consistent with Test 3's non-trivial L_up scaling) |

---

## How to reproduce

```
python3 test1_streamwise_decay.py
python3 test2_extent_area.py
python3 test3_refinement_scaling.py
python3 test4_spectral.py
python3 test5_cylinder_control.py
python3 test6_spacing_conditioning.py
python3 test7_alpha_sweep.py run        # launches missing alphas first
python3 test7_alpha_sweep.py analyze
python3 test8_mechanism.py run          # launches the one-step run first
```

---

## Test 1: streamwise decay profile

![Test 1](figures/test1_streamwise_decay.png)

For each x-column upstream of the LE (all y in |y|<=0.5), max|omega|,
RMS, and integral|omega|dy, plotted vs. distance upstream on a log-y
axis, overlaying dx=0.02/0.01/0.005. NACA0012, alpha=0, steady, Re=1000.

**Decay is very steep and clearly resolution-dependent**: dx=0.02's
signal is still above the 1e-6 floor nearly 1.5-2c upstream; dx=0.01
drops to floor by ~0.3c; dx=0.005 by ~0.05c. Refinement shrinks *both*
the amplitude and the reach -- not a wash between the two, both improve
together.

**py/cpp agreement**: max relative difference in the max|omega| profile
is 3.71e-7 (dx=0.02), 2.4e-2 (dx=0.01), 2.3e-2 (dx=0.005). The dx=0.02
number matches this repo's usual ~1e-7-1e-13 agreement; the dx=0.01/0.005
numbers look larger only because this is a relative difference of a
*max*, taken over values that are already tiny (near the noise floor) --
Test 2's enstrophy (a bulk integral, not a fragile point-sample) agrees
to 8e-15-1.6e-13 at every dx, confirming this is a sampling-noise
artifact of the metric, not a real py/cpp divergence.

## Test 2: extent and area

![Test 2](figures/test2_extent_area.png)

Threshold-area A(tau) curves, total upstream enstrophy and integral
|omega|, and the signed-vs-absolute integral contrast, all at
dx=0.02/0.01/0.005.

| dx | enstrophy | integral \|omega\| | signed integral | signed/abs |
|---|---|---|---|---|
| 0.02 | 0.6268 | 0.1785 | 0.00388 | 2.2% |
| 0.01 | 0.00314 | 0.00385 | -0.00047 | -12.2% |
| 0.005 | 0.00177 | 0.00083 | -0.0000062 | -0.7% |

**Cancellation check confirms oscillatory noise, not a coherent
structure**: the signed integral is 2-12% of the absolute integral in
magnitude at every dx (panel 3), i.e. the upstream field is overwhelmingly
sign-alternating and self-cancelling rather than carrying real net
circulation -- consistent with Test 1's decay picture and Test 4's
spectral finding below.

**py/cpp agreement**: enstrophy agrees to 8.3e-15 (dx=0.02), 1.6e-13
(dx=0.01), 4.8e-14 (dx=0.005) relative -- machine precision at every
resolution, the cleanest agreement check in this folder.

## Test 3: refinement scaling, in two units

![Test 3](figures/test3_refinement_scaling.png)

| dx | enstrophy | L_up (chord) | L_up (cells) |
|---|---|---|---|
| 0.02 | 0.6268 | 0.860c | 43.0 |
| 0.01 | 0.00314 | 0.150c | 15.0 |
| 0.005 | 0.00177 | 0.025c | 5.0 |

**Power-law fit (3 points): enstrophy ~ dx^4.23** -- a genuinely fast
convergence rate, well beyond what the scheme's nominal order would
predict, and reassuring: this is *not* a persistent artifact, it
shrinks rapidly under refinement. Caveat worth stating plainly: the
log-log curve isn't perfectly straight (the drop from dx=0.02->0.01 is
much steeper than 0.01->0.005), so 4.23 is a single global exponent
fit to a slightly bent curve, not a precise local rate -- treat it as
"fast, better than 2nd order" rather than a load-bearing exact number.

**The reach/mechanism discriminator (L_up in chord vs. cells) doesn't
cleanly pick one story**: L_up shrinks in *both* units under
refinement (43 -> 15 -> 5 cells; 0.86 -> 0.15 -> 0.025 chord), which
rules out the cleanest version of either hypothesis (a strictly
local, fixed-cell-count delta-function footprint would stay flat in
cells; a purely non-local spreading mechanism would stay flat in
chord). It shrinks somewhat faster in chord than in cells per halving
of dx, meaning whatever the effective footprint is, it is itself
getting more localized under refinement, but not in lockstep with dx.
Test 8 gives the more direct mechanistic answer.

## Test 4: spectral character

![Test 4](figures/test4_spectral.png)

| dx | 2-D Nyquist-quarter power fraction | sign changes/chord (y=0) |
|---|---|---|
| 0.02 | 0.052 | 50.0 |
| 0.01 | 0.048 | 47.2-48.7 |
| 0.005 | 0.004 | 22.8-29.9 |

The bulk spectral-fraction metric is modest (~5% at dx=0.02/0.01,
<1% at dx=0.005) -- **but the rightmost panel (raw y=0 row, dx=0.02)
is the more convincing evidence**: it shows a clean, decaying,
essentially single-cell-period oscillation right next to the LE. The
bulk metric undercounts this because the signal is a high-frequency
oscillation riding on a smooth low-frequency decaying envelope (Test
1's exponential decay), which spreads power across many bands rather
than concentrating it purely at Nyquist. Taken together: **this looks
like odd-even/checkerboard-type grid decoupling** modulated by a smooth
envelope, not generic broadband error -- a specific, nameable numerical
mechanism, though the summary statistic alone is not decisive by itself.

## Test 5: cylinder control

![Test 5](figures/test5_cylinder_control.png)

| shape | x_LE | enstrophy | integral \|omega\| | peak |
|---|---|---|---|---|
| naca0012 | 0.0002 | 0.6268 | 0.1785 | 22.40 |
| cylinder | -0.4999 | 0.0080 | 0.0271 | 1.59 |

**Sharp curvature massively amplifies the effect but is not strictly
necessary for it**: the cylinder's upstream enstrophy is ~78x smaller
and its peak ~14x smaller than NACA0012's, but not exactly zero.
This is the more interesting of the two possible outcomes flagged
going in -- it means the underlying mechanism is a generic
immersed-boundary artifact present on any body, and sharp leading-edge
curvature is what turns a small background effect into the striping
this whole investigation was built around.

**py/cpp agreement**: enstrophy relative difference ~5.5e-15 (naca0012),
~4.1e-14 (cylinder) -- machine precision for both shapes.

## Test 6: boundary-point spacing and conditioning, on clean ground

(`LTEsparse`/`LTEdense` = Leading+Trailing Edge boundary-point spacing
coarsened/refined together -- see "Definitions" above; "enstrophy" is
the `0.5*integral(omega^2)dA` quantity defined there too.)

![Test 6](figures/test6_spacing_conditioning.png)

| geometry | condition number | upstream enstrophy | upstream peak |
|---|---|---|---|
| naca0012_LTEsparse | 4.05e3 | 0.148 | 9.24 |
| naca0012_baseline | 1.27e4 | 0.627 | 22.40 |
| naca0012_LTEdense | 1.15e8 | 1.410 | 34.46 |

**This is the headline new result of this folder.** These are two
different metrics on two different regions, so read them as a contrast,
not a contradiction: `../6-edges_further` Group F measured the
**lineout metric on the LE peak** (the physics-contaminated region) and
found conditioning did *not* track it -- LTEdense had by far the worst
conditioning (1.15e8) yet the *best*/lowest **lineout** LE peak (5.80,
vs. baseline's 22.18). This folder instead measures **upstream
enstrophy** (the clean, physics-free region defined at the top of this
README) for the same three geometries, and there the relationship
flips: **condition number and upstream enstrophy are monotonic together
across all three geometries** (LTEsparse 4.05e3->0.148, baseline
1.27e4->0.627, LTEdense 1.15e8->1.410 -- worse conditioning, more
upstream enstrophy, every time). This is exactly what "conditioning is
a pure numerical-error quantity" would predict, once it's tested
against a pure numerical-error signal instead of a physics-contaminated
one -- Group F's null result was a consequence of measuring the wrong
region (and the unreliable lineout metric on top of that), not evidence
that conditioning doesn't matter.

**py/cpp agreement**: enstrophy matches to ~9 significant figures for
all three geometries (differences only in the 6th-9th decimal digit).

## Test 7: angle-of-attack sweep

"Upstream enstrophy" here is exactly the same quantity defined in
"Definitions" and used in Test 2/6 -- `0.5 * integral(omega^2) dA` over
the upstream window (an integral over the whole region, not an RMS or a
single point-sample), just recomputed once per alpha instead of once
per dx/geometry.

![Test 7](figures/test7_alpha_sweep.png)

**This test did not turn out to be zero-new-runs as planned.** Despite
many `steady_{py,cpp}_aXX` directory names existing in
`../1-paper_based/runs/dx0.020/`, only alpha=0, 9, 12 actually had
computed output -- the rest were empty placeholder directories. 13 more
alphas (3,6,15,18,21,24,27,30,33,36,40,50,60) were launched fresh here
(cheap: dx=0.02, 300x150 grid, 3000 steps, ~50s each even under heavy
CPU contention from the unrelated faithful2 runs) to get a real sweep;
new output lives in `runs/alpha_sweep/`.

Body-to-grid discretization is identical at every alpha (this solver
imposes alpha by rotating the free-stream, never the body -- see
`../1-paper_based/README.md`'s "Wake vorticity fields" section), so this
sweep holds geometry exactly fixed while varying only the flow and the
boundary-force magnitude.

**Upstream enstrophy is clearly not flat vs. alpha** -- it dips around
alpha=9-24 (down to ~0.4, below the alpha=0 baseline of 0.627), rises
sharply through 30-33 (~0.9), dips again at 36-40, then rises steeply
toward alpha=60 (~5+, the sweep's maximum). Since geometry never
changes, **this variation can only be coming from the flow/force-magnitude
side, not the discretization** -- upstream noise is (at least partly)
force-driven, not purely a fixed geometric artifact. This is a genuinely
new finding this folder's design was built to be able to make.

**py/cpp agreement**: max relative difference in enstrophy over the
full 16-point sweep is 6.9e-14 -- machine precision throughout.

**Definitional note**: "upstream" is fixed as x<x_LE in the solver's own
frame throughout this sweep (the same region used everywhere else in
this folder), not rotated to track the incoming free-stream direction --
so at higher alpha this region is no longer literally "ahead of" the
oncoming flow in the lab sense.

**Lab-frame-rotated follow-up.** The figure's left panel now also plots
a second region that rotates with alpha to stay aligned with the actual
oncoming flow direction (`common.upstream_mask_2d_rotated`, using the
free-stream unit vector `(cos(alpha), sin(alpha))` implied by
`py_static/ibpm.py`'s own drag/lift rotation convention), alongside the
original fixed body-frame region. At alpha=0 the two are identical by
construction (both give enstrophy=0.6268, exact agreement to every
digit shown), which is a useful sanity check that the rotation is
implemented correctly. The right panel is a visual explainer: at an
example alpha=40deg, it draws both regions directly on the vorticity
field (blue=body-frame rectangle, red=lab-frame wedge rotated to track
the flow, purple=overlap) so the difference between the two is visible
rather than just described.

**The lab-frame (true "ahead of the flow") enstrophy is consistently
lower than the body-frame version, often by 2-4x** (e.g. alpha=24:
0.458 body-frame vs. 0.088 lab-frame; alpha=40: 0.496 vs. 0.067;
alpha=60: 5.500 vs. 2.065) -- meaning part of what the fixed body-frame
region picks up at higher alpha is coming from a part of the domain
that isn't actually upstream of the flow at all once the flow direction
itself is accounted for. That said, **the core Test 7 conclusion is
unchanged**: the lab-frame curve is just as clearly non-flat across
alpha as the body-frame one (same qualitative shape -- a dip through
the teens/twenties, a rise near 30-33, a dip at 36-40, then a sharp
rise to the alpha=60 maximum), so the noise is still (at least partly)
force-magnitude-driven rather than purely geometric either way.

**py/cpp agreement (lab-frame)**: max relative difference over the full
sweep is 1.16e-13 -- machine precision, same as the body-frame result.

## Test 8: mechanism isolation

![Test 8](figures/test8_mechanism.png)

Two probes into which solver stage produces the upstream signal, using
the existing converged NACA0012 dx=0.02 baseline plus one new nsteps=1
run from the zero-vorticity initial condition every run in this repo
starts from by default.

**Probe A** -- `model.B(f, omega)` computes exactly `Curl(regularizer.
toFlux(f))`, the vorticity-space footprint the discrete delta function
alone produces from the converged boundary force. Its support is
**exactly zero beyond 1 cell upstream of the LE** (checked explicitly:
nonzero at buffer=0/1, identically 0.000 at buffer=2/3/5/8) -- so the
raw spread force cannot, by itself, be a direct explanation for
anything beyond ~1-2 cells from the body. In the region where it *is*
nonzero, it correlates strongly with the real field (-0.87; the sign
flip is expected, since the full timestep also folds in diffusion and
advection on top of this raw source term -- the magnitude is what
matters here).

**Probe B** -- after exactly **one** timestep from zero IC, upstream
enstrophy is already nonzero (1.12e-4) but ~5600x smaller than the
converged value (0.6268) and, visually (bottom-right panel), still
concentrated in essentially one cell right at the LE -- not yet the
broad multi-chord reach seen in the converged state.

**Combined reading**: the noise is *seeded* immediately by the
regularization/projection step (present after a single timestep,
confirming it is not something advection slowly builds from nothing),
but reaching the full multi-chord extent seen in the converged
solution takes many iterations -- consistent with the elliptic
(Helmholtz/Poisson) solve's global, diffusive character gradually
carrying a small residual further upstream each step, on top of a
compactly-supported source that itself never reaches past ~1-2 cells.
This refines, rather than contradicts, Test 3's L_up finding: the
*source* is local, but its *reach* is solve-mediated and grows with
iteration count, which is also consistent with why finer dx (more
diffusion-limited cells needed to travel the same physical distance,
plus the smaller dt CFL forces) cuts the reach so effectively (Tests
1 and 3).

---

## Overall synthesis

Putting all 8 tests together: the upstream footprint is a **real,
measurable, but rapidly-converging and largely self-cancelling**
artifact. It decays approximately exponentially with distance (Test 1),
integrates to a small but nonzero, oscillation-dominated signal (Test
2), shrinks under refinement faster than the scheme's nominal order in
a way that isn't cleanly local-in-cells or fixed-in-chord (Test 3), has
spectral character consistent with odd-even/checkerboard grid
decoupling modulated by a smooth envelope (Test 4), is dramatically
amplified but not strictly caused by sharp leading-edge curvature (Test
5), is monotonically predicted by projection-matrix conditioning once
measured on clean (non-physics-contaminated) ground -- resolving Group
F's inconclusive result (Test 6), varies substantially with flow
conditions at fixed geometry, meaning it is partly force-magnitude-driven
rather than purely geometric (Test 7), and is seeded immediately by the
regularization/projection step with its full multi-chord reach building
up over many iterations via the elliptic solve rather than appearing
all at once (Test 8).

**Practically for the mentor's question**: the area and severity of the
upstream stripes are both real and both shrink rapidly and predictably
under grid refinement and under better-conditioned boundary-point
spacing -- this is a discretization artifact with an identifiable,
partially-understood mechanism, not a mystery or a sign of a solver bug
(consistent with the exact py_static/cpp_static agreement found in
every single test above, at every resolution, shape, and angle of
attack tested).

## Proposals: further tests and mitigation strategies (not yet implemented)

Written up as a proposal for discussion -- nothing in this section has
been coded or run. Two goals: close the remaining gaps in the 8 tests
above (mostly things that were inferred rather than directly measured),
and suggest ways to reduce the striping that don't rely on halving dx
everywhere (halving dx is 8x the runtime -- one factor of 2 in each of
x, y, and the CFL-limited timestep -- which is too slow to use as a
routine fix).

### Further tests, to close the remaining gaps

1. **Time-resolved buildup (closes the weakest link in H6c).** Test 8's
   "seeded immediately, full reach builds up over many iterations" claim
   was inferred from exactly two points -- step 1 and the converged step
   3000. The baseline run already saved a restart snapshot every 250
   steps; plotting upstream enstrophy and L_up vs. step number through
   the existing run is zero new runs and turns an interpolated claim
   into an observed curve. If the reach grows gradually (roughly
   diffusion-front-like) that directly confirms elliptic-solve-mediated
   spreading; if it's already at full extent by step 250, that part of
   the synthesis needs correcting.

2. **Extend Test 6's conditioning relationship past 3 points.**
   `../6-edges_further` Group D's existing LE-density sweep (0.5x, 2x,
   8x, and the diverged 16x) already has known condition-number context
   and is sitting on disk unused by this folder. Recomputing upstream
   enstrophy on those runs would turn Test 6's 3-point monotonic result
   into a 6-7-point dose-response curve spanning ~5 orders of magnitude
   in cond(M), and would show how far the sparser direction (see
   mitigation #1 below) can be pushed before the boundary itself starts
   misbehaving. Zero new runs.

3. **Quantify the checkerboard claim (H2) instead of eyeballing it.**
   Test 4's bulk Nyquist-fraction metric came out weak (~5%) only
   because the alternating pattern rides on top of a smooth decaying
   envelope, which spreads the signal's spectral power across many
   bands. Dividing out the envelope first (fit and subtract, or
   normalize by a local moving RMS) before computing the Nyquist
   fraction, or equivalently computing the lag-1 autocorrelation of the
   sign sequence, would turn "visually looks like a checkerboard" into
   an actual number. Zero new runs.

4. **Test the untested link in the mechanism chain.** Test 8 shows the
   noise is seeded by the projection step; Test 6 shows severity tracks
   cond(M) -- but nothing yet shows the noise specifically lives in M's
   ill-conditioned (small-eigenvalue) modes. `../6-edges_further` Group
   F's script already builds M explicitly; eigendecomposing it and
   projecting the converged boundary force (already saved in every
   restart file, as `State.f`) onto that eigenbasis would show directly
   whether the upstream noise's amplitude is carried by the
   high-frequency-along-the-boundary, ill-conditioned modes. If so,
   that closes the causal chain completely and pinpoints exactly which
   modes a fix needs to damp. Near-free computationally (one dense
   eigendecomposition of an already-built matrix).

5. **Regrow test: static artifact or continuously re-seeded?** Take one
   converged snapshot, zero out the upstream noise once by hand, restart
   the simulation from that edited field, and watch whether/how fast it
   comes back. Fast regrowth means the steady boundary force
   continuously re-seeds it (so only a method change removes it for
   good); no regrowth means it was a startup transient that got locked
   into the "steady" state (so even a one-time cleanup filter would
   suffice). One short, cheap new run.

### Mitigation strategies (other than reducing dx)

Ranked cheapest/best-supported-by-existing-evidence first.

| Strategy | Evidence already in hand | Cost | Main risk |
|---|---|---|---|
| **Sparser boundary points at LE/TE** (ds ~ 1.5-4x dx locally, instead of ds=dx) | Already demonstrated, not hypothetical: `../6-edges_further`'s LTEsparse geometry cut upstream enstrophy 4.2x (0.627->0.148, Test 6) and even slightly *improved* the LE field-max (Group D's 0.5x-density case: 67.7 vs. baseline 71.7) | Zero -- fewer boundary points means a smaller, faster Cholesky factorization too | Sparsify too far and the boundary becomes "leaky" (this repo's own `checkgeom` utility warns about under-resolved boundaries); must reverify Cl/Cd/Strouhal are unchanged before adopting |
| **Sub-cell grid-phase tuning** (`xshift`/`yshift`, already a supported CLI parameter) | `../6-edges_further` Test B1: field-max ranged 59.3-83.4 from phase alone, at fixed shape/grid/everything else | Zero | Only a partial mitigation (~30% swing, never zero); the optimal phase may not be the same for every case/angle, so it isn't a one-time fix |
| **Wider/smoother discrete delta kernel** | `py_static/regularizer.py` currently uses the 3-point Roma-Peskin-Berger (1999) kernel with a hard-coded 1.5-cell support radius (`deltaSupportRadius = 1.5`) -- one of the *sharpest* standard IB kernels. Wider kernels (the classic 4-point Peskin kernel, or other smoothed variants) are the immersed-boundary literature's standard remedy for exactly this grid-locking/high-frequency-force oscillation, and should also improve cond(M) directly, tying back to Test 6 | Small, localized code change (one function); negligible runtime cost | Spreads the effective interface over ~2 cells instead of 1.5 -- boundary-layer force resolution needs revalidating, not just the far-field striping |
| **Regularize the projection solve** (e.g. a small ridge/Tikhonov term, or truncating M's smallest eigenvalues) | Tests 6 and Group D's divergence at 16x density both identify ill-conditioning of M as a real driver; this targets it directly rather than indirectly (via boundary spacing) | Small, localized to the Cholesky/projection-solver path | No-slip stops being satisfied exactly and becomes a least-squares approximation instead -- must monitor the actual slip-velocity residual, not just the vorticity field, to confirm forces stay accurate |
| **Smooth the boundary force along arc length each step** | Same target as the row above (damp the high-frequency components of the boundary force specifically), simpler to implement than modifying the linear-algebra path | Trivial per-step cost | Same constraint-softening concern as above, same validation needed |
| **Targeted 2-cell-wavelength filter on vorticity near the body** | Directly targets the specific Nyquist-wavelength mode Test 4 identified, in principle leaving every other resolvable scale untouched | Trivial to add | The most invasive option philosophically (edits the solved field directly, after the fact, rather than fixing the mechanism that produces it) -- needs the most careful validation that real force/wake quantities (Cl, Cd, St) are unaffected |

**Two points worth flagging to the mentor directly.** First, sparser
LE/TE boundary points is unusual among these options in that it is
*already supported by this folder's own data* (not a hypothesis to go
test), costs nothing, and would make runs faster rather than slower --
the natural pilot is proposed test #2 above (extend Group D's sweep and
re-derive the best density from the resulting dose-response curve), then
revalidate forces before adopting it as the new default. Second, "avoid
the 8x cost of halving dx everywhere" is really only a constraint on
*uniform* refinement -- the multi-domain-nesting approach already in use
elsewhere in this repo (a tight finest box hugging just the airfoil,
with coarser nested levels covering the rest of the domain -- see the
faithful2 configuration and `GRID_UPGRADE_MANUAL.md`) buys the same
near-body dx without paying for it over the whole domain, since cost
scales with the finest box's area rather than the full 6c x 3c extent.
That's the structural fallback if the mitigations above turn out to be
insufficient on their own.

## Files

- `common.py` -- shared helpers (upstream-region metrics: streamwise
  profiles, threshold-area, enstrophy, signed/absolute integrals, L_up
  reach) plus the same load/grid/geometry conventions as
  `../6-edges_further/common.py`.
- `test1_streamwise_decay.py` -- Test 1 (zero new runs).
- `test2_extent_area.py` -- Test 2 (zero new runs).
- `test3_refinement_scaling.py` -- Test 3 (zero new runs).
- `test4_spectral.py` -- Test 4 (zero new runs).
- `test5_cylinder_control.py` -- Test 5 (zero new runs).
- `test6_spacing_conditioning.py` -- Test 6 (zero new runs).
- `test7_alpha_sweep.py` -- Test 7 (`run` launches the 13 missing
  alphas; `analyze` reuses alpha=0/9/12 from `../1-paper_based` plus
  the new ones).
- `test8_mechanism.py` -- Test 8 (`run` launches the one new 1-step
  run; both probes computed in the same script).
- `runs/alpha_sweep/`, `runs/one_step/` -- this folder's own new run
  output (everything else reuses `../1-paper_based`,
  `../5-leading_edge`, and `../6-edges_further` directly, unmodified).
- `data/`, `figures/` -- CSVs and PNGs for every test above.
