# Is SD8000's coarse-grid blow-up chaos, or a numerical instability?

Follow-up to [`../7-chaos_sensitivity/`](../7-chaos_sensitivity/README.md),
which showed that SD8000's coarse grid (dx=0.04) sometimes goes fully to
NaN within a few thousand steps, and that the exact blow-up step is
scrambled by imperceptible (1e-8 to 1e-5 relative) Re perturbations --
taken there as evidence of chaotic sensitivity ("butterfly effect").

That conclusion has a gap worth closing: **a NaN blow-up is not what
genuine chaos does.** Real chaotic turbulence lives on a *bounded*
attractor -- the solution stays finite forever, just unpredictable in
detail. Going fully to NaN is the textbook signature of a *numerical*
(discretization) instability. The two are easy to conflate because BOTH
make blow-up timing sensitive to tiny perturbations -- sensitivity alone
doesn't distinguish them. This directory runs three tests that do.

## Test 1: does refining dt (dx=0.04 fixed) delay the blow-up?

[`run_dt_refinement.py`](run_dt_refinement.py) reruns SD8000's coarse
case (dx=0.04, alpha=-0.81, Re=60800) at dt=0.01 (the original), 0.005,
and 0.0025, **5 repeats each** (C++, identical command every time -- the
repeats isolate genuine dt-dependence from `FFTW_EXHAUSTIVE`'s
run-to-run replanning noise, already documented in
`../7-chaos_sensitivity/README.md` as enough of a round-off-level
perturbation to reshuffle blow-up timing on its own). All runs held to
the same physical end time (t=40).

**If this were a CFL/temporal instability, smaller dt should delay or
remove the blow-up. It doesn't.**
[`dt_refinement_blowup.png`](figures/dt_refinement_blowup.png) /
[`data/dt_refinement_summary.txt`](data/dt_refinement_summary.txt): **all
15/15 runs blow up**, and the mean blow-up time is *not* monotonic in dt
-- dt=0.01: t=30.3±2.5, dt=0.005: t=25.2±1.8 (earlier, not later),
dt=0.0025: t=34.1±1.7. If anything the middle dt blows up soonest. This
rules out a straightforward temporal-stability (CFL-type) explanation:
refining the timestep neither delays nor removes the instability, it
just resamples where the (still noisy) blow-up time happens to land.

## Test 2: vorticity wavenumber spectrum just before blow-up

[`compute_spectrum_and_energy.py`](compute_spectrum_and_energy.py) reuses
field snapshots already saved by `../7-chaos_sensitivity/` (no new runs)
-- the last snapshot before each recorded blow-up (SD8000 py, SD8000
cpp, SD7003 py; SD7003 cpp never blew up through 4000 steps, so it's
excluded), 2D-FFTs the vorticity field, and radially bins it into a
power-spectral-density vs. wavenumber curve.

**Result:** [`spectrum_prebreakup.png`](figures/spectrum_prebreakup.png)
-- all three cases show the spectrum **turning up again as it approaches
the grid's Nyquist wavenumber** (dashed line), after first dipping at
intermediate k. A genuine (even if under-resolved) turbulent cascade
decays smoothly with k; energy piling back up right at the smallest
scale the grid can represent is the textbook signature of an **aliasing
instability** -- high-wavenumber content that has nowhere to go (this is
a DNS-style solver with no subgrid/dissipation model by design, see
`../LSAT-SD7003/6-explicit_dissipation/README.md`) folding back and
accumulating at the grid scale instead of cascading away.

## Test 3: does the growth to blow-up look bounded or runaway?

Same script, same reused snapshots: domain-integrated kinetic energy
(`0.5*InnerProduct(q,q)`) and enstrophy (`sum(omega^2)*dx^2`) tracked
from t=0 to just past each case's blow-up step.

**Result:** [`energy_growth.png`](figures/energy_growth.png) -- both KE
and enstrophy grow **monotonically and visibly accelerate** (on a log
scale) all the way to the NaN step, in all three cases, with no
saturation, oscillation, or plateau beforehand. A bounded chaotic
attractor would show fluctuating, statistically-stationary energy; this
is unmistakably a runaway.

## Test 4 (context, not rerun here): does refining dx remove it?

Not a new run -- this is exactly what `../run_all_airfoils.py` /
`../run_all_airfoils_cpp.py`'s existing grid-convergence sweep already
shows, and what `../LSAT-SD7003/README.md` / `../LSAT-SD8000/README.md`'s
grid-convergence figures document: dx=0.02 and dx=0.01 (both airfoils,
both implementations) run cleanly to completion with finite, physically
reasonable Cd/Cl -- **no blow-up at any finer resolution that's been
tried.** Only the coarsest level, dx=0.04, blows up. That contrast is the
other half of this test's logic: refining **dx** removes the instability
entirely; refining **dt** (test 1, above) does not touch it. That
combination -- spatially cured, not temporally cured -- is exactly what
an under-resolved *spatial* discretization/aliasing problem looks like,
and is not consistent with a temporal-stability limit.

## Conclusion

Combining all four: SD8000's (and SD7003's) dx=0.04 blow-up is better
explained as an **under-resolved spatial (aliasing) numerical
instability** than as genuine chaotic physics:

- it is cured by refining dx, not by refining dt (tests 1 and 4);
- its pre-blowup spectrum piles up at the grid's own Nyquist scale, the
  aliasing fingerprint, not a decaying physical cascade (test 2);
- its energy/enstrophy growth is a monotonic runaway to the NaN step,
  not bounded fluctuation around an attractor (test 3).

This refines (doesn't fully overturn) `../7-chaos_sensitivity/`'s
finding: the *sensitivity* documented there (imperceptible Re
perturbations scrambling blow-up timing, fresh reruns flipping which
implementation dies first) is real and reproduced again here (test 1's
own repeats scatter by several time units at fixed dt). But sensitivity
to perturbations is a property BOTH genuine chaos and a marginally
unstable numerical scheme share -- it was never enough on its own to
tell them apart, and tests 1-3 here point specifically at the
under-resolved-discretization explanation rather than bounded chaotic
turbulence. Practically: dx=0.04 is simply too coarse for this
configuration (consistent with `../5-grid_refine/`'s and
`../LSAT-SD7003/6-explicit_dissipation/`'s independent finding that the
same coarse grid produces broadband speckle elsewhere in this suite,
too) -- not a property of the underlying physics, and not a Python/C++
porting defect (both implementations show the identical pattern).

## Files

| File | What it is |
|---|---|
| [`run_dt_refinement.py`](run_dt_refinement.py) | Launches the 15 dt-refinement runs (3 dt x 5 repeats), parallel via `subprocess.Popen`. |
| [`gen_dt_refinement_figs.py`](gen_dt_refinement_figs.py) | Reads `_run_data/dt*_rep*/run.force`, writes `dt_refinement_blowup.png` + `data/dt_refinement_summary.txt`. |
| [`compute_spectrum_and_energy.py`](compute_spectrum_and_energy.py) | Reads `../7-chaos_sensitivity/_run_data/` (no rerun), writes `spectrum_prebreakup.png` + `energy_growth.png`. |
| `_run_data/dt<dt>_rep<NN>/` | Raw force traces for the 15 dt-refinement runs (`.bin`/`.cholesky` gitignored as elsewhere in this repo). |
| `figures/`, `data/` | Generated outputs listed above. |

Regenerate with:
```bash
python3 SURF_test/airfoils/8-dt_refinement_and_spectra/run_dt_refinement.py        # ~140s, needs build/ibpm
python3 SURF_test/airfoils/8-dt_refinement_and_spectra/gen_dt_refinement_figs.py
python3 SURF_test/airfoils/8-dt_refinement_and_spectra/compute_spectrum_and_energy.py
```
