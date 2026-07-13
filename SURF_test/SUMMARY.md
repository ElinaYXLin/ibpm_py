# Summary: the "weird" vorticity field investigation, port fidelity, and genuinely low-Re results

This document ties together a multi-stage investigation into why
SD7003's vorticity field looked broadband/speckled instead of clean and
coherent, whether `py/ibpm.py` is a faithful port of the C++ reference
solver, and what happens at Reynolds numbers actually in the "hundreds"
range. Every claim below is backed by a script and a figure/data file
somewhere under `SURF_test/`; this is the map, not the primary evidence.

## 1. The original question

A mentor flagged that `SURF_test/airfoils/LSAT-SD7003/2-c++included/flow_evolution.png`
(Re=61,100) showed broadband, grid-scale vorticity speckle rather than
the clean, coherent vortex structures a textbook figure would show, and
suggested the cause might be Reynolds number -- IBPM-type immersed-
boundary solvers are supposed to suit low-to-medium Re well.

## 2. Ruling out candidate causes, one at a time

| # | Candidate cause | Where tested | Result |
|---|---|---|---|
| 1 | Wrong initial conditions | `airfoils/LSAT-SD7003/README.md`, `3-small_dt/` | Ruled out -- impulsive-start IC and Re=61,100 already match the UIUC LSAT dataset's own stated test condition |
| 2 | Far-field domain too small (`ngrid=1`) | `airfoils/LSAT-SD7003/3-ngrid=3/` | Ruled out -- the solver's own multi-domain scheme doesn't clean it up (and introduces a new instability past t~20) |
| 3 | Under-resolved (dx=0.02) | `airfoils/LSAT-SD7003/3-dx0.01/` | Ruled out -- finer dx doesn't clean it up; peak vorticity grows sharper without converging |
| 4 | The specific airfoil, or Re~61k specifically | `airfoils/LSAT-ClarkY/`, `airfoils/LSAT-GM15/` (Re=60,700 and Re=40,600) | Ruled out -- both still speckle identically; a 1.5x change in Re isn't enough |
| 5 | **Resolution relative to Re** (delta ~ c/sqrt(Re) vs. grid spacing dx) | `airfoils/LSAT-SD7003/4-Re_sweep/`, `airfoils/LSAT-SD8000/4-Re_sweep/`, `vortall/2-Re_sweep/`, `*/5-grid_refine/`, `vortall/3-grid_refine/` | **Confirmed** -- see below |

## 3. The answer: it *was* Reynolds number, just not the range tested

Two independent Reynolds-number sweeps, run in opposite directions on
two different geometries, converge on the same transition zone:

- **`airfoils/LSAT-SD7003/4-Re_sweep/`, `airfoils/LSAT-SD8000/4-Re_sweep/`**
  (Re swept DOWN from ~61k/60.8k to 200): clean, coherent wake through
  Re~200-1000, transitional waviness at Re~5000-10000, full broadband
  speckle by Re~20000-40000.
- **`vortall/2-Re_sweep/`** (Re swept UP from the cylinder's clean
  Re=100 baseline): clean through Re~1000, transitional at ~3000,
  speckled by ~10000 -- the same transition zone, found completely
  independently, on a completely different body shape (blunt cylinder
  vs. thin cambered airfoil).

**Grid-refinement at a fixed, transitional Re=5000**
(`airfoils/LSAT-{SD7003,SD8000}/5-grid_refine/`, `vortall/3-grid_refine/`)
confirmed the mechanism directly: the coarsest grid (dx=0.04) aliases a
genuine, organized shear-layer instability into broadband speckle;
dx=0.02 or finer resolves the *same physics* as ordered vortex structures
(visible Kelvin-Helmholtz-type roll-up), not noise. Domain-RMS vorticity
stays flat or drops as dx refines while peak vorticity grows sharper --
finer grids resolve real structure, they don't add spurious energy.

**Conclusion: the mentor's original intuition was correct.** This
solver's `ngrid=1`, dx=0.02, no-subgrid-dissipation configuration is
clean up to roughly Re~1000-3000 and broadband-speckled well above that.
SD7003/SD8000 (Re~61k/61k) and even ClarkY/GM15 (Re~41-61k) never came
close to testing low enough Re to see this -- section 5 below closes the
loop by testing Re actually in the hundreds.

## 4. Python-vs-C++ fidelity: exact agreement in the deterministic regime, expected chaos-divergence above it

Every experiment above was run in **both** `py/ibpm.py` and C++
`build/ibpm`, and every comparison figure (`re_sweep_comparison.png`,
`grid_refine_comparison.png` in each experiment folder) shows them side
by side, with `fidelity_summary.txt` quantifying the agreement:

| Regime | Re range tested | py-vs-cpp agreement |
|---|---|---|
| Laminar / transitional (deterministic) | 200 - 10,000 (airfoils); 100 - 3,000 (cylinder) | **0.00-0.01% relative pointwise difference** -- domain-RMS vorticity matches to 4+ significant figures; visually pixel-identical in every comparison figure |
| Broadband / chaotic | 20,000 - 40,000 (airfoils); 10,000 (cylinder) | 13-143% relative pointwise difference, but domain-RMS/max stay in the same ballpark |

The high-Re divergence is **not a fidelity defect**. It's the same
chaos-amplification phenomenon `vortall/1-baseline/README.md` already
documented for its own statistically-converged Re=100 shedding case: once
a flow is genuinely chaotic, two independently-computed trajectories
(even from bit-identical initial conditions, using two different
codebases with different operation orderings) diverge in *instantaneous
phase* exponentially fast, while remaining statistically equivalent (same
shedding period, same broadband character, same amplitude envelope). This
was already well-documented in `airfoils/LSAT-SD8000/2-c++included/port_fidelity_diagnostic.png`;
this investigation reconfirms it holds at every new Re/dx point tested,
including exact (0.00%) agreement at every grid-refinement level, and at
every Re up to the point flows actually turn chaotic.

**Conclusion: the Python port has verified, essentially perfect fidelity
to the C++ reference wherever the underlying physics is deterministic**
(the overwhelming majority of cases in this whole test suite), and
diverges from C++ only exactly where two independent runs of the *same*
C++ binary would also diverge from each other -- a property of chaotic
dynamics, not of the port.

## 5. Two airfoils genuinely in the "hundreds" Re range: `SURF_test/low_re/`

To directly answer the mentor's original suggestion (rather than relying
on the Re-sweep's endpoints alone), two airfoils were run from scratch at
Re=500 -- literally in the "hundreds" range:

- **NACA0012** ("easy": thin, symmetric, no camber) -- new to this repo,
  coordinates converted from UIUC's airfoil database.
- **SD7003** ("hard": cambered, laminar-separation-bubble-prone) -- the
  *same* airfoil that originally prompted the mentor's question, now at
  Re=500 instead of its usual ~61,100.

**No experimental wind-tunnel validation data exists at this Re** (see
`low_re/README.md` for the search): UIUC LSAT's practical floor is
Re~40,000-60,000, and the broader low-Re literature confirms Re=40k-60k
is close to a hard floor for reliable force-balance/momentum-method
measurement in a conventional wind tunnel. Genuinely low-Re (hundreds)
Cl/Cd data, where it exists at all, is CFD/DNS-computed, not
experimental, and not published in a machine-tabulated format. SD7003
(a non-canonical, less-studied low-Re airfoil) has no such computational
reference either, so it stays a Python-vs-C++ fidelity + qualitative
flow-field check; NACA0012 does have one, quantitatively -- see 5b below.

**Result**: both airfoils, at Re=500, show a smooth, coherent, laminar
wake sheet -- no speckle at any of the 7 snapshots from t=0 to t=30, in
both implementations, visually indistinguishable between py/ibpm.py and
C++ build/ibpm at every timestep (`low_re/NACA0012/flow_evolution.png`,
`low_re/SD7003/flow_evolution.png`). This directly confirms, on a
brand-new airfoil as well as the original one, exactly what
`airfoils/LSAT-SD7003/4-Re_sweep/` already found from a coarser Re grid.

## 5b. Quantitative low-Re validation against a NON-LSAT dataset: `low_re/NACA0012/`

The flow-field check above is qualitative. To also validate
*quantitatively* at Re in the hundreds, `low_re/NACA0012/` (this content
originally lived at `airfoils/Lockard-NACA0012/`, merged here since it's
the same airfoil/Re as the qualitative check above -- see that folder's
README for the full history) compares NACA0012 drag against published
**computational** benchmarks (the only kind that exist this low):

| Re | Cd(alpha=0) reference | Source |
|---|---|---|
| 500 | 0.1762 / 0.1759 / 0.178 | Lockard et al. / Wu et al. / Nita et al. (LBM) |
| 1000 | 0.119 / 0.119 / ~0.12 | Di Ilio et al. (HLBM) / (XFOIL) / Kurtulus |

Result: at dx=0.02, py and C++ agree to machine precision (Cd=0.1891 both,
at Re=500 alpha=0), ~7% above the benchmark band -- and grid-convergence
(dx=0.04 -> 0.02 -> 0.01) drives Cd monotonically *toward* the reference,
confirming the offset is the expected immersed-boundary resolution effect,
not a modeling error. This remains the suite's only non-LSAT
(computational-reference) validation.

## 6. Where everything lives

```
SURF_test/
  vortall/
    1-baseline/       cylinder Re=100 vs. VORTALL.mat (original validation)
    2-Re_sweep/        cylinder Re swept UP (100->10000), py vs. cpp
    3-grid_refine/     cylinder dx refined at fixed Re=5000, py vs. cpp
  airfoils/           (folders named <dataset>-<airfoil>)
    LSAT-SD7003/, LSAT-SD8000/   UIUC LSAT wind-tunnel validation (Re~61k/60.8k)
      1-orig/, 2-c++included/, 3-small_dt/, 3-ngrid=3/, 3-dx0.01/   (SD7003 only)
      4-Re_sweep/      Re swept DOWN (61k-> ... ->200), py vs. cpp
      5-grid_refine/   dx refined at fixed Re=5000, py vs. cpp
      6-explicit_dissipation/   (SD7003 only) post-hoc filter demo
    LSAT-ClarkY/, LSAT-GM15/     Re~41-61k, different-airfoil / lower-Re controls
  low_re/
    NACA0012/  genuinely low-Re (Re=500) airfoil: flow-field + py-vs-cpp
      check, PLUS quantitative Cl/Cd polar + grid convergence vs. published
      CFD drag benchmarks (formerly airfoils/Lockard-NACA0012/)
    SD7003/    genuinely low-Re (Re=500) airfoil, flow-field + py-vs-cpp check only
```

Each folder's own `README.md` has the full detail; `airfoils/README.md`
and `vortall/README.md` are the per-area index files this summary draws
from.
