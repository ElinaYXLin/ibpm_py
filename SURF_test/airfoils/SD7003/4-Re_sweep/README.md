# SD7003 Reynolds-number sweep (E1 + E2)

Pushes Re DOWN from SD7003's usual ~61,100 (`../2-c++included/`) toward
the cylinder's clean Re=100 baseline (`../../../vortall/1-baseline/`),
at fixed dx=0.02, alpha=4.60 (same geometry/domain/alpha as the existing
flowfield case). C++ `build/ibpm` only -- see
[`../../run_airfoil_re_sweep.py`](../../run_airfoil_re_sweep.py)'s
docstring for why (this is a resolution/Re *behavior* question, not a
re-check of port fidelity, already established at the baseline Re
elsewhere in this suite). Figures from
[`../../gen_airfoil_re_sweep_figs.py`](../../gen_airfoil_re_sweep_figs.py).

## Result: a clean, sharp transition -- and it answers the mentor's original question

`re_sweep_comparison.png` shows both vorticity and velocity-magnitude
fields at t=30, for Re = 200, 500, 1000, 5000, 10000, 20000, 40000:

- **Re=200, 500, 1000: clean, coherent wake sheet.** No broadband
  speckle at all -- visually indistinguishable in character from
  `../../../vortall/1-baseline/`'s clean Re=100 cylinder result.
- **Re=5000, 10000: transitional.** Fine ripples/waviness appear along
  the wake sheet, growing with Re, but the flow is still largely
  organized (not yet broadband noise).
- **Re=20000, 40000: full broadband speckle**, matching what
  `../2-c++included/flow_evolution.png` (Re=61,100) already showed --
  this is where SD7003's usual runs sit.

`domain_rms_vs_Re.png` confirms this quantitatively: domain-RMS
vorticity grows slowly from Re=200-1000, then accelerates sharply past
Re~1000-5000 -- consistent with a delta/dx~1 threshold (delta ~
c/sqrt(Re) is the near-wall vorticity layer thickness; at dx=0.02,
delta=dx predicts Re~(c/dx)^2 = 2500, right in the observed transition
zone).

**This directly answers the mentor's original suggestion** ("switch to a
low-Re airfoil") -- it was right in spirit, just needed a much lower Re
than ClarkY/GM15 (Re=40,600-60,700) tested: genuinely low Re (a few
hundred to ~1000) DOES give a clean result on this exact airfoil, exact
solver, exact `ngrid=1` domain that looked "weird" at Re~61k.
[`../../../vortall/2-Re_sweep/`](../../../vortall/2-Re_sweep/) finds the
same transition zone from the opposite direction (pushing the clean
cylinder baseline's Re UP), independently confirming the same Re~1000-5000
threshold on a completely different geometry.
