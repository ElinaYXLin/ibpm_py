# SD7003 Reynolds-number sweep (E1 + E2)

Pushes Re DOWN from SD7003's usual ~61,100 (`../2-c++included/`) toward
the cylinder's clean Re=100 baseline (`../../../vortall/1-baseline/`),
at fixed dx=0.02, alpha=4.60 (same geometry/domain/alpha as the existing
flowfield case). Originally run C++-only (this was framed as a
resolution/Re *behavior* question, not a port-fidelity recheck); **Python
was added afterward at every Re value** (see
[`../../run_airfoil_re_sweep_py.py`](../../run_airfoil_re_sweep_py.py))
so `re_sweep_comparison.png` now shows py/ibpm.py vs. C++ build/ibpm side
by side at every Re, and `fidelity_summary.txt` quantifies the agreement.
Figures from
[`../../gen_airfoil_re_sweep_figs.py`](../../gen_airfoil_re_sweep_figs.py).

## Fidelity result: exact agreement through Re=10000, then chaos-divergence (not a bug)

`fidelity_summary.txt`: domain-RMS vorticity matches to 4+ significant
figures (0.00-0.01% relative pointwise difference) between py/ibpm.py and
C++ build/ibpm at every Re from 200 through 10000 -- the two
implementations are, for practical purposes, computing identical
trajectories in this deterministic/laminar-to-transitional regime. At
Re=20000 and 40000 (both fully in the broadband-speckle regime), the
pointwise difference jumps to 12-143% even though domain-RMS/max stay in
the same ballpark -- this is the SAME chaos-amplification phenomenon
`../../../vortall/1-baseline/README.md` already documents for its own
statistically-converged Re=100 case ("match in periodic amplitude/
frequency, not instantaneous phase"): once the flow is genuinely chaotic,
floating-point-level differences between the two implementations'
operation order get amplified exponentially, so pointwise agreement is
not expected or meaningful there -- only statistical agreement is, and
that still holds (see `re_sweep_comparison.png`'s Re=20000/40000 rows:
same broadband character, different instantaneous micro-detail).

## Result: a clean, sharp transition -- and it answers the mentor's original question

`re_sweep_comparison.png` shows py/ibpm.py vs. C++ build/ibpm vorticity
fields side by side at t=30, for Re = 200, 500, 1000, 5000, 10000, 20000, 40000:

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
