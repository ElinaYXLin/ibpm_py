# Clark-Y airfoil vs. UIUC LSAT experiment

See [`../README.md`](../README.md) for why this directory exists (mentor
question about whether a lower-Re/different airfoil avoids the broadband
vorticity speckle documented for SD7003/SD8000) and full data provenance.

Clark-Y: flat-bottomed, historically the most widely-used airfoil in
aviation (dating to the 1920s). Picked as the "easy" case here -- simple
geometry, no laminar-separation-bubble design intent (unlike SD7003), and
tested by UIUC LSAT at Re≈60,700 -- deliberately the **same Re ballpark as
SD7003**, to isolate "does a completely different, far more mainstream
airfoil geometry avoid the speckle at the same Re SD7003 was run at?" as
its own question (see `../LSAT-GM15/` for the "lower Re, same-ish airfoil
complexity" question instead).

## `1-orig/`

The only content generated so far for this airfoil (mirrors what
`airfoils/LSAT-SD7003/2-c++included/` shows, but there's no earlier
Python-only `1-orig/` predecessor for this airfoil to distinguish from,
so this *is* the first pass here): `polar_comparison.png`, `drag_polar.png`,
`grid_convergence.png`, `flow_evolution.png`, `summary.txt` -- Python
`py/ibpm.py` vs. C++ `build/ibpm` vs. UIUC LSAT experiment, same style/code
as `airfoils/`'s `2-c++included/` (see
[`../../gen_clarky_gm15_report.py`](../../gen_clarky_gm15_report.py) /
[`../../gen_clarky_gm15_flowfield_figs.py`](../../gen_clarky_gm15_flowfield_figs.py)).

## Results

**Lift matches the experimental polar well**, same as SD7003/SD8000 (see
`1-orig/polar_comparison.png`, left panel). **Drag is overpredicted at
higher alpha** by both implementations, by a similar amount to each other
(e.g. at alpha=10.26°: py Cd=0.181, cpp Cd=0.194, both vs. exp Cd=0.035) --
same documented cause as SD7003 (no turbulence/transition model at this
Re). Python and C++ agree closely at every alpha and grid level (see
`1-orig/summary.txt`).

**The vorticity field (`1-orig/flow_evolution.png`) shows the same
broadband grid-scale speckle as SD7003/SD8000**, in both implementations
-- the headline (negative) result this pair of airfoils was run to check;
see `../README.md`'s "Results" section.

## One instability found and fixed: coarse-grid (dx=0.04) divergence

The grid-convergence sweep's coarsest level (dx=0.04, only 51 boundary
points around this thin airfoil's ~2.05c perimeter) **diverged to NaN
around t≈19.8-19.9 at dt=0.01, in both `py/ibpm.py` and C++
`build/ibpm`** -- the same rapid-blowup signature (force trace jumping
from O(1) to O(1e30+) within 2-3 steps) documented for SD7003's dx=0.01
fine grid and for the `ngrid=3` investigation in `airfoils/LSAT-SD7003/3-ngrid=3/`.
Fixed the same way: halving dt to 0.005 (doubling nsteps to 6000 to reach
the same t=30) made both implementations stable for the full run, with no
further changes. Likely cause: at only 51 points around this airfoil's
tighter leading-edge curvature than SD7003, the coarsest resolution's
boundary representation locally under-resolves the leading edge enough to
tighten the effective CFL limit past what dt=0.01 tolerates -- consistent
with the "leading-edge suction peak locally accelerates the flow" CFL
mechanism `airfoils/LSAT-SD7003/README.md` documents for its own dx=0.01 case,
here appearing at the *coarse* end instead because of Clark-Y's blunter
leading-edge geometry combined with very coarse point spacing. The
resulting coarse-grid Cl/Cd carries a large std. dev. (Cl std≈1.1, both
implementations) reflecting a genuinely much more unsteady solution at
this resolution -- consistent with, not contradicting, the grid-convergence
story SD7003/SD8000 already tell (unsteadiness shrinks sharply as dx
refines; see `1-orig/grid_convergence.png`).
