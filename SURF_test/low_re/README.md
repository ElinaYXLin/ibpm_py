# Low-Reynolds-number airfoil validation: ClarkY and GM15

## Why this directory exists

A mentor asked, after the broadband grid-scale vorticity speckle
documented throughout `high_re/SD7003/README.md` and `high_re/SD8000/README.md`
didn't clear up under either `ngrid=3` (multi-domain far field) or `dx=0.01`
(finer resolution) -- see `high_re/SD7003/3-ngrid=3/` and
`high_re/SD7003/3-dx0.01/` -- whether switching to a genuinely
low-Reynolds-number airfoil to start with would help, since IBPM-type
immersed-boundary solvers are supposed to be well-suited to low-to-medium Re.

SD7003/SD8000 were tested at Re≈61,100/60,800 -- already commonly called
"low Reynolds number" in the aerodynamics literature, but this directory
asks the question two ways at once (see the two folders below), so the
airfoil-identity and Re-magnitude variables aren't conflated into a single
change.

## The two airfoils

Both come from the same UIUC LSAT (Low-Speed Airfoil Tests) source and
clean-tabulated `.DRG`/`.LFT` text format as SD7003/SD8000 -- see
"Data provenance" below.

- **[`ClarkY/`](ClarkY/)** -- the Clark-Y airfoil, Re≈60,700 (**same Re
  ballpark as SD7003**, deliberately -- this isolates "does a different,
  far more mainstream airfoil geometry avoid the speckle at the same Re?"
  as its own question). Clark-Y is a flat-bottomed, historically the most
  widely-used airfoil in aviation (dating to the 1920s, used on countless
  early aircraft designs and still a standard reference/teaching airfoil
  today) -- the "easy" case: simple geometry, generally attached/benign
  flow behavior relative to a laminar-separation-bubble-prone design like
  SD7003.
- **[`GM15/`](GM15/)** -- the GM15 airfoil (Gilbert Morris, F1C-class free
  flight design), Re≈40,600 -- **genuinely, substantially lower than
  SD7003's 61,100** (roughly 1.5x lower), isolating "does meaningfully
  lower Re avoid the speckle?" as its own question. GM15 is a small-model
  airfoil, not a famous name -- see "Why not a more famous airfoil at this
  Re" below for why a widely-used low-Re airfoil with a genuinely lower
  clean experimental Re bin than SD7003 turned out to be hard to find.
  Nonlinear, more strongly cambered Cl-alpha behavior at this Re than
  Clark-Y makes it the "hard" case of this pair.

## Results: same speckle, at both airfoils, at both Re levels

**The broadband grid-scale vorticity speckle documented for SD7003/SD8000
appears identically in both `ClarkY/1-orig/flow_evolution.png` and
`GM15/1-orig/flow_evolution.png`**, in both `py/ibpm.py` and C++
`build/ibpm`, filling the same far-field regions that should be
undisturbed uniform flow. Neither a completely different, much
more benign, historically ubiquitous airfoil (Clark-Y) at SD7003's own Re,
nor a genuinely lower Re (GM15, Re=40,600) cleared it up.

This is consistent with (not a new finding contradicting) the mechanism
already documented in `high_re/SD7003/README.md`'s "Limitations" section
and `high_re/SD7003/3-ngrid=3/summary.txt`'s diagnostic: the speckle comes
from running this solver's `ngrid=1` single-uniform-grid configuration
over a 6-chord domain with no explicit subgrid dissipation -- a property
of the **domain/resolution configuration shared by every case in this
whole test suite**, not of any particular airfoil's geometry or Reynolds
number. Changing the airfoil or lowering Re by 1.5x doesn't touch either
of those two variables.

**The integrated force coefficients still validate well, exactly as they
do for SD7003/SD8000**: `Cl(alpha)` tracks the UIUC LSAT experimental
polar closely for both airfoils, in both implementations (see
`ClarkY/1-orig/polar_comparison.png`, `GM15/1-orig/polar_comparison.png`);
`Cd` is systematically overpredicted at higher alpha, same as SD7003/SD8000,
for the same physically-explainable reason (no turbulence/transition
model, so no laminar-separation-bubble reattachment physics at this Re).
Python and C++ agree closely with each other at every alpha and every
`dx` level (see each `1-orig/summary.txt`). One instability was found and
fixed the same way SD7003's fine grid was: **ClarkY's coarse grid
(dx=0.04) diverged to NaN at dt=0.01 in both implementations, fixed by
halving dt to 0.005** -- see `ClarkY/README.md` for detail.

## Why not a more famous airfoil at a genuinely lower Re

I looked for a widely-used airfoil (like NACA 0012/4412, S1223, or E387)
with a clean, LSAT-format experimental Re bin meaningfully below SD7003's
61,100, across all 5 UIUC LSAT volumes (which nominally span
Re=30,000-500,000) plus SoarTech 8. In practice, essentially every famous
airfoil's *lowest* clean/untripped LSAT bin sits at Re≈59,000-61,000 --
apparently close to a practical floor for that wind tunnel's clean-flow
measurement quality. The only entries found with clean bins meaningfully
below that (Re≈38,600-40,700) were individual small-model-builder
submissions -- GM15, A18, MA409 -- not widely-cited airfoils. GM15 was
picked as the most complete of these (7 Reynolds-number bins, 8-11 alphas
per bin, vs. 1-2 alphas for some of the others).

## Data provenance

- **Coordinates**: `ClarkY/clarky.dat.txt` converted from
  `https://m-selig.ae.illinois.edu/ads/coord/clarky.dat` (Lednicer format
  -- separate upper/lower surface blocks from the leading edge -- reordered
  into the Selig closed-loop format `py/ibpm.py`'s pipeline expects, same
  as `sd7003.dat.txt`). `GM15/gm15.dat.txt` is
  `https://m-selig.ae.illinois.edu/ads/coord/gm15sm.dat` (already
  Selig-loop format, used as-is).
- **Performance data**: `ClarkY/ClarkY.DRG.txt`/`.LFT.txt` and
  `GM15/GM15.DRG.txt`/`.LFT.txt` are extracted directly from UIUC LSAT
  Volume 1's `DRAG01.TXT`/`LIFT01.TXT` (GM15's "clean" builder-G.-Morris
  block) and Volume 3's `DRAG03.TXT`/`LIFT03.TXT` (Clark-Y's "clean"
  builder-J.-Robertson block) --
  `https://m-selig.ae.illinois.edu/pd/pub/lsat/volume0{1,3}/`. These
  combined-volume files are themselves concatenations of the original
  per-airfoil `.DRG`/`.LFT` files (visible from the `::::::::` filename
  banners and `File ####.DRG created ...` footers between blocks) --
  confirming, incidentally, that this is also almost certainly where
  SD7003/SD8000's own `.DRG`/`.LFT` files came from (`high_re/SD7003/README.md`
  had flagged this as unconfirmed).

## Setup

Same methodology as `high_re/SD7003/README.md`: body = raw UIUC
coordinates resampled to uniform arc-length spacing matched to each run's
grid `dx`, quarter-chord (`0.25, 0`) rotation center, angle of attack via
`py.ibpm`'s `-alpha` flag, domain `length=6, xoffset=-2, yoffset=-1.5`,
`ngrid=1`, `dt=0.01` (`dt=0.005` where noted above/in `ClarkY/README.md`),
Cl/Cd time-averaged over the last 60% of a 3000-step (t=30) run. See
[`../run_low_re_airfoils.py`](../run_low_re_airfoils.py) /
[`../run_low_re_airfoils_cpp.py`](../run_low_re_airfoils_cpp.py) (polar +
grid convergence), [`../run_low_re_flowfield.py`](../run_low_re_flowfield.py) /
[`../run_low_re_flowfield_cpp.py`](../run_low_re_flowfield_cpp.py)
(flowfield case), and [`../gen_low_re_report.py`](../gen_low_re_report.py) /
[`../gen_low_re_flowfield_figs.py`](../gen_low_re_flowfield_figs.py) (figures --
same result-presentation code/style as `high_re/`'s
`gen_airfoil_report_v2.py`/`gen_flowfield_figs_v2.py`, retargeted to
`1-orig/` for these two airfoils).
