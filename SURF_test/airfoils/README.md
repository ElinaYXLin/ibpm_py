# Airfoil validation: SD7003, SD8000, ClarkY, GM15

This directory holds all airfoil validation cases for `py/ibpm.py` vs.
`build/ibpm` vs. UIUC LSAT (Low-Speed Airfoil Tests) wind-tunnel data:

| Airfoil | Re | Why picked |
|---|---|---|
| [`SD7003/`](SD7003/) | ≈61,100 | Cambered, laminar-separation-bubble-prone low-Re airfoil; the original case, extensively wind-tunnel tested at UIUC |
| [`SD8000/`](SD8000/) | ≈60,800 | Sibling Selig-Donovan design, run for direct comparison against SD7003 |
| [`ClarkY/`](ClarkY/) | ≈60,700 | Flat-bottomed, historically the most widely-used airfoil in aviation; same Re ballpark as SD7003 -- isolates "does a different, more mainstream airfoil geometry help at the same Re" |
| [`GM15/`](GM15/) | ≈40,600 | Small free-flight-model airfoil; genuinely, substantially lower Re than SD7003 -- isolates "does meaningfully lower Re help" |

(ClarkY/GM15 were originally in a separate `low_re/` directory, distinct
from SD7003/SD8000's `high_re/`; both were merged here once the
lower-Re-alone hypothesis was ruled out -- see "Mentor question" below.)

## Mentor question: why does the vorticity field look "weird"?

A mentor flagged that the vorticity field in SD7003's flow evolution
(`SD7003/2-c++included/flow_evolution.png`) shows broadband grid-scale
speckle, instead of the clean, coherent vortex structures textbooks show.
The investigation proceeded in stages, each ruling out a candidate cause,
each documented where it happened:

1. **Wrong initial conditions?** No -- `SD7003/README.md`/`SD7003/3-small_dt/`
   confirmed the impulsive-start IC and Re=61,100 already match the UIUC
   LSAT dataset's own stated test condition.
2. **Far-field domain too small (`ngrid=1`)?** No -- `SD7003/3-ngrid=3/`
   tested the solver's own documented multi-domain far-field scheme;
   speckle unchanged (and a new instability appeared past t≈20).
3. **Under-resolved (`dx=0.02`)?** No -- `SD7003/3-dx0.01/` showed finer
   dx does not clean it up either; peak vorticity grows sharper without
   converging (see `SD7003/3-ngrid=3/instability_diagnostic.png` panel C).
4. **The specific airfoil, or its Re specifically (≈61k)?** No -- ClarkY
   (same Re, different airfoil) and GM15 (different airfoil, genuinely
   lower Re=40,600) both still speckle identically -- see each folder's
   `1-orig/flow_evolution.png`.

That last (negative) result led to a fifth stage: comparing against
`SURF_test/vortall/`, this solver's own canonical clean benchmark (cylinder
at Re=100, matches the published `VORTALL.mat` reference dataset with a
textbook von Kármán street, **no speckle at all**, same solver, same
`ngrid=1` configuration). If the solver/domain isn't the cause, and
neither is the specific airfoil or a 1.5x change in Re, what differs
between the clean Re=100 cylinder and the speckled Re≈40-61k airfoils?

**Answer: resolution relative to the viscous length scale, not Re as an
independent variable.** The near-wall vorticity layer thickness scales as
δ ~ c/√Re. At Re=100 (cylinder), δ≈0.1c is ~5 grid cells wide at dx=0.02 --
comfortably resolved. At Re≈40-61k (airfoils), δ≈0.004-0.005c is well
under *one* grid cell at the same dx -- the boundary layer itself is
aliased by the grid, which is what generates the broadband noise (no
subgrid/turbulence model exists to compensate, by design, since this is a
DNS-style solver). Re=40k vs. 61k (a 1.5x change) doesn't cross this
threshold; the experiments below (numbered `4-`/`5-`/`6-` in each airfoil
folder, `2-`/`3-` in `../vortall/`) test this directly by sweeping Re and
dx far enough to see the clean-to-speckled transition happen.

## Numbering convention

Each airfoil folder numbers its subdirectories in the order experiments
were run; folder names describe what changed relative to the baseline.
`SD7003/` and `SD8000/` have pre-existing `1-`/`2-`/`3-` folders from
earlier stages (see each folder's own README); new resolution/Re
experiments continue from `4-`. `ClarkY/`/`GM15/` only have `1-orig/`, so
new experiments there would start at `2-` (not run in this pass -- see
each folder's own README for why).

See each subfolder's own `README.md` for full detail; this file is the
map, not the destination.

## Why not a more famous airfoil at a genuinely lower Re (ClarkY/GM15 choice)

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

## Data provenance (ClarkY/GM15; see each of SD7003/SD8000's own README for theirs)

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
  SD7003/SD8000's own `.DRG`/`.LFT` files came from (`SD7003/README.md`
  had flagged this as unconfirmed).

## Setup (ClarkY/GM15's `1-orig/`)

Same methodology as `SD7003/README.md`: body = raw UIUC coordinates
resampled to uniform arc-length spacing matched to each run's grid `dx`,
quarter-chord (`0.25, 0`) rotation center, angle of attack via
`py.ibpm`'s `-alpha` flag, domain `length=6, xoffset=-2, yoffset=-1.5`,
`ngrid=1`, `dt=0.01` (`dt=0.005` for ClarkY's coarse grid -- see its
README), Cl/Cd time-averaged over the last 60% of a 3000-step (t=30) run.
Driven by [`../run_clarky_gm15.py`](../run_clarky_gm15.py) /
[`../run_clarky_gm15_cpp.py`](../run_clarky_gm15_cpp.py) (polar + grid
convergence), [`../run_clarky_gm15_flowfield.py`](../run_clarky_gm15_flowfield.py) /
[`../run_clarky_gm15_flowfield_cpp.py`](../run_clarky_gm15_flowfield_cpp.py)
(flowfield case), and [`../gen_clarky_gm15_report.py`](../gen_clarky_gm15_report.py) /
[`../gen_clarky_gm15_flowfield_figs.py`](../gen_clarky_gm15_flowfield_figs.py)
(figures -- same result-presentation code/style as SD7003/SD8000's
`gen_airfoil_report_v2.py`/`gen_flowfield_figs_v2.py`).
