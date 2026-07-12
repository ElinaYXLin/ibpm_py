# Genuinely low-Reynolds-number (hundreds) airfoils: NACA0012 and SD7003

**See [`../SUMMARY.md`](../SUMMARY.md) for the full consolidated
findings** across this directory, `../airfoils/`, and `../vortall/`.

## Why this directory exists (again)

`../airfoils/LSAT-SD7003/4-Re_sweep/` experiments already established that this
solver is clean at Re~200-1000 and only speckles above Re~5000-10000
(see `../airfoils/README.md`'s "Mentor question" section). This
directory follows up on that with two dedicated, from-scratch runs at a
single representative Re actually **in the hundreds** (Re=500) -- the
literal range the mentor originally asked about -- rather than reusing
that sweep's numbers alone, and pairs an "easy" airfoil (new to this
repo) with a "hard" one (SD7003 itself) at the same Re.

## Why no experimental wind-tunnel validation data exists here

Searched for real Cl/Cd polar data at Re in the hundreds (both UIUC
LSAT-style wind-tunnel sources and the broader low-Re aerodynamics
literature). **None exists in a usable tabulated form.** UIUC's own LSAT
facility's practical floor is Re~40,000-60,000 (see
`../airfoils/README.md`'s "Why not a more famous airfoil" section) --
below that, standard force-balance lift measurement and momentum-method
drag measurement both become unreliable in a conventional wind tunnel.
Published low-Re literature (Selig/Deters/Williamson 2011,
"Wind Tunnel Testing Airfoils at Low Reynolds Numbers"; Princeton's
low-Re facility) confirms the practical experimental floor sits around
Re=40,000-60,000 even at specialized low-turbulence low-Re tunnels. Below
that, essentially all published Cl/Cd data is CFD/DNS-computed (e.g.
NACA 0012 has been benchmarked computationally at Re=10-500 in the CFD
literature), not measured -- and even those computational references
report bare Cl(alpha) numbers in papers, not the machine-readable
tabulated format `parse_uiuc.py` expects.

**Consequence: no *experimental* reference polar exists for either
airfoil at this Re** -- unlike `../airfoils/`, which validates against
real UIUC LSAT wind-tunnel data throughout. But a *computational*
(CFD/LBM) drag benchmark does exist for NACA0012 specifically (the most
commonly-used low-Re CFD validation airfoil) -- see below. SD7003 has no
published reference at all, computational or experimental, at Re=500 (it
is a far less common, non-canonical low-Re test case), so it stays a
Python-vs-C++ fidelity + qualitative flow-field check only.

## The two airfoils, both at Re=500

- **[`NACA0012/`](NACA0012/)** -- the single most standard airfoil in
  aerodynamics: thin (12% thickness), symmetric (zero camber), no
  laminar-separation-bubble design intent. The "easy" case. Coordinates
  from `https://m-selig.ae.illinois.edu/ads/coord/n0012.dat` (Lednicer
  format, converted to the Selig closed-loop format this pipeline
  expects, same conversion as `../airfoils/LSAT-ClarkY/`). Two things
  live here: a qualitative flow-field run at alpha=5 deg (see
  [`run_naca0012.py`](run_naca0012.py), `flow_evolution.png`), and a
  **quantitative** Cl/Cd polar + grid-convergence study at alpha=0-10 deg,
  Re=500/1000, validated against published CFD drag benchmarks (Lockard,
  Wu, Nita, Di Ilio, Kurtulus) -- see
  [`NACA0012/README.md`](NACA0012/README.md) for the full writeup
  (`polar_comparison.png`, `grid_convergence.png`,
  `fidelity_summary.txt`; driven by
  [`../run_naca0012_polar.py`](../run_naca0012_polar.py) /
  [`NACA0012/gen_naca0012_report.py`](NACA0012/gen_naca0012_report.py)).
  This is the same content that used to live at
  `airfoils/Lockard-NACA0012/`, merged here so all low-Re NACA0012
  results are in one place.
- **[`SD7003/`](SD7003/)** -- the same cambered, laminar-separation-
  bubble-prone airfoil that originally prompted the mentor's question
  (`../airfoils/LSAT-SD7003/`), now shown at Re=500 instead of its usual
  ~61,100. The "hard" case -- and choosing the *same* airfoil that looked
  "weird" at high Re, rather than a new complex geometry, closes the
  loop directly: this is proof that SD7003 itself, unmodified, is clean
  at low Re. Reuses the Re=500 run already produced by
  `../airfoils/run_airfoil_re_sweep.py`/`run_airfoil_re_sweep_py.py`
  (part of `../airfoils/LSAT-SD7003/4-Re_sweep/`) rather than duplicating the
  simulation -- only the flow_evolution figure is generated fresh here
  (see [`gen_low_re_figs.py`](gen_low_re_figs.py)). No reference polar
  exists for SD7003 at this Re (see above), so this one stays qualitative.

## Result: both clean, both implementations agree

`NACA0012/flow_evolution.png` and `SD7003/flow_evolution.png`: both
airfoils, at Re=500, show a smooth, coherent, laminar wake sheet with no
broadband speckle at every snapshot from t=0 to t=30, in **both**
py/ibpm.py and C++ build/ibpm, visually indistinguishable between the two
implementations at every timestep. This directly confirms, on a brand-new
airfoil (NACA0012) as well as the original one (SD7003), what
`../airfoils/LSAT-SD7003/4-Re_sweep/` already found: **the "weird" broadband vorticity
speckle is a Re~1000-5000+ phenomenon, not present at Re=500, regardless
of airfoil shape or camber complexity.**
