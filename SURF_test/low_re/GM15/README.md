# GM15 airfoil vs. UIUC LSAT experiment

See [`../README.md`](../README.md) for why this directory exists (mentor
question about whether a lower-Re/different airfoil avoids the broadband
vorticity speckle documented for SD7003/SD8000) and full data provenance,
including why a more widely-known airfoil couldn't be found at this Re
level.

GM15: a small free-flight-model airfoil (Gilbert Morris, F1C class),
tested by UIUC LSAT at Re≈40,600 -- **genuinely, substantially lower than
SD7003's 61,100** (not just a different airfoil at the same Re, as
`../ClarkY/` is). Picked as the "hard" case of this pair: more strongly
nonlinear/cambered Cl-alpha behavior than Clark-Y at this Re (see
`1-orig/polar_comparison.png`).

## `1-orig/`

Same content/style as `../ClarkY/1-orig/` and `high_re/`'s
`2-c++included/`: `polar_comparison.png`, `drag_polar.png`,
`grid_convergence.png`, `flow_evolution.png`, `summary.txt` -- Python
`py/ibpm.py` vs. C++ `build/ibpm` vs. UIUC LSAT experiment (see
[`../../gen_low_re_report.py`](../../gen_low_re_report.py) /
[`../../gen_low_re_flowfield_figs.py`](../../gen_low_re_flowfield_figs.py)).

## Results

**Lift matches the experimental polar well** across the tested range
(-3.62° to 9.94°), for both implementations (see
`1-orig/polar_comparison.png`, left panel) -- despite this being the
lowest Re anywhere in this whole test suite. **Drag is overpredicted at
higher alpha**, same pattern and same cause (no transition model) as every
other airfoil here; e.g. at alpha=9.94°: py Cd=0.206, cpp Cd=0.217, vs.
exp Cd=0.026. Python and C++ agree closely with each other at every alpha
and grid level -- no instabilities were encountered anywhere in this
airfoil's sweep (unlike `../ClarkY/`'s coarse-grid case), despite the
lower Re.

**The vorticity field (`1-orig/flow_evolution.png`) still shows the same
broadband grid-scale speckle as SD7003/SD8000 and ClarkY**, in both
implementations, at this substantially lower Re -- the headline (negative)
result this airfoil was specifically run to check. See `../README.md`'s
"Results" section for why: the speckle traces to the shared `ngrid=1`
single-grid / no-subgrid-dissipation domain configuration used everywhere
in this suite, not to any airfoil's Reynolds number.
