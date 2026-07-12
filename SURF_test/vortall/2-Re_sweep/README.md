# Cylinder Reynolds-number sweep (E1 + E2/E3)

Pushes Re UP from `../1-baseline/`'s clean Re=100 (matches `VORTALL.mat`)
toward the airfoils' usual Re~40-61k, at fixed geometry/grid
(`examples/cylinder.geom`, nx=450, ny=200, dx=0.02, same domain as
`1-baseline/`). Originally C++-only; Python added afterward at every Re
(see [`../run_cylinder_re_sweep.py`](../run_cylinder_re_sweep.py) /
[`../run_cylinder_re_sweep_py.py`](../run_cylinder_re_sweep_py.py)).
Figures from
[`../gen_cylinder_re_sweep_figs.py`](../gen_cylinder_re_sweep_figs.py).

**Re>=1000 needed a smaller dt**: at dt=0.02 (the Re=100/500 baseline
dt), Re=1000/3000/10000 all diverged to NaN within the first ~15-30
steps -- not a late-time chaotic blowup, a genuine CFL violation from the
very start (the impulsive-start transient's sharp initial vorticity
gradient becomes too violent for this dt once viscous damping drops
enough). Fixed with dt=0.005 (4x smaller, nsteps=6000 for the same t=30)
for those three Re values, in both implementations.

## Fidelity result: exact agreement through Re=3000, chaos-divergence at Re=10000

`fidelity_summary.txt`: py/ibpm.py and C++ build/ibpm match to 0.00%
relative pointwise difference at Re=100, 500, 1000, and 3000. At Re=10000
(fully speckled/turbulent-looking), the pointwise difference is 93% even
though domain-RMS is close (4.48 vs 4.39) -- the same chaos-amplification
pattern documented in `../1-baseline/README.md` and
`../../airfoils/LSAT-SD7003/4-Re_sweep/README.md`, not a fidelity defect: once
the flow is genuinely chaotic, both implementations remain statistically
correct but diverge in instantaneous phase.

## Result: converges on the same transition zone as the airfoil sweep, from the opposite direction

`re_sweep_comparison.png`:

- **Re=100, 500: clean**, coherent von Kármán-style wake (matches
  `1-baseline/`'s validated Re=100 result).
- **Re=1000: still clean**, though a growing recirculation
  bubble/instability is visible forming.
- **Re=3000: clear transition to broadband speckle.**
- **Re=10000: fully speckled/turbulent-looking**, similar broadband
  character to the airfoils' usual Re~40-61k runs.

`domain_rms_vs_Re.png` confirms this quantitatively (RMS vorticity: 1.53
-> 2.35 -> 2.60 -> 3.16 -> 4.39, Re=100 to 10000 -- growing, accelerating
past Re~1000-3000).

**This transition zone (Re~1000-3000 here) lines up closely with
`../../airfoils/LSAT-SD7003/4-Re_sweep/`'s own transition zone (Re~1000-5000),
found by sweeping a completely different geometry (thin cambered airfoil
vs. blunt cylinder) from the opposite direction** (down from SD7003's
usual ~61k, vs. up from this cylinder's usual Re=100). Two independent
geometries, two independent sweep directions, same answer: this solver's
`ngrid=1`/no-subgrid-dissipation configuration produces clean results up
to roughly Re~1000-3000 at dx=0.02, and broadband speckle well above
that -- consistent with a boundary-layer-thickness-vs-grid-spacing
(delta/dx) threshold, not a property of any particular airfoil or of Re
in isolation.
