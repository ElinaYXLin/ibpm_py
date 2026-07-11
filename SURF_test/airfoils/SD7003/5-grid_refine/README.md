# SD7003 grid refinement at fixed Re=5000 (E4)

At a FIXED Reynolds number (Re=5000, a point in `../4-Re_sweep/`'s
transitional zone -- reused directly from that sweep's dx=0.02 run
rather than rerun), refines dx across coarse/medium/fine (0.04/0.02/0.01,
this suite's usual three convergence levels) to separate two variables
that `4-Re_sweep/` alone conflates: is it Re, or is it resolution (at a
given Re), that controls whether the field looks clean or speckled?
C++ only, same rationale as the Re sweep. See
[`../../run_airfoil_grid_refine.py`](../../run_airfoil_grid_refine.py) /
[`../../gen_airfoil_grid_refine_figs.py`](../../gen_airfoil_grid_refine_figs.py).

## Result: resolution matters, but not the way "just refine dx" suggests

`grid_refine_comparison.png`, at fixed Re=5000:

- **dx=0.04 (coarse): broadband, grid-scale speckle** spread through the
  whole wake and beyond -- the classic "weird" look.
- **dx=0.02 (medium): a clean, smooth, coherent wake sheet** with a
  gently growing instability -- no speckle.
- **dx=0.01 (fine): the SAME instability, but now resolved as an
  organized wave train rolling up into discrete vortices** (visible
  Kelvin-Helmholtz-type shear-layer roll-up) -- also not speckled, just
  more detail than the medium grid could show.

So refining dx at Re=5000 does NOT show a simple monotonic
"coarser=messier, finer=cleaner" story -- the coarsest level aliases a
genuine, organized shear-layer instability into broadband noise (grid
too coarse to represent the wave, so it folds into noise); adequate
resolution (dx=0.02 or finer) reveals that same physics as ordered
vortex dynamics, not noise. Domain-RMS vorticity stays roughly flat
across all three dx (2.84 -> 2.75 -> 2.68), while max|omega| grows
(68 -> 171 -> 201) -- consistent with `../6-explicit_dissipation/`'s and
`../3-dx0.01/`'s finding elsewhere in this suite: finer grids resolve
sharper peaks in genuine flow structure, they don't add spurious energy.

**Combined with `../4-Re_sweep/`'s finding**: the controlling parameter
is resolution *relative to* the Reynolds number (delta/dx, where delta ~
c/sqrt(Re)) -- at Re=61,100 (SD7003's usual Re), even dx=0.01 is nowhere
near fine enough to resolve the ~0.004c boundary layer, which is why
`../3-dx0.01/` (dx=0.01 at Re=61,100) still speckles. At Re=5000, dx=0.02
is already enough to resolve the (much thicker, ~0.014c) boundary layer
cleanly.
