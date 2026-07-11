# Cylinder grid refinement at fixed Re=5000 (E4)

At FIXED Re=5000 (a point in `../2-Re_sweep/`'s transitional zone),
refines dx across coarse/medium/fine (0.04/0.02/0.01) -- same question as
`../../airfoils/SD7003/5-grid_refine/`, on the cylinder instead. Circle
geometry regenerated fresh at each dx via `circle_n` (point spacing
tracks dx at every level, same convention `make_airfoil_raw.py` uses for
the airfoils). See
[`../run_cylinder_grid_refine.py`](../run_cylinder_grid_refine.py) /
[`../gen_cylinder_grid_refine_figs.py`](../gen_cylinder_grid_refine_figs.py).

**Also needed smaller dt**: at Re=5000, both dx=0.02 (dt=0.02) and
dx=0.01 (dt=0.01) diverged within the first ~15-18 steps -- same
early-blowup CFL signature as `../2-Re_sweep/`'s Re>=1000 cases, but this
blunt body's impulsive-start vorticity layer needed an even smaller dt
than the equivalent-Re airfoil case did. Fixed with dt=0.005/nsteps=6000
(medium) and dt=0.0025/nsteps=12000 (fine), same t=30.

## Result: same pattern as the airfoils -- coarse aliases, fine resolves

`grid_refine_comparison.png`:

- **dx=0.04 (coarse): broadband speckle** throughout the wake -- the
  familiar "weird" look, same character as `../2-Re_sweep/`'s Re=10000
  panel despite this being a 2x-lower Re.
- **dx=0.02 (medium): clean, organized recirculation** -- large coherent
  vortex structures, no speckle.
- **dx=0.01 (fine): the same organized structures**, now with additional
  resolved small-scale shear-layer roll-up near the body -- more detail,
  not more noise.

Domain-RMS vorticity drops from coarse to medium/fine (4.03 -> 3.46 ->
3.43) while max|omega| grows (84 -> 123 -> 162) -- same signature the
airfoil grid-refine experiments show: the coarse grid aliases genuine
shear-layer/wake instability into broadband noise; adequate resolution
resolves it as organized vortex dynamics with sharper (not spurious)
peaks. Confirms `../../airfoils/SD7003/5-grid_refine/`'s finding is not
airfoil-specific -- it reproduces on a completely different body shape.
