# SD8000 Reynolds-number sweep (E1 + E2)

Same methodology as `../../SD7003/4-Re_sweep/README.md` (see that file
for the full explanation) -- Re sweep DOWN from SD8000's usual ~60,800 at
fixed dx=0.02, alpha=5.36.

## Result: same transition, same conclusion, on the sibling airfoil

`re_sweep_comparison.png` / `domain_rms_vs_Re.png` show the identical
pattern SD7003 showed: clean coherent wake at Re=200-1000, transitional
waviness at Re=5000-10000, full broadband speckle by Re=20000-40000
(matching `../2-c++included/flow_evolution.png`'s Re=60,800 result).
Domain-RMS vorticity: 1.27 (Re=200) -> 1.56 -> 1.85 -> 2.81 -> 3.70 ->
4.81 -> 6.08 (Re=40000) -- same accelerating-past-Re~1000-5000 shape as
SD7003's own curve. Running the sweep on both SD7003 and SD8000 (rather
than just one) confirms this is a property of the shared solver/domain
configuration acting on any airfoil at these Re, not something specific
to SD7003's particular geometry.
