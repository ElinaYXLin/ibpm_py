# SD8000 grid refinement at fixed Re=5000 (E4)

Same methodology as `../../SD7003/5-grid_refine/README.md` -- see that
file for the full explanation.

## Result: same pattern as SD7003

`grid_refine_comparison.png`: coarse (dx=0.04) aliases the transitional
shear-layer instability into broadband speckle; medium (dx=0.02) and
fine (dx=0.01) both resolve it as an organized, growing wave that rolls
up into discrete vortices further downstream, not noise. Domain-RMS
vorticity stays roughly flat across dx (2.84 -> 2.81 -> 2.74) while
max|omega| grows (84 -> 166 -> 190) -- same signature as SD7003: finer
grids resolve sharper peaks in genuine flow structure, not spurious
energy. Confirms the `4-Re_sweep/` finding is not SD7003-specific.
