# Cylinder flow validation and resolution/Re experiments

- **[`1-baseline/`](1-baseline/)** — the original validation: cylinder at
  Re=100 vs. the published `VORTALL.mat` reference dataset. Clean,
  coherent von Kármán vortex street, matching the reference almost
  exactly (`1-baseline/vorticity_comparison_3way.png`). This is this
  solver's own canonical "clean" result, and the reference point for the
  resolution/Re experiments below — see `../airfoils/README.md`'s
  "Mentor question" section for why it matters (SD7003/SD8000/ClarkY/GM15
  all show broadband vorticity speckle at Re≈40-61k; this cylinder does
  not, at Re=100).
- **[`2-Re_sweep/`](2-Re_sweep/)** — pushes Re UP from the clean Re=100
  baseline toward the airfoils' usual Re~40-61k. Result: clean through
  Re~1000, transitional at Re~3000, fully speckled by Re~10000 — a
  resolution/Re threshold, not a property specific to airfoils.
- **[`3-grid_refine/`](3-grid_refine/)** — at FIXED Re=5000, refines dx
  (0.04/0.02/0.01). Result: coarse dx aliases genuine shear-layer
  instability into broadband speckle; dx=0.02 or finer resolves the same
  physics as organized vortex structures, not noise.

## Bottom line

Two independent experiments here — a Re sweep (`2-Re_sweep/`) and a
grid-refinement study at fixed Re (`3-grid_refine/`) — both converge on
the same answer as the equivalent experiments in
`../airfoils/SD7003/4-Re_sweep/` and `../airfoils/SD7003/5-grid_refine/`
(run on a completely different geometry, sweeping Re from the opposite
direction): the broadband speckle documented throughout this test suite
is a resolution-vs-Reynolds-number threshold effect (boundary-layer
thickness delta ~ c/sqrt(Re) relative to grid spacing dx), not a defect
tied to any specific airfoil, and not something a 1.5x change in Re
(SD7003's 61k vs. GM15's 40k) could ever have revealed. See
`../airfoils/README.md` for the full mentor-question narrative this
completes.

See each subfolder's own `README.md` for full detail; this file is the
map, not the destination.
