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

See each subfolder's own `README.md` for full detail; this file is the
map, not the destination.
