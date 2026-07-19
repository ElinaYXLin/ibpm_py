# paper_figs/

Six small crops crops of Kurtulus (2019)'s own wake-vorticity figures, at
the same 3 angles (α₀=0°, 9°, 12°) `gen_kurt_figs.py`'s `wake_contours()`
plots for `py_static`/`cpp_static`, so `wake_steady.png` and
`wake_f4hz.png` can show what's actually being compared against instead of
only ibpm's own output.

- `steady_a{00,09,12}.png` — rows from the paper's Figure 2 (instantaneous
  vorticity, non-oscillating NACA0012, t=100s).
- `f4hz_a{00,09,12}.png` — rows from the paper's Figure 6 (instantaneous
  vorticity, f=4Hz pitching, t=36s).

Regenerate with `python3 ../extract_paper_figs.py` (crops the source PDF
by hardcoded pixel offsets — see that script for exact provenance and
caveats). Used here strictly for side-by-side scientific comparison against
this repo's own results, with the source explicitly labeled in every figure
that uses them.
