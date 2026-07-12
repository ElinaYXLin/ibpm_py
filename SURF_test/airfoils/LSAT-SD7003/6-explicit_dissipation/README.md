# SD7003: does explicit dissipation clean up the speckle? (E5)

**Not a solver modification.** `py/ibpm.py` and `build/ibpm` are
DNS-style solvers with no subgrid/turbulence model, by design -- faithful
ports of the unmodified `cwrowley/ibpm` method (confirmed byte-identical
to upstream, see `../README.md`'s "Confirmed" section). Actually adding a
dissipation term would mean forking the timestepper itself (touching the
RK3 substep / elliptic-solve chain), a much larger and riskier change
than any other experiment in this suite, and the result would no longer
be "the same solver" this whole test suite otherwise validates the
fidelity of.

Instead, `filter_demo.png` (from
[`../gen_dissipation_demo.py`](../gen_dissipation_demo.py)) is a cheap,
honest proxy: it takes the ALREADY-COMPUTED, speckled Re=40000 vorticity
snapshot (`../4-Re_sweep/_run_data_cpp/Re40000/run03000.bin`) and applies
a Gaussian spatial filter post-hoc, purely as a visualization -- "if this
field had subgrid dissipation damping grid-scale content, roughly this
is what it would look like." **This does not rerun the physics and
proves nothing about solver correctness.**

## Result

Filtering at increasing sigma (0, 1, 2, 4 grid cells) progressively
reveals large-scale coherent vortex structures underneath the speckle,
and monotonically reduces both metrics:

| sigma (cells) | RMS \|omega\| | max \|omega\| |
|---|---|---|
| 0 (raw) | 8.10 | 218.6 |
| 1 | 3.07 | 52.8 |
| 2 | 2.17 | 30.2 |
| 4 | 1.45 | 14.5 |

This is consistent with (not new evidence beyond) `../4-Re_sweep/`'s and
`../5-grid_refine/`'s finding: the speckle is high-spatial-frequency
content sitting on top of genuine larger-scale vortex dynamics. A real
subgrid model would very likely produce a visually similar cleanup (that
is what such models are for) -- but confirming that properly would
require actually implementing and validating one in the solver, which is
out of scope here. The resolution/Re story (`../4-Re_sweep/`,
`../5-grid_refine/`, `../../../vortall/2-Re_sweep/`,
`../../../vortall/3-grid_refine/`) already fully explains the phenomenon
without needing to invoke a missing dissipation term.
