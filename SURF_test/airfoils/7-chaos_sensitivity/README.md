# Chaos sensitivity: Re=200 control, extended coarse-grid runs, and a perturbation ensemble

Follow-up to `2-c++included/port_fidelity_diagnostic.py`'s finding that
SD8000's coarse grid (dx=0.04) shows py/cpp disagreeing sharply (even
opposite-signed time-averaged Cl) at the tail of a 3000-step run, traced
there to a single late-run numerical blow-up in C++ dominating a naive
time-average. This directory pushes that investigation further in four
directions, all proposed (not run) in an earlier pass and executed here:

1. **A clean control case (Re=200, no experimental reference needed)** --
   confirms py and cpp agree *exactly* (not just statistically) when the
   flow is genuinely steady/laminar, so the coarse-grid disagreement is
   not a general property of this solver, only of the unstable regime.
2. **Extending SD7003 and SD8000's coarse-grid (dx=0.04) runs from 3000 to
   4000 steps** -- tests whether Python's trajectory (bounded within the
   original 3000-step window) also eventually blows up given more time.
3. **Phase-space (Cl vs Cd) plots**, contrasting the clean Re=200 case
   against the chaotic coarse-grid case.
4. **A tiny-Re-perturbation ensemble** (16 runs, C++ only, Re varied by
   only 5e-8 to 1e-5 relative) with a histogram of blow-up step, directly
   testing whether the instability's timing is chaotically sensitive
   (butterfly-effect) rather than a fixed, deterministic artifact.

## Files

| File | What it is |
|---|---|
| `run_chaos_sensitivity.py` | Launches the 4 extended-to-4000-step runs (SD7003/SD8000 x py/cpp) and the 16-run perturbation ensemble (C++, SD8000 config), all in parallel via `subprocess.Popen` (one OS process per run -- this machine has 10 cores; all 20 runs completed in ~72s wall-clock). |
| `gen_chaos_sensitivity_figs.py` | Reads `_run_data/` (this directory's new runs) plus the existing `4-Re_sweep/Re200` data (no rerun needed for that part) and writes every figure below. |
| `re200_comparison.png` | Cl(t), Cd(t), and phase-space Cl-vs-Cd for both airfoils at Re=200, py vs cpp overlaid. |
| `phase_space_coarse_ext4000.png` | Phase-space Cl vs Cd for the extended coarse-grid runs, both airfoils. |
| `blowup_histogram.png` | Histogram of blow-up timestep across the 16-run Re-perturbation ensemble. |
| `SD8000_ext4000_py_snapshots.png`, `SD8000_ext4000_cpp_snapshots.png`, `SD7003_ext4000_py_snapshots.png` | Vorticity field snapshots at the last several stable steps (every 25 steps) leading up to each run's blow-up, plus the first post-blow-up snapshot. |
| `_run_data/` | Raw output: `<airfoil>_ext4000_<impl>/` (the 4000-step extensions, restart every 25 steps) and `SD8000_perturb_<NN>/` (the 16-run ensemble, force traces only). `.bin`/`.cholesky` files are gitignored, as elsewhere in this repo; `.force`/`.cmd`/`run_log.txt` are committed. |

Regenerate with:
```bash
python3 SURF_test/airfoils/7-chaos_sensitivity/run_chaos_sensitivity.py   # ~70-90s, needs build/ibpm
python3 SURF_test/airfoils/7-chaos_sensitivity/gen_chaos_sensitivity_figs.py
```

## Results

### 1. Re=200: exact agreement, clean phase space

`max|Cl_py - Cl_cpp| = 0.000e+00`, `max|Cd_py - Cd_cpp| = 0.000e+00`,
across all 3001 recorded steps, for **both** SD7003 and SD8000 (existing
restart data from `4-Re_sweep/Re200`, no rerun needed). Phase space
(`re200_comparison.png`, right column) shows a clean inward spiral
settling onto a single fixed point -- fully steady/laminar, nothing to
amplify. This is the deterministic baseline the unstable coarse-grid case
below should be read against.

### 2. Extending to 4000 steps: the outcome flipped, with nothing deliberately changed

The original coarse runs used `-restart 0` (no checkpoints saved), so
extending required a fresh 0->4000 rerun of each -- same Re, same alpha,
same everything explicit as the original 3000-step runs.

| case | blow-up step |
|---|---|
| SD7003 py | 3959 |
| SD7003 cpp | no blow-up through 4000 |
| SD8000 py | 2895 |
| SD8000 cpp | 2970 |

Previously (`2-c++included/port_fidelity_diagnostic.py`), only **C++**
blew up (at step 3000) while Python looked bounded through that window.
Here, in a fresh rerun with nothing explicitly changed, **Python** blows
up *earlier* than C++ for SD8000, and for SD7003 it's Python that blows
up while C++ stays bounded through 4000 steps -- the opposite pattern.
This is itself evidence for the chaos-sensitivity explanation: this
codebase's `FFTW_EXHAUSTIVE` planning is documented elsewhere
(`SURF_test/cost/1-multi-core/README.md`) as not perfectly reproducible
run-to-run (it re-times and can select a different, numerically-equivalent
plan each time), which is enough of a last-bit perturbation to reshuffle
which implementation crosses this resolution's instability threshold
first.

### 3. Phase space: chaotic attractor, not a porting difference

`phase_space_coarse_ext4000.png`: both py and cpp trace tangled,
non-repeating trajectories occupying a similar overall region and shape
-- not identical paths, but visually the same *kind* of attractor. This
contrasts directly with Re=200's clean spiral-to-a-point, and is the
qualitative signature of two independent samples of the same chaotic
dynamics, not two different algorithms.

### 4. Vorticity snapshots: the domain is already saturated well before the force-trace "blow-up"

For all three blow-up cases, the snapshots leading up to blow-up (every
25 steps) show the vorticity field **already fully broadband-saturated
across the entire domain** several time-units beforehand -- not a
localized event. The field then goes **fully to NaN** within 1-2 recorded
snapshots (confirmed directly, `np.isnan(w).all() == True`) after the
force-trace threshold crossing. The visible "blow-up" in Cl/Cd is the
tail end of a much longer-building, whole-domain degradation.

### 5. Perturbation ensemble: 14/16 blow up, spread over 1805 steps

16 C++-only runs, SD8000's coarse config, Re perturbed by only
**5e-8 to 1e-5 relative** (the 6th-8th significant digit of Re=60800):
**14 of 16 blow up**, at steps ranging from **2105 to 3910**; **2 of 16**
never blow up through 4000 steps. See `blowup_histogram.png`. No
discernible relationship between perturbation sign/magnitude and blow-up
timing -- an imperceptible change to one parameter scrambles both
*whether* and *when* the instability crosses its threshold, which is
about as direct a demonstration of chaotic sensitivity as this kind of
test can produce.

## Conclusion

Combined with `2-c++included/port_fidelity_diagnostic.py`'s original
finding, this closes the loop: the coarse-grid (dx=0.04) disagreement
between py and cpp is fully explained by both implementations sharing the
same genuinely chaotic, marginally-unstable dynamics at this
under-resolved configuration -- not a porting defect. Re=200 (steady,
non-chaotic) shows exact agreement; the unstable regime shows the
textbook signature of chaos (chaotic-attractor phase space, blow-up
timing scrambled by imperceptible parameter changes, and even a fresh
rerun with nothing deliberately changed reshuffling which implementation
happens to cross the threshold first).
