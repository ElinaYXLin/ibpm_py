# Multi-Domain Grid Upgrade Manual

Instructions manual for implementing grid-shape and coarsening-ratio upgrades to the
IBPM multi-domain solver in `py/`. Written for a Claude (Sonnet/Opus) implementation
session. **This document contains no code — it specifies what to build, where, the
math, the invariants to preserve, and the acceptance tests.**

---

## 0. Executive summary of feasibility

| Idea | Verdict | Why |
|---|---|---|
| Rectangular (non-square) grid levels | **Already supported** | `nx ≠ ny` is legal today; cells stay square (`dx = length/nx`), the domain is a rectangle. All levels share the same `nx × ny` and aspect ratio. |
| Pixelated U-shaped levels | **Not doable** in this architecture | The per-level fast solver is an FFTW DST-I (sine transform), which requires a tensor-product rectangle. Every transfer operator (`coarsify`, `getBC`, the BC helpers in `vector_operations.py`) assumes four straight edges. A U-shaped Poisson solve needs domain decomposition or an iterative solver — a research project, and almost certainly slower than one extra rectangular DST level at these sizes. |
| Auto-detect body shape → grid plan | **Doable, easy, high value** | Pure pre-processing. No solver changes. See Phase 1. |
| Coarsening ratio 3 | **Not recommended** | Odd ratio breaks the even/odd node-coincidence structure that `getBC`, `NxExt`, and the shift quantization are built on. Far more invasive than ratio 4, for less benefit. |
| Coarsening ratio 4 | **Doable, but bounded payoff** | Keeps node coincidence. Full spec in Phase 3. Speedup is capped (see §0.1) and interface accuracy degrades, so this is lower priority than Phases 1–2. |
| Per-level independent rectangle extents | **Doable, large refactor** | The honest version of "shape-adaptive levels." Design sketch in Phase 4; do only if Phases 1–2 prove insufficient. |

### 0.1 Why ratio 4 buys less than it looks like it should

Every level has the **same** `nx × ny` cell count (one array of shape
`(ngrid, nx−1, ny−1)` in `Scalar`), so per-timestep cost of the multi-domain
machinery is proportional to `ngrid`. The outermost domain extent is
`length · R^(ngrid−1)`. To cover an outer extent `D` starting from finest extent `L`:

    ngrid = 1 + ceil( log_R (D / L) )

Example: `L = 2c`, `D = 32c` → R=2 needs 5 levels, R=4 needs 3 levels.
That removes 40% of the *elliptic/transfer* work — but the timestep `dt` is still
CFL-limited by the finest `dx`, the immersed-boundary projection cost is unchanged,
and the finest-level work is unchanged. Realistic end-to-end speedup: **1.3–1.7×**.

Meanwhile the resolution jump at each interface goes from 2× to 4×: coarse-to-fine
boundary interpolation error grows ~(R·dx)², and vortices crossing an interface see
a 4× dissipation jump (wake structures partially reflect/smear). The Kurtulus
validation work in `SURF_test/low_re/NACA0012/kurt_comp/` already showed Strouhal
sensitivity to `ngrid`/domain choices — R=4 will make that worse, so it must ship
behind a preset with R=2 remaining the default.

### 0.2 What actually helps at dx ≈ 0.0015c (do these first)

At `dx = 0.0015c` with a body-fitting finest box (say 1.5c × 0.5c), you need roughly
`nx ≈ 1000, ny ≈ 350`. At that size the dominant *avoidable* costs in this port are:

1. **FFTW planning waste** — every `EllipticSolver2d` (one per level, per Poisson
   solver, per Helmholtz solver, per RK substep) builds its own `NativeDST2D` with
   `FFTW_EXHAUSTIVE` on the *identical* transform size `(nx−1) × (ny−1)`.
   At nx ≈ 1000+ each exhaustive plan can take minutes; there may be a dozen of them.
   One shared plan (or cached wisdom) eliminates all but the first. Zero numerical risk.
2. **Python-loop hotspots in the transfer operators** — `Scalar.getBC`
   ([scalar.py:176](py/scalar.py)) gathers coarse-grid values with per-element list
   comprehensions calling `self(lev+1, i, j)` one point at a time. That is O(nx)
   interpreted-Python calls per boundary edge, per level, per elliptic solve, per
   substep — order 10⁵ Python calls per timestep at nx ≈ 2000. Vectorizing to numpy
   slices is a mechanical change with a large payoff.
3. **Grid placement** — the existing `xshift`/`yshift` parameters
   ([grid.py:201](py/grid.py:201)) already let coarse levels sit off-center (e.g. bias
   the coarse levels downstream for the wake while the fine box hugs the body). Most
   runs don't exploit this; the Phase-1 planner should.

These are Phases 1–2 and they are the recommended path. Phase 3 (ratio 4) is
optional after that.

---

## 1. Architecture primer (read before touching anything)

Facts the implementer must internalize, with source locations:

- **Grid** ([py/grid.py](py/grid.py)): scalar parameters only. `dx = length/nx`;
  cells are square; `Dx(lev) = dx·2^lev` (line 155, `1 << lev`); level `lev` extent
  is `2^lev · length`, positioned by `_getXOffset/_getYOffset` (lines 242–248, uses
  `(1 << lev) − 1`). The fine domain occupies the middle half of the next coarser
  level, modulated by shift: `NxExt = nx/4 · (1 − xShift)` = number of coarse cells
  left of the fine domain (line 111). Index maps `c2f`/`f2c` hard-code ×2 and //2
  (lines 121–140). `resize` asserts `nx % 4 == 0 and ny % 4 == 0` (line 79);
  `setXShift` asserts `xShift·nx % 4 == 0` (line 204).
- **Scalar** ([py/scalar.py](py/scalar.py)): one array `(ngrid, nx−1, ny−1)` holding
  node-centered values at *all* levels (line 73). `coarsify()` (line 85) restricts
  fine → coarse with the 9-point full-weighting stencil (weights 1/4, 1/8, 1/16 —
  the R=2 hat function). `getBC(lev, bc)` (line 176) builds the fine-level boundary
  values from the next coarser level: fine boundary nodes with even index coincide
  with coarse nodes (copied), odd-index nodes are the average of their two even
  neighbors. **This even/odd split is the deepest factor-2 assumption in the code.**
- **Flux** ([py/flux.py](py/flux.py)): edge fluxes per level; no coarsify of its own.
- **EllipticSolver** ([py/elliptic_solver.py](py/elliptic_solver.py)): the
  multi-domain driver — coarsify RHS, solve coarsest with zero BCs, walk down levels
  pulling BCs from the level above via `getBC`. Line 86: `dx = self._dx * (1 << lev)`.
- **EllipticSolver2d** ([py/elliptic_solver_2d.py](py/elliptic_solver_2d.py)): one
  uniform-rectangle solve via FFTW DST-I with eigenvalue division. Requires the
  level to be a full rectangle of uniform square cells — this is what kills U-shapes.
  Each instance owns a `NativeDST2D` planned with `FFTW_EXHAUSTIVE`
  ([py/_fftw_native.py](py/_fftw_native.py)).
- **vector_operations.py**: `Curl`, `Laplacian`, flux↔velocity conversions do the
  inter-level coupling through `getBC`; `CrossProduct` calls `Scalar.coarsify`
  (line 731). Audit these whenever a transfer operator changes.
- **Regularizer / geometry**: the immersed body lives **entirely on the finest
  level** ([py/regularizer.py](py/regularizer.py), "Checks only the finest grid
  level"). The IB delta function assumes isotropic spacing `grid.Dx()` — this is why
  anisotropic cells (dx ≠ dy) are off the table.
- **CLI** ([py/ibpm.py:182–189](py/ibpm.py:182)): the eight knobs a planner must
  emit: `nx, ny, ngrid, length, xoffset, yoffset, xshift, yshift`.

**Grep checklist for factor-2 assumptions** (run before and after Phase 3; every hit
must be either generalized or justified in a comment):
`"1 << lev"`, `"// 2"`, `"* 2"`, `"// 4"`, `"% 4"`, `"0.25"`, `"0.125"`, `"0.0625"`,
`"[::2]"`, `"[1::2]"` across `py/grid.py`, `py/scalar.py`, `py/flux.py`,
`py/elliptic_solver.py`, `py/vector_operations.py`, `py/bc.py`, `py/state.py`,
`py/output_*.py`.

---

## Phase 0 — Baseline harness (do this before any change)

**Goal:** a reproducible before/after check so every later phase can prove it changed
nothing (or changed only what it claims).

1. Pick two reference runs: (a) a small stationary-cylinder or NACA case,
   `ngrid ≥ 3`, ~200 steps; (b) one existing Kurtulus pitching case from
   `SURF_test/low_re/NACA0012/kurt_comp/` at reduced resolution, ~1000 steps.
2. Record, per run: force history (CL, CD time series), final-state checksum
   (e.g. L2 norms of omega per level), and wall-clock broken down by phase.
3. Add coarse timing instrumentation (a context-manager timer or cProfile wrapper —
   whichever exists most naturally) around: FFTW plan construction, per-level DST
   solves, `Scalar.getBC`, `Scalar.coarsify`, the projection (Cholesky/CG) solve,
   and the nonlinear term. Emit a one-line summary table at end of run.
4. Save both baselines under a `benchmarks/` directory with the exact command lines.

**Acceptance:** re-running a baseline twice gives identical force histories
(deterministic), and the timing table clearly identifies the top 3 costs.
Per the project's conventions, plot the timing breakdown (bar chart) and the CL/CD
traces, and save them alongside the data.

---

## Phase 1 — Auto grid planner ("plan mode")

**Goal:** a new standalone module (suggested: `py/grid_planner.py`) plus a CLI
entry point that reads the geometry (and motion, if any), detects the body's swept
envelope, and emits an optimal set of the eight existing grid knobs. **No solver
changes.** This delivers the user-facing "auto-detect the shape and make a plan"
feature — for a short-and-wide airfoil it will naturally produce a short-and-wide
finest rectangle, using capability the solver already has.

### 1.1 Inputs

- Geometry file (same format `Geometry.load` reads).
- Motion description: either "stationary" or the motion parameters already in the
  input deck (pitching amplitude/axis, plunge amplitude, etc.).
- Targets: finest `dx` (e.g. 0.0015), outer domain extent `D` (e.g. 30c),
  margins — upstream, transverse, downstream (defaults: 0.75c, 0.75c, 2.0c;
  downstream larger because the near wake must stay in the finest box).
- Coarsening ratio R (default 2; accepts 4 after Phase 3 exists).

### 1.2 Algorithm

1. Load the geometry; get body points.
2. **Swept envelope:** if the body moves, sample the motion over one full period
   (≥ 36 phases), apply each transform to the body points (reuse the existing
   `RigidBody`/`Motion` machinery — e.g. `geom.moveBodies(t)` at sampled times),
   and take the union bounding box. For stationary bodies this is just the bbox.
3. Pad the envelope with the three margins → finest-domain rectangle
   `[x0, x1] × [y0, y1]`.
4. `nx = ceil((x1−x0)/dx)` rounded **up** to the divisibility requirement
   (multiple of 4 for R=2; multiple of 2R in general). Same for `ny`.
   Prefer values where `nx−1` and `ny−1` have small prime factors (the DST-I of
   size n costs like an FFT of size 2(n+1)); when rounding up, scan the next few
   admissible values and pick the most FFT-friendly.
5. `length = nx·dx`; `xoffset = x0` (re-center the slack from rounding);
   `yoffset` likewise.
6. `ngrid = 1 + ceil(log_R(D / max(length, ny·dx)))`, then report the achieved
   outer extents in both directions.
7. **Shifts:** choose `xshift` so the body sits ~25–30% from the upstream edge of
   the *outermost* domain (wake bias), quantized to the legal step
   (`xshift·nx ≡ 0 mod 4` for R=2). `yshift = 0` unless the motion envelope is
   vertically asymmetric.
8. Sanity checks: warn if Lagrangian point spacing on the body differs from `dx` by
   more than ~20% (reuse the logic in [py/checkgeom.py](py/checkgeom.py), and
   recommend regenerating the geometry at spacing ≈ dx); warn if the finest level
   contains > ~8M nodes; error if the body envelope isn't strictly inside the
   finest rectangle including margins at every sampled phase.

### 1.3 Outputs

- Printed plan: the eight knob values, node counts, per-level dx and extents,
  memory estimate, and the exact ready-to-paste command line / input-deck fragment.
- **A figure** (matplotlib, saved PNG): nested level rectangles drawn to scale with
  the body silhouette and its swept envelope overlaid. This is mandatory — it is
  the whole point of "plan" mode, and it follows the project convention that every
  analysis output ships with a graph.
- Optional `--apply` flag that writes/patches an input deck.

### 1.4 Acceptance

- Planner on the existing Kurtulus NACA0012 pitching case reproduces (or improves
  on) the hand-chosen grid from `kurt_comp` and its figure clearly shows the swept
  envelope inside the finest box.
- Planner output for a stationary cylinder matches a hand computation.
- A run launched with planner-emitted knobs starts and passes the Grid asserts.

**Effort estimate:** small. One module + tests + figure. No solver risk.

---

## Phase 2 — Performance engineering for very small dx

**Goal:** remove the two measured hotspots from §0.2 without changing any numbers.

### 2.1 FFTW plan sharing / wisdom

Current behavior: every `EllipticSolver2d.__init__` builds a fresh `NativeDST2D`
with `FFTW_EXHAUSTIVE` for the same `(nx−1, ny−1)` shape. (The port notes in
[py/elliptic_solver_2d.py](py/elliptic_solver_2d.py) explain this was a deliberate
authenticity choice — that tradeoff is wrong at nx ≈ 1000, so it becomes a
configurable behavior with the fast path as default.)

Steps:
1. Add a module-level plan cache in [py/_fftw_native.py](py/_fftw_native.py) keyed
   by transform shape: first request plans, subsequent requests reuse. Reference
   counting or "never free until process exit" are both acceptable; document which.
2. Add FFTW wisdom persistence: export wisdom to a cache file (respect an env var
   for the path) after first planning; import before planning. Second process run
   plans in ~0 time.
3. Make the planner rigor a config option (`exhaustive` | `measure` | `estimate`),
   default `measure`, with the old behavior available for authenticity comparisons.
4. Keep the change strictly inside `_fftw_native.py` — callers must not change.

**Acceptance:** force histories from both Phase-0 baselines are **bit-identical**
(plan choice never affects DST math up to FFTW's deterministic output for a given
plan; if bit-identity fails across plan algorithms, require agreement < 1e-13 and
document). Startup time at nx ≥ 1024, ngrid = 5 drops by >10×; show a before/after
timing table.

### 2.2 Vectorize the transfer operators

1. `Scalar.getBC` ([py/scalar.py:176](py/scalar.py:176)): replace the per-element
   list comprehensions (the `self(lev+1, …)` gathers for bottom/top/left/right even
   nodes) with direct numpy slicing on `self._data[lev+1]`. Mind the (1,1) index
   offset documented in the port notes (element (i,j) lives at `[i−1, j−1]`) and
   the boundary rows/columns that are *not stored* (value 0 on the coarse level's
   own outer boundary) — the current scalar `__call__` handles those implicitly;
   the vectorized version must reproduce zeros there explicitly.
2. Audit and vectorize any remaining per-point Python loops in
   [py/vector_operations.py](py/vector_operations.py) inter-level paths
   (`FluxToXVelocity`, `XVelocityToFlux`, `CrossProduct` overlap fill) and in
   `Scalar.coarsify`'s ghost/edge handling.
3. Do **not** restructure the algorithms — same stencils, same traversal, just
   array-at-a-time.

**Acceptance:** both baselines agree with Phase-0 references to < 1e-12 in force
histories (ideally bit-identical; summation-order changes may prevent that — if so,
state it). `getBC` disappears from the top-10 profile at nx = 1024. Plot the new
timing breakdown next to the old one.

### 2.3 Optional extras (only if still needed after 2.1–2.2)

- Enable FFTW threads for the per-level DSTs (worthwhile above ~1500²).
- The JAX/GPU port the codebase is annotated for (every file carries
  "JAX-readiness" notes) — the biggest absolute lever, but a separate project;
  do not fold it into this work.

**Effort estimate:** small-to-medium. High payoff, near-zero numerical risk.

---

## Phase 3 — Configurable coarsening ratio R ∈ {2, 4} (preset)

**Goal:** generalize the level-nesting ratio from the hard-coded 2 to an even
integer R, exposed as a new input knob (suggested name `gridratio`, default 2,
validated to {2, 4}). Odd ratios (3) are explicitly out of scope — they break the
node-parity structure and the shift quantization for marginal benefit.

**Ship with a startup warning when R = 4:** larger inter-level resolution jump;
validate wake quantities against an R = 2 run before trusting results.

### 3.1 Math spec

Let R be the ratio, fine spacing dx, coarse spacing R·dx. The fine domain occupies
the middle 1/R of the coarse domain (before shift).

- **Divisibility:** `nx % (2R) == 0` and `ny % (2R) == 0` (for R=2 this is the
  existing `% 4`). Reason: coarse cells outside the fine domain per side is
  `nx·(R−1)/(2R)`, which is an integer when `nx = 2R·m`.
- **NxExt** (coarse cells left of fine domain):
  `NxExt = nx·(R−1)/(2R) · (1 − xShift)`, rounded as the current code rounds.
  Check: R=2 gives `nx/4·(1−xShift)` — matches [grid.py:111](py/grid.py:111).
- **Shift quantization:** `xShift · nx · (R−1) mod 2R == 0`
  (R=2 reduces to the existing `xShift·nx % 4 == 0`).
- **Index maps:** `c2f: ii = (i − NxExt)·R`; `f2c: i = ii // R + NxExt`.
- **Spacing/offsets:** `Dx(lev) = dx · R^lev`; in `_getXOffset/_getYOffset` replace
  `(1 << lev) − 1` with `R^lev − 1` (the centered-nesting derivation
  `Σ L(R^{k+1}−R^k)/2 = L(R^lev−1)/2` reproduces the current formula at R=2,
  including the shift term).
- **Restriction (coarsify):** general full-weighting hat kernel. 1-D weights over
  fine offsets k ∈ {−(R−1), …, R−1}:  `w_k = (R − |k|) / R²`  (sums to 1).
  2-D stencil is the tensor product; for R=2 this is exactly the current
  1/4–1/8–1/16 9-point stencil, for R=4 it is a 7×7 stencil. Coarse interior nodes
  strictly inside the fine domain are `i = NxExt+1 … NxExt + nx/R − 1` (and
  analogously in y); each maps to fine node `ii = (i − NxExt)·R`.
  **Invariant: weights sum to exactly 1** (constant fields are preserved) — assert
  this in a unit test, not in production code.
- **Prolongation for BCs (getBC):** fine boundary node `ii` coincides with a coarse
  node when `ii % R == 0` (copy). For `ii = R·q + k`, `k = 1 … R−1`: linear
  interpolation `u = ((R−k)/R)·u_coarse(q) + (k/R)·u_coarse(q+1)`. R=2 reduces to
  the current copy-evens / average-odds. Same treatment on all four edges; corners
  come out consistent automatically because corner fine nodes have `ii % R == 0`
  at both extents (guaranteed by the divisibility rule).

### 3.2 File-by-file checklist

Work through in this order; run the R=2 regression after each file.

1. **[py/grid.py](py/grid.py)** — add `_ratio` (constructor arg, default 2, plumb
   through `resize`); update the `% 4` assert, `NxExt`/`NyExt`, `c2f`/`f2c`,
   `Dx(lev)`, `setXShift`/`setYShift` asserts, `_getXOffset`/`_getYOffset`,
   `isEqualTo`. Every `1 << lev` and `// 4` in this file must go.
2. **[py/scalar.py](py/scalar.py)** — `coarsify` (general kernel + generalized loop
   bounds), `getBC` (general prolongation; coordinate with the Phase-2.2
   vectorization — do Phase 2.2 first so this lands on vectorized code).
3. **[py/elliptic_solver.py](py/elliptic_solver.py)** — `dx = self._dx * (1 << lev)`
   → ratio power (line 86).
4. **[py/vector_operations.py](py/vector_operations.py)** — audit every inter-level
   path (`_curl_scalar_to_flux`'s `getBC` use is fine once getBC is general; check
   `FluxToXVelocity`/`YVelocityToFlux`/`CrossProduct` for literal `//2`, `[::2]`,
   halving of extents).
5. **[py/flux.py](py/flux.py), [py/bc.py](py/bc.py)** — expected unchanged (sized
   by nx/ny only); verify by grep, note the verification in the PR description.
6. **[py/state.py](py/state.py) / restart I/O** — the restart format stores grid
   parameters; add ratio with a backward-compatible default of 2 when absent.
   Old restart files must load unchanged.
7. **[py/ibpm.py](py/ibpm.py)** — new `gridratio` parameter next to lines 182–189,
   validated ∈ {2, 4}; echo it in the startup banner block (lines 271–278);
   pass to `Grid`.
8. **[py/grid_planner.py]** (from Phase 1) — accept R, apply the generalized
   divisibility and level-count formulas.
9. Re-run the **grep checklist** from §1 across the whole `py/` tree; every
   remaining hit must be justified.

### 3.3 Tests and acceptance

- **R=2 regression:** both Phase-0 baselines bit-identical (nothing about R=2
  semantics may change).
- **Unit — restriction:** coarsify of a constant field is exact at every level for
  R ∈ {2, 4}; coarsify of a linear field is exact away from domain edges.
- **Unit — prolongation:** getBC of a linear coarse field reproduces the linear
  function on the fine boundary exactly, R ∈ {2, 4}.
- **Integration — Poisson:** manufactured smooth RHS with compact support in the
  finest domain; solve with (R=2, ngrid=5) and (R=4, ngrid=3) chosen to give the
  same outer extent; compare against a single-domain fine reference. R=4 error may
  be larger but must converge at second order as dx is refined.
- **Physics — Kurtulus case:** one pitching NACA0012 case from `kurt_comp` run both
  ways at matched outer extent; report St, mean/peak CL, CD side by side **with the
  comparison plot**. Decide from the data whether R=4 is acceptable for wake work;
  record the verdict in the README next to the plot.
- **Performance:** timing table showing the ngrid reduction and net speedup, so the
  bounded-payoff claim in §0.1 is confirmed or corrected with measurements.

**Effort estimate:** medium. Mechanical but wide; the danger is a missed factor-2
site producing silently wrong far-field BCs — this is what the manufactured-solution
Poisson test exists to catch.

---

## Phase 4 (stretch, only if justified by data) — Per-level rectangle extents

The honest version of "shape-adaptive levels": let each level ℓ have its own
`(nx_ℓ, ny_ℓ)` rectangle (cells still square, ratio still R, every fine-box corner
constrained to lie on coarse-grid nodes with the parity the transfer operators
need). A thin finest strip around an airfoil inside progressively squarer coarse
levels captures nearly everything a U-shape would, while every level stays a
rectangle the DST solver can handle.

Consequences to be aware of before choosing to do this:

- `Scalar`/`Flux` storage changes from one `(ngrid, nx−1, ny−1)` array to a list of
  per-level arrays — touching indexing, `__getitem__`, arithmetic operators, every
  function in `vector_operations.py`, `InnerProduct` (energy weighting), restart
  I/O, and all Tecplot/probe output. This is the single largest refactor in this
  manual and invalidates the "same-shape every level" assumption baked into ~10
  files.
- `getBC`/`coarsify` need generalized offsets (`NxExt` becomes per-level and
  per-side) — the Phase-3 math carries over with `NxExt_ℓ,side` as free integers
  instead of the centered formula.
- The elliptic solver needs one DST plan per distinct level shape (plan cache from
  Phase 2.1 handles this for free).

**Gate:** only start Phase 4 if, after Phases 1–3, profiling shows the finest level
still dominated by cells far from the body/wake (e.g. > 60% of finest-level area
outside the region needing fine dx). Otherwise the added complexity is not paid for.

---

## Rejected ideas — do not implement, and why (for the record)

- **U-shaped (non-convex) levels:** incompatible with the DST-I fast solver
  (tensor-product rectangles only); all four transfer/BC code paths assume straight
  edges; the immersed body must be fully inside the finest region anyway, so the
  useful U-geometries are rare; a non-rectangular Poisson solve requires an
  iterative/domain-decomposition method that would cost more than the ≤ ~40% of one
  level's cells a U carve-out saves. Phase 4 is the sanctioned substitute.
- **Ratio 3 (or any odd ratio):** breaks even/odd node coincidence in `getBC`, the
  `% 4` family of asserts, and shift quantization; every parity-based slice
  (`[::2]`, `[1::2]`) would need three-phase equivalents. All cost, little gain
  over R=4.
- **Anisotropic cells (dx ≠ dy):** the DST solver and curl/Laplacian generalize
  easily, but the IB regularizer's discrete delta function assumes isotropic
  spacing (`grid.Dx()` is *the* spacing throughout
  [py/regularizer.py](py/regularizer.py)), and anisotropic cells degrade vortical
  accuracy exactly where this solver is used. Rectangular *domains* (already
  supported) give the shape benefit without this.

---

## Implementation order and definition of done

1. Phase 0 (baselines + profiler) — prerequisite for everything.
2. Phase 1 (planner) — independent, ships alone, immediate user value.
3. Phase 2 (FFTW + vectorization) — the real answer to "extremely small dx".
4. Phase 3 (ratio preset) — optional, after 2, behind `gridratio` with default 2.
5. Phase 4 — only if the Phase-4 gate condition is met, as its own project.

Done means: all acceptance tests above pass; R=2 default behavior bit-identical to
pre-change; every phase's PR includes the before/after timing table and the
required figures; the grep checklist comes back clean or annotated.
