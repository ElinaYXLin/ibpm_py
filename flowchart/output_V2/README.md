# V2: a "more user-friendly" execution flowchart

`main_execution_flowchart_v2.png` is built from a mentor's own draft outline
(`flowchart_draft.docx`), not from re-tracing the code from scratch. It's a
**pruned, 3-level-deep** version of [`../output/main_execution_flowchart.png`](../output/main_execution_flowchart.png)
("V1") — same underlying run (`py/ibpm.py`'s `main()`), but ~29 boxes
instead of V1's ~60+, and call depth capped at 3 layers instead of V1's 8.
Regenerate with `python3 flowchart/generate_execution_flowchart_v2.py`.

**Numbering** (matches the draft, and is *also* how deep each box goes into
the call stack): `1, 2, 3, ...` = a top-level step, in the actual source
order of `ibpm.py`'s `main()`. `a, b, c, ...` = one layer of function-call
depth to the right of its parent number. `i, ii, iii, ...` = one more layer
of call depth to the right of its parent letter. Every box still carries a
`file:line` citation (small text above the box), same convention as V1.

**Box text is the mentor's draft, verbatim** — including its original
wording and (where present) typos — **except** the six boxes that were
originally tagged `(unsure)` in the draft. Those six were corrected against
the actual source and now show the corrected explanation instead (the
`(unsure)` tag itself doesn't appear in the diagram). This file explains
each correction, since if the mentor was unsure here, other readers of the
diagram probably would be too.

## The 6 corrections (boxes originally marked "(unsure)")

### Box 4 — "Create base flow"

**Draft**: "Create base flow (unsure)."
**Corrected to**: "Create the free-stream (\"base\") flow, at the given
angle of attack."
**Why**: `BaseFlow(grid, magnitude, alpha)` ([ibpm.py:298](../../py/ibpm.py#L298))
isn't the flow field being solved for — it's the **prescribed uniform
free-stream flow** (magnitude 1, at the `-alpha` angle of attack) that gets
superimposed on the perturbation vorticity field the solver actually
computes. For unsteady cases (a pitching/plunging reference frame), this
same object also carries the base flow's own prescribed motion.

### Box 6 — "Initialize state and load initial conditions"

**Draft wording was already accurate** — kept close to verbatim, just
tightened. `x = State(grid, ...)` creates the state (vorticity ω, flux q,
boundary force f), zeros it, then loads it from `-ic <file>` if one was
given (otherwise zero IC); for a linearized run it can also subtract a base
flow to form a perturbation state, and for an SFD run it also loads the
filtered state. See [ibpm.py:361-369](../../py/ibpm.py#L361-L369).

### Box 7 — "Move bodies (elaborate on why)"

**Draft**: "Move bodies (unsure, elaborate on why we need to move bodies)."
**Corrected to**: "Move bodies to their position at the current time."
**Why**: the source has a comment directly above this call explaining the
ordering ([ibpm.py:358-360](../../py/ibpm.py#L358-L360)): *"still need to
initialize model, but wait until after loading the initial condition, so we
know what the initial time is, for moving the bodies."* In other words: the
geometry object was built from the static `.geom` file (a body's *reference*
shape/position), and if that body has a prescribed motion, its *actual*
position depends on time — but the actual starting time isn't known until
*after* the initial condition is loaded (an IC file can set a nonzero
starting time, and `-resettime` can zero it again). So bodies are moved to
their correct position for `x.time` right after the IC is settled, and
*before* the model/solver get initialized (steps 8-9), so those see the
correct initial body position.

### Box 8 — "Initialize the model"

**Draft**: "Initialize the model (unsure, what exactly is the model?)."
**Corrected to**: "Initialize the model." (text unchanged — the *box* was
already right; the confusion was about what "the model" refers to, answered
here instead of in the box itself, to keep the box short.)
**Why**: "the model" = the `NavierStokesModel` object built in step 5a. It
bundles the physics: the grid, the geometry, the Reynolds number, and the
base flow, and exposes the methods the solver calls every step
(`updateOperators`, `refreshState`, `getConstraints`, the nonlinear term
`N(x)`, etc.). `model.init()` itself
([navier_stokes_model.py:66-72](../../py/navier_stokes_model.py#L66-L72))
does one thing: the first time it's called, it triggers the regularizer's
initial setup (see box 10 below for what the regularizer is).

### Box 10 — "Update operators, regularizers, and move the flow"

**Draft**: "Update operators, regularizers, and move the flow (unsure, what
are operators? what is moving the flow?)."
**Corrected to**: "Update operators: re-sync the base flow / body positions
/ regularizer to the current time."
**Why**: `model.updateOperators(t)`
([navier_stokes_model.py:110-116](../../py/navier_stokes_model.py#L110-L116))
does up to two things, each conditional: (1) if the **base flow** is
time-dependent (an unsteady/pitching free-stream), it calls
`baseFlow.moveFlow(t)` — this is "moving the flow" the draft's question was
about, it's the base flow's own prescribed velocity being advanced to time
`t`, not the fluid velocity field. (2) if the **geometry** is time-dependent
(a moving body), it calls `geometry.moveBodies(t)` *and*
`regularizer.update()`. The "regularizer" is the interpolation/regularization
operator that maps between the fluid's Eulerian grid and the body's
Lagrangian boundary points (spreading forces from the boundary onto the
grid, and interpolating velocities from the grid onto the boundary) — it has
to be recomputed whenever the body's position changes, which is why it's
updated in lockstep with `moveBodies`.

### Box 15bii — "Calculate the 'a' and 'b' terms"

**Draft**: "Calculate a and b (constraints) terms (unsure, what are a and b
terms)."
**Corrected to**: "Calculate the \"a\" and \"b\" right-hand-side terms for
the constrained (projection) solve."
**Why**: both are named directly in the source's own comments
([ib_solver.py:190,205](../../py/ib_solver.py#L190)): *"Evaluate
Right-Hand-Side (a) for first equation of ProjectionSolver"* and
*"Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver."*
Concretely: **a** combines viscous diffusion (`Laplacian(ω)`, scaled by
time-scheme coefficients) with the convective/nonlinear term computed in
step 15a — the right-hand side of the semi-implicit vorticity update. **b**
is `model.getConstraints()`, whose own docstring
([navier_stokes_model.py:99-104](../../py/navier_stokes_model.py#L99-L104))
says it plainly: *"the velocity of the bodies minus the base flow
velocity"* — the no-slip constraint the boundary force has to satisfy that
step.

### Box 15bv — "If SFD solver, integrate filtered state ωhat"

**Draft**: "If SFDSolver, integrate filtered state ωhat (unsure, why?)."
**Corrected to**: text unchanged (already accurate) — the "why" is answered
here.
**Why**: "SFD" = Selective Frequency Damping
([ib_solver.py:324](../../py/ib_solver.py#L324)), a technique for finding
an unstable steady-state solution by feeding back a damping force
proportional to the gap between the instantaneous state and a low-pass
*filtered* version of it (`ωhat`): the nonlinear term picks up an extra
`-chi*(ω - ωhat)` term
([ib_solver.py:347-352](../../py/ib_solver.py#L347-L352), the `N` override
inside `class SFDSolver`). `ωhat` isn't a
fixed reference — it's defined by its own relaxation ODE,
`d(ωhat)/dt = (ω - ωhat)/Delta`, so it has to be time-integrated forward
right alongside the real state every substep — that's this extra step. It
only applies when `-model sfd` is selected; every other solver type skips
it entirely.

## Other things flagged, NOT auto-corrected (your call)

The draft's instruction was to fix "(unsure)" boxes automatically and flag
everything else for a decision, so these are reported but the boxes above
still show the draft's original text:

- **Boxes 5a/5b spell "Navier-Stokes" as "Naviar-Stokes"** (also in box 11's
  "Naviar-Stokes model"), and box 5b spells "Linearized" as "Linearied."
  Minor typos, not flagged as `(unsure)`, so left as-is in the diagram.
- **Box 17 merges two calls with different cadence.** The draft's text is
  "Use logger to output results, and clean up the logger," placed as the
  last step of "SECTION 4: solver loop." In the source, though, this is
  actually two separate calls: `logger.doOutput(q_potential, x)`
  ([ibpm.py:462](../../py/ibpm.py#L462)) runs **every iteration** of the
  step loop (it's what "output results" refers to), while
  `logger.cleanup()` ([ibpm.py:480](../../py/ibpm.py#L480)) runs **once**,
  *after* the loop finishes entirely. Bundling them into one box makes it
  look like both happen every step. The diagram's loop-back arrow points
  from box 17 back to box 15, which is only strictly accurate for the
  `doOutput` half of box 17 — the `cleanup` half happens on the way out,
  after the last iteration. Splitting box 17 into two ("17. per-step:
  output results" inside the loop, "18. once, after the loop: clean up the
  logger" outside it) would remove this ambiguity, but that changes your
  numbering scheme, so it's flagged here rather than done automatically.
- **Section 4 has no explicit box for the loop header itself**
  (`for i in 1..numSteps:`, [ibpm.py:449](../../py/ibpm.py#L449)) — the
  section title text carries that instead. Consistent with the rest of the
  draft's pruning choices, just noting it since V1 does have a dedicated
  box for it.

## What's simplified relative to V1 (by design, not an error)

- No `ProjectionSolver.solve()` 6-step internal breakdown under box 15biii
  (V1 expands this one layer further — see V1's "Ainv/C/Minv/B/Ainv/assign"
  chain).
- No `createAllSolvers()`/`CholeskySolver`/`ConjugateGradientSolver`
  breakdown under box 5b, and no `NavierStokesModel.__init__` breakdown
  under box 5a.
- No per-output-type breakdown under boxes 14a (V1 shows
  `OutputTecplot`/`OutputRestart`/`OutputForce`/`OutputEnergy` as four
  separate boxes; V2 shows the one loop that dispatches to all of them).
- No separate `geometry.moveBodies(t) -> RigidBody.moveBody(t)` expansion
  under box 7 (V2 states the reason in prose here instead, per the
  "(unsure)" correction above).

None of these are wrong to omit — they're exactly the detail V1 already
has, for anyone who wants it. V2 is meant to be read on its own as a
simpler map of the same run.
