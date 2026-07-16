# V2: a "more user-friendly" execution flowchart

`main_execution_flowchart_v2.png` is built from a draft outline
(`flowchart_draft.docx`), not from re-tracing the code from scratch. It's a
**pruned, 3-level-deep** version of [`../output/main_execution_flowchart.png`](../output/main_execution_flowchart.png)
("V1") — same underlying run (`py/ibpm.py`'s `main()`), but ~29 boxes
instead of V1's ~60+, and call depth capped at 3 layers instead of V1's 8.
Regenerate with `python3 flowchart/generate_execution_flowchart_v2.py`.

**Depth is shown by column position and rightward arrows only** — no
number/letter/roman-numeral labels on the boxes. Column 0 is a top-level
step, in the actual source order of `ibpm.py`'s `main()`. Column 1 is one
layer of function-call depth to the right of its column-0 parent. Column 2
is one more layer to the right of its column-1 parent. Every box carries a
`file:line` citation, positioned close above the box it cites, with enough
clearance that it never overlaps the box placed above it in the same
column.

**Box text is the draft's own wording, verbatim** — **except** the six
boxes that were originally tagged `(unsure)` in the draft. Those six were
corrected against the actual source and now show the corrected explanation
instead (the `(unsure)` tag itself doesn't appear in the diagram). This
file explains each correction, since if you were unsure here, other readers
of the diagram probably would be too. A few other inaccuracies were found
in the rest of the draft (not marked `(unsure)`) and have also been fixed
in the diagram — listed at the end of this file.

## The 6 corrections (boxes originally marked "(unsure)")

### "Create the free-stream ('base') flow, at the given angle of attack"

**Draft**: "Create base flow (unsure)."
**Why corrected**: `BaseFlow(grid, magnitude, alpha)`
([ibpm.py:298](../../py/ibpm.py#L298)) isn't the flow field being solved
for — it's the **prescribed uniform free-stream flow** (magnitude 1, at the
`-alpha` angle of attack) that gets superimposed on the perturbation
vorticity field the solver actually computes. For unsteady cases (a
pitching/plunging reference frame), this same object also carries the base
flow's own prescribed motion.

### "Initialize state and load initial conditions"

**Draft wording was already accurate** — kept close to verbatim, just
tightened. `x = State(grid, ...)` creates the state (vorticity ω, flux q,
boundary force f), zeros it, then loads it from `-ic <file>` if one was
given (otherwise zero IC); for a linearized run it can also subtract a base
flow to form a perturbation state, and for an SFD run it also loads the
filtered state. See [ibpm.py:361-369](../../py/ibpm.py#L361-L369).

### "Move bodies to their position at the current time"

**Draft**: "Move bodies (unsure, elaborate on why we need to move bodies)."
**Why corrected**: the source has a comment directly above this call
explaining the ordering
([ibpm.py:358-360](../../py/ibpm.py#L358-L360)): *"still need to initialize
model, but wait until after loading the initial condition, so we know what
the initial time is, for moving the bodies."* In other words: the geometry
object was built from the static `.geom` file (a body's *reference*
shape/position), and if that body has a prescribed motion, its *actual*
position depends on time — but the actual starting time isn't known until
*after* the initial condition is loaded (an IC file can set a nonzero
starting time, and `-resettime` can zero it again). So bodies are moved to
their correct position for `x.time` right after the IC is settled, and
*before* the model/solver get initialized, so those see the correct initial
body position.

### "Initialize the model"

**Draft**: "Initialize the model (unsure, what exactly is the model?)."
Box text is unchanged — it was already right; the confusion was about what
"the model" refers to, answered here instead of in the box, to keep the box
short.
**What "the model" is**: the `NavierStokesModel` object built in the "Build
model and solver" step. It bundles the physics: the grid, the geometry, the
Reynolds number, and the base flow, and exposes the methods the solver
calls every step (`updateOperators`, `refreshState`, `getConstraints`, the
nonlinear term `N(x)`, etc.). `model.init()` itself
([navier_stokes_model.py:66-72](../../py/navier_stokes_model.py#L66-L72))
does one thing: the first time it's called, it triggers the regularizer's
initial setup (see the next correction for what the regularizer is).

### "Update operators: re-sync the base flow / body positions / regularizer to the current time"

**Draft**: "Update operators, regularizers, and move the flow (unsure, what
are operators? what is moving the flow?)."
**Why corrected**: `model.updateOperators(t)`
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

### "Calculate the 'a' and 'b' right-hand-side terms for the constrained (projection) solve"

**Draft**: "Calculate a and b (constraints) terms (unsure, what are a and b
terms)."
**Why corrected**: both are named directly in the source's own comments
([ib_solver.py:190,205](../../py/ib_solver.py#L190)): *"Evaluate
Right-Hand-Side (a) for first equation of ProjectionSolver"* and
*"Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver."*
Concretely: **a** combines viscous diffusion (`Laplacian(ω)`, scaled by
time-scheme coefficients) with the convective/nonlinear term computed one
step earlier — the right-hand side of the semi-implicit vorticity update.
**b** is `model.getConstraints()`, whose own docstring
([navier_stokes_model.py:99-104](../../py/navier_stokes_model.py#L99-L104))
says it plainly: *"the velocity of the bodies minus the base flow
velocity"* — the no-slip constraint the boundary force has to satisfy that
step.

### "If SFD solver: also integrate the filtered state ωhat"

**Draft**: "If SFDSolver, integrate filtered state ωhat (unsure, why?)."
Box text is unchanged — it was already accurate; the "why" is answered
here.
**Why**: "SFD" = Selective Frequency Damping
([ib_solver.py:324](../../py/ib_solver.py#L324)), a technique for finding
an unstable steady-state solution by feeding back a damping force
proportional to the gap between the instantaneous state and a low-pass
*filtered* version of it (`ωhat`): the nonlinear term picks up an extra
`-chi*(ω - ωhat)` term
([ib_solver.py:347-352](../../py/ib_solver.py#L347-L352), the `N` override
inside `class SFDSolver`). `ωhat` isn't a fixed reference — it's defined by
its own relaxation ODE, `d(ωhat)/dt = (ω - ωhat)/Delta`, so it has to be
time-integrated forward right alongside the real state every substep —
that's this extra step. It only applies when `-model sfd` is selected;
every other solver type skips it entirely.

## Other errors found and fixed (not marked "(unsure)" in the draft)

- **"Navier-Stokes" / "Linearized" typos**: the draft consistently spelled
  these correctly in the boxes that ended up in the final wording used
  here, so no change was needed there in the end.
- **The per-step output call and the one-time cleanup call were bundled
  into a single box under the loop-back arrow.** The draft's original text
  was "Use logger to output results, and clean up the logger," placed as
  the last step of "SECTION 4: solver loop." In the source, though, this is
  actually two separate calls with different cadence:
  `logger.doOutput(q_potential, x)` ([ibpm.py:462](../../py/ibpm.py#L462))
  runs **every iteration** of the step loop, while `logger.cleanup()`
  ([ibpm.py:480](../../py/ibpm.py#L480)) runs **once**, *after* the loop
  finishes entirely. Bundling them made it look like both happen every
  step, and made the loop-back arrow's target ambiguous. **Fixed**: split
  into two boxes — "Output results for this step" (462, inside the loop,
  the loop-back arrow now correctly wraps around it) and "Clean up the
  logger (once, after the loop)" (480, drawn after the loop-back, outside
  the repeated cycle).
- **No explicit box for the loop header itself**
  (`for i in 1..numSteps:`, [ibpm.py:449](../../py/ibpm.py#L449)) — only
  the section title carried that context, so the loop-back arrow had
  nowhere correct to point. **Fixed**: added a "For each step, 1 to
  numSteps:" box at the top of Section 4; the loop-back arrow now targets
  it directly, matching how V1 handles the same loop.

## What's simplified relative to V1 (by design, not an error)

- No `ProjectionSolver.solve()` 6-step internal breakdown under "Solve
  using the projection solver, given a and b" (V1 expands this one layer
  further — see V1's "Ainv/C/Minv/B/Ainv/assign" chain).
- No `createAllSolvers()`/`CholeskySolver`/`ConjugateGradientSolver`
  breakdown under "Build a solver...", and no `NavierStokesModel.__init__`
  breakdown under "Build a Navier-Stokes model".
- No per-output-type breakdown under "Call each entry that needs to be
  outputted..." (V1 shows `OutputTecplot`/`OutputRestart`/`OutputForce`/
  `OutputEnergy` as four separate boxes; V2 shows the one loop that
  dispatches to all of them).
- No separate `geometry.moveBodies(t) -> RigidBody.moveBody(t)` expansion
  under "Move bodies..." (V2 states the reason in prose here instead, per
  the correction above).

None of these are wrong to omit — they're exactly the detail V1 already
has, for anyone who wants it. V2 is meant to be read on its own as a
simpler map of the same run.
