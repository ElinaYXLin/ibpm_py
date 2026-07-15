# IBPM Python port

This is a Python port of the Immersed Boundary Projection Method (IBPM)
code in the parent directory (`src/`). It solves the 2D incompressible
Navier-Stokes equations around complex geometries using the fast
immersed-boundary projection method of Taira & Colonius (2007) / Colonius &
Taira (2008), with the same multi-domain far-field boundary condition
scheme as the original C++ code.

The port is line-by-line faithful to the C++ source: class names, method
names, and file layout mirror `src/*.h` / `src/*.cc` as closely as Python
allows (see the `NOTE(port)` comments sprinkled throughout for the handful
of places where C++ and Python semantics genuinely diverge). Nothing here
has been "improved" or restructured relative to the original algorithm.

## Requirements

- Python 3.9+
- `numpy`
- `scipy` (used by `elliptic_solver_2d.py`, for the sine-transform Poisson solve)

Install both with:

```bash
pip install -r py/requirements.txt
```

No compilation step is needed — unlike the C++ version, there's nothing to
`make`. Just run the scripts directly with Python.

## File organization

All ~45 modules live flat inside `py/`, one file per class/topic, e.g.
`py/state.py` <-> `src/State.h`/`src/State.cc`. This mirrors `src/`, which
is *also* one flat directory of ~40 files with no subfolders — so the flat
layout here isn't a Python-specific choice, it's carried over from the
original C++ project's own convention.

We considered pulling the handful of general-purpose helper modules (e.g.
`utils.py`, `parm_parser.py`, `logger.py`, `scheme.py`) into a `py/utils/`
subfolder, but decided against it:

- **It breaks the 1:1 mapping to `src/`.** Right now, given a C++ file you
  can find its Python port by name alone (`Grid.h` -> `grid.py`). Moving
  files into subfolders would break that correspondence for no benefit,
  making it harder to cross-check the port against the original.
- **There's no clean "utility vs. core" line to draw.** Files like
  `vector_operations.py` or `elliptic_solver.py` are both "core physics"
  *and* general-purpose numerical utilities used all over the codebase;
  a `utils/` bucket would end up being an arbitrary, ill-defined split
  rather than a meaningful grouping.
- **It's not solving an actual problem.** All imports are already
  explicit (`from .grid import Grid`), so there's no import-cycle or
  namespace-collision issue a subfolder would fix — 45 files in one
  directory is still easy to navigate/search, and the original C++ project
  manages 80 files in `src/` the same way.
- **It has real, non-trivial cost.** Every relative import
  (`from .x import Y`) in every moved file, and every file that imports
  *from* a moved file, would need updating. That's a mechanical but
  error-prone change across dozens of files, with real risk of quietly
  breaking the cross-validation tests in `py/tests/cross_validation/`,
  for a purely cosmetic reorganization.

If the project grows substantially beyond its current scope (e.g. a real
package split with a public API vs. internals), it may be worth revisiting
this — but for a faithful, already-flat 1:1 port of a flat C++ codebase,
introducing a `utils/` subfolder now would add churn and risk without a
concrete benefit.

## Which file do I run?

There are two runnable programs, both with a `main()` function, mirroring
the two C++ executables in `src/`:

| To do this...                                    | Run this                  | C++ equivalent    |
|---------------------------------------------------|----------------------------|--------------------|
| Run a full flow simulation                         | `py/ibpm.py`               | `src/ibpm.cc`       |
| Check/preview a geometry file before simulating it | `py/checkgeom.py`          | `src/checkgeom.cc`  |

**Start with `checkgeom.py`** to sanity-check your geometry file, then use
`ibpm.py` to actually run the simulation.

Both are invoked as Python modules from the repository root (the directory
that contains `py/`), **not** by running the `.py` file directly, since
they rely on relative imports within the `py` package:

```bash
cd /path/to/ibpm-master

# check a geometry file and write a Tecplot preview of it
python3 -m py.checkgeom -geom examples/cylinder.geom -o check.plt

# run a simulation
python3 -m py.ibpm -geom examples/cylinder.geom
```

Every option has a `-h` flag that prints full usage, exactly like the C++
binaries did:

```bash
python3 -m py.ibpm -h
python3 -m py.checkgeom -h
```

## Running a simulation (`py.ibpm`)

The most basic invocation just points at a geometry file:

```bash
python3 -m py.ibpm -geom examples/cylinder.geom
```

This uses all the defaults: a 200x200 grid, zero initial condition, 250
timesteps, the nonlinear model, and an RK3 time-stepping scheme. It writes
its output into the current directory:

- `ibpm00100.bin`, `ibpm00200.bin`, ... — binary restart files (readable
  again by `py.ibpm -ic ...`, see below)
- `ibpm00100.plt`, `ibpm00200.plt`, ... — ASCII Tecplot files with `u`,
  `v`, and vorticity fields, for visualization
- `ibpm.force` — lift/drag history, one line per timestep
- `ibpm.cmd` — the exact command-line arguments used, so the run can be
  reproduced later
- `*.cholesky` — cached Cholesky factorization of the projection operator
  (computed once, reused on subsequent runs with the same grid/geometry/dt)

Some of the most commonly-adjusted options:

| Flag          | Meaning                                              | Default |
|---------------|-------------------------------------------------------|---------|
| `-geom`       | geometry file to load                                 | `<name>.geom` |
| `-name`       | run name, used as a prefix for all output files        | `ibpm`  |
| `-outdir`     | directory to write output into                         | `.`     |
| `-nx -ny`     | grid resolution                                         | `200 200` |
| `-ngrid`      | number of nested grid levels (multi-domain far field)  | `1`     |
| `-Re`         | Reynolds number                                          | `100`   |
| `-dt`         | timestep                                                | `0.02`  |
| `-nsteps`     | number of timesteps to run                              | `250`   |
| `-scheme`     | `euler`, `ab2`, `rk3`, or `rk3b`                        | `rk3`   |
| `-model`      | `nonlinear`, `linear`, `adjoint`, `linearperiodic`, `sfd` | `nonlinear` |
| `-ic`         | restart file to use as the initial condition (instead of zero) | (none) |
| `-tecplot`    | write a Tecplot file every N steps (`0` = never)       | `100`   |
| `-restart`    | write a restart file every N steps (`0` = never)       | `100`   |
| `-force`      | write forces every N steps (`0` = never)                | `1`     |
| `-energy`     | write kinetic energy every N steps (`0` = never)        | `0`     |

Run `python3 -m py.ibpm -h` for the complete list (there are also options
for angle of attack, unsteady base flow, SFD stabilization parameters, and
the linear/adjoint/periodic model variants).

## Choosing a geometry

Geometries are plain-text `.geom` files — see `examples/cylinder.geom` for
the simplest possible example:

```
# Cylinder, diameter 1, with 160 points

body Cylinder
  circle_n 0 0 0.5 160
end
```

A `.geom` file contains one or more `body ... end` blocks. Inside a body
block, you build up the boundary points with one or more of:

| Command                                   | Meaning |
|--------------------------------------------|---------|
| `circle x y radius dx`                     | circle centered at `(x,y)`, spaced roughly `dx` apart |
| `circle_n x y radius n`                    | same, but with exactly `n` points |
| `line x0 y0 x1 y1 dx`                      | straight segment, spaced roughly `dx` apart |
| `line_n x0 y0 x1 y1 n`                     | same, but with exactly `n` points |
| `line_aoa length xC yC angle n`            | a line of given length/angle-of-attack about center `(xC,yC)`, with `n` points |
| `point x y`                                | a single point |
| `raw filename`                             | load a raw list of `x y` points (one pair per line) from `filename` |
| `center x y`                               | set the body's center of rotation |
| `name ...`                                 | give the body a name |
| `motion <type> <params...>`                | attach a prescribed motion (see below) |

You can combine several of these within one `body` block to build up more
complex shapes, and you can have multiple `body ... end` blocks in one file
for multi-body geometries.

To use your own geometry, either point `-geom` at a new file:

```bash
python3 -m py.ibpm -geom my_airfoil.geom
```

or edit `examples/cylinder.geom` directly / save a copy under a new name.
**Always check a new geometry with `py.checkgeom` first** — it's much
cheaper than a full run, and lets you catch a too-sparse or malformed
boundary before wasting time on a simulation that won't converge:

```bash
python3 -m py.checkgeom -geom my_airfoil.geom -nx 200 -ny 200 -o my_airfoil_check.plt
```

Open the resulting `.plt` file in Tecplot (or another Tecplot-ASCII-capable
viewer) and look at the boundary point spacing relative to the background
grid — if the points look "leaky" (spaced farther apart than the grid
cells), increase the point count (`circle_n`'s last argument, or decrease
`dx` for `circle`/`line`).

### Moving bodies

To simulate a body in motion (pitching, plunging, etc.), add a `motion`
line inside the `body` block, e.g.:

```
body Airfoil
  line_aoa 1.0 0.0 0.0 0.0 200
  motion pitchplunge 0.1 1.0 0.0 0.1 1.0 1.5707963
end
```

Supported motion types: `fixed`, `fixedvel`, `pitchplunge`,
`sigmoidalstep`, `lagstep1`, `lagstep2`, `eldredge`, `eldredgecombined2`,
`eldredge1`, `eldredge2`, `motionfile`, `motionfileperiodic` — each takes a
different number of numeric parameters; see the corresponding
`py/<motion_name>.py` file's docstring/`__init__` for what they mean (they
mirror `src/<MotionName>.h` exactly).

## Known issue

`py/cholesky_solver.py` (ported in an earlier session) currently has a bug
that can make `solver.init()` fail with `AssertionError: assert s > 0` for
some grid/geometry/Re/dt combinations. If you hit this, it's a known,
already-flagged issue — not something wrong with your command line or
geometry file.
