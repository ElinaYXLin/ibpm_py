# IBPM Python Port — Porting & Validation Guide

## Quick Start: Validation Workflow

### 1. Generate Reference Data (Already Done)

The C++ binary has been run and reference restart files are in `examples/`:

```bash
ls -lh examples/ibpm*.bin     # 11 files: timesteps 0-200
ls -lh examples/ibpm*.force   # Force data
```

These represent a full simulation (200 timesteps, dt=0.01, Re=100) on a 200×200 grid with a cylinder.

### 2. Understand the Binary Format

See `VALIDATION_SETUP.md` for the complete format specification.

TL;DR: Each `.bin` file contains:
- Grid parameters (6 numbers)
- Flux q: (nx+1)×ny X-fluxes + nx×(ny+1) Y-fluxes
- Vorticity ω: interior points only (nx-1)×(ny-1)
- Forces f: numPoints × 2 boundary forces
- Timestep and time

### 3. Read & Validate Files in Python

```python
from validate_restart_reader import IBPMRestart, compare_restarts

# Load a file
r = IBPMRestart('examples/ibpm00000.bin')
r.info()  # Print summary

# Extract data
q_x, q_y = r.get_flux_xy()       # (ngrid, nx+1, ny), (ngrid, nx, ny+1)
omega = r.get_omega()            # (ngrid, nx, ny) with zero boundaries
forces = r.get_forces()          # (numPoints, 2)

# Compare two files
result = compare_restarts('examples/ibpm00000.bin', 'examples/ibpm00100.bin', tol=1e-10)
print(f"Match: {result['match']}")  # Should be False (data changed)
```

## Porting Strategy: 5 Stages

Each stage builds on the previous and validates against reference data.

### Stage 1: Data Structures (Week 1)

**Goal**: Port Grid, State, Scalar, Flux classes. No computation yet.

**Validation**:
1. Create a Python State, load reference data into it
2. Save to binary, re-read, verify byte-for-byte match

**Files to port**:
- `src/Grid.h/cc` → `python/ibpm/grid.py`
- `src/Scalar.h/cc` → `python/ibpm/scalar.py`
- `src/Flux.h/cc` → `python/ibpm/flux.py`
- `src/State.h/cc` → `python/ibpm/state.py`
- `src/BoundaryVector.h/cc` → `python/ibpm/boundary.py`

**Test example**:
```python
from ibpm import Grid, State, Scalar, Flux, BoundaryVector
from validate_restart_reader import IBPMRestart

# Load reference
ref = IBPMRestart('examples/ibpm00000.bin')

# Create Python state
grid = Grid(nx=200, ny=200, ngrid=1, length=4.0, xOffset=-2.0, yOffset=-2.0)
state = State(grid, numPoints=160)

# Load reference data
state.omega = ref.get_omega()[0]
q_x, q_y = ref.get_flux_xy()
state.q = Flux.from_xy(q_x[0], q_y[0])  # (to be implemented)
state.f = BoundaryVector(ref.get_forces())
state.time = ref.data['time']
state.timestep = ref.data['timestep']

# Save and re-read
state.save('test_state.bin')
state2 = State.load('test_state.bin')

# Verify
assert np.allclose(state.omega, state2.omega)
assert np.allclose(state.get_forces(), state2.get_forces())
print("✓ State I/O validation passed")
```

### Stage 2: Basic Linear Algebra (Week 2)

**Goal**: Port Regularizer, EllipticSolver, ProjectionSolver.

**Validation**:
1. Load reference flux/vorticity
2. Apply operation, compare to reference

**Files to port**:
- `src/Regularizer.h/cc` → `python/ibpm/regularizer.py`
- `src/EllipticSolver2d.h/cc` → `python/ibpm/elliptic_solver.py` (uses scipy.fft.dst)
- `src/ProjectionSolver.h/cc` → `python/ibpm/projection_solver.py`

**Test example**:
```python
# Load reference initial and final states
ref0 = IBPMRestart('examples/ibpm00000.bin')
ref1 = IBPMRestart('examples/ibpm00001.bin')

# Create Python solver components
grid = Grid(...)
geom = Geometry.load('examples/cylinder.geom')

# Test regularizer: force spreading
regularizer = Regularizer(grid, geom)
forces = ref0.get_forces()
flux_from_f = regularizer.to_flux(forces)
# Compare to reference (no reference for this yet — need to compute numerically)

# Test elliptic solver: vorticity → streamfunction
omega_ref = ref0.get_omega()[0]
poisson = PoissonSolver(grid)
psi = poisson.solve(omega_ref)
# Check that curl(psi) ≈ omega, compare to reference flux

print("✓ Linear algebra validation passed")
```

### Stage 3: Geometry & I/O (Week 2)

**Goal**: Port Geometry, Motion, Output classes.

**Validation**:
1. Load geometry from .geom file
2. Verify point coordinates match expectations
3. Output Tecplot format (visual inspection)

**Files to port**:
- `src/Geometry.h/cc` → `python/ibpm/geometry.py`
- `src/RigidBody.h/cc` → `python/ibpm/rigid_body.py`
- `src/Motion*.h` → `python/ibpm/motion.py`
- `src/OutputTecplot.h/cc` → `python/ibpm/output.py`

**Test example**:
```python
geom = Geometry.load('examples/cylinder.geom')
print(f"Loaded {geom.getNumBodies()} bodies")
print(f"Total points: {geom.getNumPoints()}")  # Should be 160

# Verify point coordinates are reasonable
points = geom.getPoints()
assert points.shape == (160, 2)
assert np.all(np.sqrt(points[:, 0]**2 + points[:, 1]**2) > 0.49)
assert np.all(np.sqrt(points[:, 0]**2 + points[:, 1]**2) < 0.51)
print("✓ Geometry validation passed")
```

### Stage 4: NavierStokes Model (Week 3)

**Goal**: Port NavierStokesModel with B/C operators.

**Validation**:
1. Load reference forces and vorticity
2. Apply B operator: should produce correct forcing term
3. Apply C operator: should produce correct boundary velocities

**Files to port**:
- `src/NavierStokesModel.h/cc` → `python/ibpm/navier_stokes.py`
- `src/BaseFlow.h/cc` → `python/ibpm/base_flow.py`

**Test example**:
```python
model = NavierStokesModel(grid, geom, Reynolds=100)

# Test B operator
f_ref = ref.get_forces()
omega_from_f = model.B(f_ref)
# Verify this is the negative divergence of the stress tensor

# Test C operator
q_ref = ref.get_flux()[0]  # Level 0
f_from_q = model.C(q_ref)
# Verify this interpolates flux to boundary velocities

print("✓ NavierStokes model validation passed")
```

### Stage 5: Time Integration (Week 4)

**Goal**: Port IBSolver, Scheme, ProjectionSolver coupling.

**Validation**:
1. Load reference state at t=0
2. Advance one timestep
3. Compare to reference at t=dt
4. Iteratively advance N steps, check forces and energy

**Files to port**:
- `src/IBSolver.h/cc` → `python/ibpm/ib_solver.py`
- `src/Scheme.h` → `python/ibpm/scheme.py`

**Test example**:
```python
# Load reference
ref0 = IBPMRestart('examples/ibpm00000.bin')
ref1 = IBPMRestart('examples/ibpm00001.bin')

# Create Python solver
model = NavierStokesModel(grid, geom, Reynolds=100)
solver = NonlinearIBSolver(grid, model, dt=0.02, scheme='RK3')
solver.init()

# Load initial state
state = State(grid, geom.getNumPoints())
state.omega = ref0.get_omega()[0]
q_x0, q_y0 = ref0.get_flux_xy()
state.q = Flux.from_xy(q_x0[0], q_y0[0])
state.f = ref0.get_forces()
state.timestep = 0
state.time = 0.0

# Advance one step
solver.advance(state)

# Compare to reference
q_x_ref, q_y_ref = ref1.get_flux_xy()
omega_ref = ref1.get_omega()[0]

q_x_out, q_y_out = state.q.to_xy()
assert np.allclose(state.omega, omega_ref, rtol=1e-8), "Vorticity mismatch"
assert np.allclose(q_x_out, q_x_ref, rtol=1e-8), "Q_x mismatch"
assert np.allclose(q_y_out, q_y_ref, rtol=1e-8), "Q_y mismatch"

print(f"✓ Single timestep validation passed")
print(f"  Time: {state.time:.6f} (expected 0.02)")
print(f"  Timestep: {state.timestep} (expected 1)")

# Full integration test
for i in range(10):
    solver.advance(state)
    ref_final = IBPMRestart(f'examples/ibpm00{(i+1)*20:03d}.bin')
    # ... compare ...
```

## Reference Data Structure

```
examples/
├── ibpm00000.bin          # t=0, step=0
├── ibpm00001.bin          # t=0.02, step=1
├── ibpm00002.bin          # t=0.04, step=2
├── ...
├── ibpm00200.bin          # t=4.0, step=200
├── ibpm.force             # Drag/lift vs time
└── cylinder.geom          # Geometry file
```

**Grid parameters** (constant across all files):
- nx=200, ny=200 (grid cells)
- ngrid=1 (single domain)
- dx=0.02 (grid spacing)
- x0=-2.0, y0=-2.0 (domain origin)
- Length=4.0 (domain width)

**Simulation parameters**:
- dt=0.02 (timestep)
- Re=100 (Reynolds number)
- Scheme=RK3 (3-step Runge-Kutta)
- Model=nonlinear
- 160 boundary points (cylinder)

## Testing Utilities

### Run Validation Tests
```bash
python3 test_validation_harness.py        # Run all tests
python3 test_validation_harness.py summary # Print data summary
```

### Compare Two Files
```python
from validate_restart_reader import compare_restarts
result = compare_restarts('examples/ibpm00000.bin', 'examples/ibpm00100.bin', tol=1e-10)
print(result)
```

### Extract and Plot Data
```python
from validate_restart_reader import IBPMRestart
import matplotlib.pyplot as plt

r = IBPMRestart('examples/ibpm00100.bin')
omega = r.get_omega()[0]

plt.figure(figsize=(10, 8))
plt.contourf(omega, levels=20)
plt.colorbar(label='Vorticity')
plt.title(f'Vorticity at t={r.data["time"]:.2f}')
plt.savefig('vorticity.png')
```

## Common Pitfalls

1. **Index mismatch**: C++ uses 1-indexed interior loops; Python should use 0-indexed or be explicit about padding
2. **Flux layout**: X-fluxes and Y-fluxes are stored sequentially, not interleaved
3. **Omega boundaries**: Only interior points are stored in binary; boundaries are implicitly zero
4. **Array copies**: Assignments create copies; use in-place operations where performance matters
5. **FFT normalization**: SciPy's DCT has different conventions than FFTW; check carefully

## Performance Notes

- Initial port should prioritize correctness, not speed
- Use NumPy/SciPy for FFTW (scipy.fft.dst) and linear algebra (scipy.linalg, scipy.sparse)
- Avoid Python loops for grid operations; vectorize with NumPy
- Cholesky factorization should use scipy.linalg.cholesky, not manual implementation

## Next Steps

1. **Choose starting point**: Stage 1 (data structures) is lowest risk
2. **Set up testing**: Each ported module gets a test file comparing to reference
3. **Iterate**: Port one component, validate, commit, repeat
4. **Document decisions**: Record why C++ patterns were translated a particular way

