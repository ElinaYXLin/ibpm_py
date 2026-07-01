# IBPM Validation & Testing Framework

## Overview

This directory contains tools for validating the Python port of the IBPM (Immersed Boundary Projection Method) C++ solver. The validation framework is built on **binary restart file comparison** — we run the C++ solver as ground truth and validate each ported component against the saved restart files.

## Files

| File | Purpose |
|------|---------|
| `validate_restart_reader.py` | Python reader for IBPM binary restart files; comparison utilities |
| `test_validation_harness.py` | Test harness with pre-loaded reference data and validation tests |
| `VALIDATION_SETUP.md` | Complete specification of binary file format |
| `PORTING_GUIDE.md` | Step-by-step guide for porting and validating each component |
| `examples/ibpm*.bin` | Reference restart files from C++ solver (11 files: timesteps 0–200) |
| `examples/cylinder.geom` | Geometry file (cylinder) |

## Quick Start

### 1. Verify Reference Data
```bash
python3 test_validation_harness.py summary
```

Output shows:
- Grid configuration (200×200, dx=0.02)
- 11 timesteps (t=0 to t=4.0, dt=0.02)
- Data ranges (flux, vorticity, forces)

### 2. Read a Restart File
```bash
python3 validate_restart_reader.py examples/ibpm00000.bin
```

Output:
```
Restart file: examples/ibpm00000.bin
  Grid: nx=200, ny=200, ngrid=1
  ...
  Flux stats:
    q_x range: [0.02, 0.02]      # Uniform background flow
    q_y range: [0.0, 0.0]
  Omega stats:
    range: [0.0, 0.0]             # Zero at t=0
  Forces stats:
    range: [0.0, 0.0]
```

### 3. Load Data in Python
```python
from validate_restart_reader import IBPMRestart
import numpy as np

# Load reference
ref = IBPMRestart('examples/ibpm00000.bin')

# Extract components
print(f"Timestep: {ref.data['timestep']}")
print(f"Time: {ref.data['time']:.6f}")

# Get flux (x and y components separately)
q_x, q_y = ref.get_flux_xy()
print(f"Q_x shape: {q_x.shape}")  # (ngrid, nx+1, ny) = (1, 201, 200)

# Get vorticity
omega = ref.get_omega()
print(f"Omega shape: {omega.shape}")  # (ngrid, nx, ny) = (1, 200, 200)

# Get forces
forces = ref.get_forces()
print(f"Forces shape: {forces.shape}")  # (numPoints, 2) = (160, 2)
```

### 4. Compare Two Files
```python
from validate_restart_reader import compare_restarts

result = compare_restarts(
    'examples/ibpm00000.bin',
    'examples/ibpm00100.bin',
    tol=1e-10
)

if result['match']:
    print("Files are identical")
else:
    print("Differences found:")
    for key, val in result['differences'].items():
        print(f"  {key}: {val}")
```

## Validation Workflow for Python Port

### 1. Before Porting Anything
```bash
# Verify reference data is valid
python3 test_validation_harness.py
```

Expected output: All tests pass ✓

### 2. After Porting Data Structures (Grid, State, Flux)
```python
# Test that Python state can load/save exactly like C++
from validate_restart_reader import IBPMRestart
from ibpm import Grid, State  # (your ported modules)

ref = IBPMRestart('examples/ibpm00000.bin')
grid = Grid(**ref.data['grid_params'])
state = State(grid, ref.data['numPoints'])

# Load reference data into Python state
state.omega = ref.get_omega()[0]
q_x, q_y = ref.get_flux_xy()
state.q.from_xy(q_x[0], q_y[0])
state.f = ref.get_forces()
state.time = ref.data['time']
state.timestep = ref.data['timestep']

# Save and compare
state.save('test.bin')
ref2 = IBPMRestart('test.bin')

# Should match
assert np.allclose(state.omega, ref2.get_omega()[0])
assert np.allclose(q_x, ref2.get_flux_xy()[0])
print("✓ Data structures validation passed")
```

### 3. After Porting Solvers
```python
# Test that Python solver produces same output as C++
ref0 = IBPMRestart('examples/ibpm00000.bin')
ref1 = IBPMRestart('examples/ibpm00001.bin')

# Create Python solver and load initial state
solver = NonlinearIBSolver(grid, model, dt=0.02)
state = State.from_reference(ref0)

# Advance one step
solver.advance(state)

# Compare to reference
ref1_data = ref1.get_reference_state()
assert np.allclose(state.omega, ref1_data['omega'], rtol=1e-8)
print("✓ Single timestep validation passed")
```

## Binary File Format

### Quick Reference

Each `.bin` file contains (in order):

```
[Header]
nx (int)           # Grid cells in x
ny (int)
ngrid (int)
dx (double)        # Grid spacing
x0 (double)        # Domain origin x
y0 (double)        # Domain origin y
numPoints (int)    # Number of boundary points

[Data]
Flux q             # (ngrid * numFluxes) doubles
Omega ω            # (ngrid * (nx-1) * (ny-1)) doubles (interior only)
Forces f           # (numPoints * 2) doubles (fx, fy pairs)

[Time]
timestep (int)
time (double)
```

where `numFluxes = (nx+1)*ny + nx*(ny+1)` (X-fluxes + Y-fluxes)

**Key Points**:
- Flux storage: X-fluxes first (indices 0 to (nx+1)*ny-1), then Y-fluxes
- Omega storage: interior points only (i=1..nx-1, j=1..ny-1); boundaries are zero
- Forces storage: flat array of (fx, fy) pairs for each boundary point
- All floating-point numbers are IEEE 754 double precision (8 bytes)
- No xShift/yShift saved (see State.cc comments)

### Flux Layout Details

**X-fluxes** (nx+1)×ny:
- Row-major storage: index = i*ny + j
- Located at cell edges in x-direction

**Y-fluxes** nx×(ny+1):
- Row-major storage: index = (nx+1)*ny + i*(ny+1) + j
- Located at cell edges in y-direction

For accessing flux components:
```python
q_x, q_y = restart.get_flux_xy()
# q_x has shape (ngrid, nx+1, ny)
# q_y has shape (ngrid, nx, ny+1)
```

## Reference Data Details

### Grid Configuration
- **Domain**: 4.0 × 4.0 (length × length)
- **Cell spacing**: dx=0.02 on 200×200 grid
- **Origin**: (-2.0, -2.0) [centered at (0,0)]
- **Single level**: ngrid=1 (no multigrid)

### Simulation Configuration
- **Timestep**: dt=0.02
- **Reynolds number**: Re=100
- **Integration**: RK3 (3-step Runge-Kutta)
- **Model**: Nonlinear Navier-Stokes
- **Geometry**: Cylinder, diameter=1.0, 160 boundary points
- **Duration**: 200 timesteps (t=0 to t=4.0)

### Available Timesteps
- `ibpm00000.bin`: t=0.0, step=0
- `ibpm00001.bin`: t=0.02, step=1
- `ibpm00002.bin`: t=0.04, step=2
- ...
- `ibpm00200.bin`: t=4.0, step=200

## API Reference

### IBPMRestart Class

```python
from validate_restart_reader import IBPMRestart

# Load a file
restart = IBPMRestart('examples/ibpm00000.bin')

# Access raw data
print(restart.data['nx'])
print(restart.data['time'])

# Extract flux components
q_x, q_y = restart.get_flux_xy()
# q_x: (ngrid, nx+1, ny)
# q_y: (ngrid, nx, ny+1)

# Extract vorticity (with zero boundaries)
omega = restart.get_omega()
# shape: (ngrid, nx, ny)

# Extract forces
forces = restart.get_forces()
# shape: (numPoints, 2)

# Print summary
restart.info()
```

### compare_restarts Function

```python
from validate_restart_reader import compare_restarts

result = compare_restarts(file1, file2, tol=1e-10)
# tol: relative tolerance for floating-point comparison
# Returns dict with 'match' (bool) and 'differences' (dict)
```

### ValidationHarness Class

```python
from test_validation_harness import ValidationHarness

harness = ValidationHarness('examples')

# Print data summary
harness.print_summary()

# Run validation tests
harness.test_grid_parameters()
harness.test_data_validity()
harness.test_monotonic_time()

# Get reference data for a specific file
ref_data = harness.get_reference_state('cylinder_test00000.bin')
print(ref_data.keys())  # omega, q_x, q_y, forces, timestep, time, grid_params
```

## Troubleshooting

### "No reference files found"
- Check that `examples/ibpm*.bin` files exist
- Run from project root directory, or adjust path: `ValidationHarness('../examples')`

### "ModuleNotFoundError: No module named 'numpy'"
```bash
pip3 install numpy scipy
```

### Binary file contains NaN
- Indicates C++ solver diverged or had numerical issues
- Check simulation parameters (dt, Re, grid resolution)
- If using custom geometry, verify it's properly specified

### Validation fails after porting code
- Start with simple test: load reference → save → reload → compare
- Verify indexing: C++ uses 1-indexed interior loops, Python should use 0-indexed
- Check array shapes and slicing carefully
- Use `rtol=1e-10` for double-precision comparison

## References

- `VALIDATION_SETUP.md` — Complete binary format specification
- `PORTING_GUIDE.md` — Step-by-step porting instructions
- `src/State.cc` — Implementation of binary I/O (ground truth)
- `src/Flux.h` — Flux storage and indexing

