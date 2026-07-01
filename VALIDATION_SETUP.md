# IBPM Python Port — Validation Workflow

## Overview

This document describes the validation setup for porting the IBPM C++ solver to Python. The workflow is built on **binary restart file comparison**: we run the C++ solver as ground truth and use Python readers to validate each ported component against it.

## Generated Reference Data

### Location
```
examples/
├── cylinder_test00000.bin    # Initial state (t=0, zero)
├── cylinder_test00001.bin    # After 1 timestep
├── ... 
├── cylinder_test00010.bin    # After 10 timesteps
└── cylinder_test.force       # Forces output
```

### How to Generate
```bash
cd examples
/path/to/ibpm-master/build/ibpm \
  -name cylinder_test \
  -outdir . \
  -geom cylinder.geom \
  -nx 64 \
  -ny 64 \
  -nsteps 10 \
  -dt 0.01 \
  -Re 100 \
  -tecplot 0 \
  -restart 1 \
  -force 1
```

This creates 11 restart files (timesteps 0–10) and force data.

## Binary Format Specification

All binary files use **little-endian** IEEE 754 doubles (8 bytes) and 4-byte integers.

### Restart File Structure
Each `.bin` file contains (in order):

| Component | Type | Count | Description |
|-----------|------|-------|-------------|
| nx | int | 1 | Grid points in x |
| ny | int | 1 | Grid points in y |
| ngrid | int | 1 | Number of grid levels (1 for single grid) |
| dx | double | 1 | Grid spacing |
| x0 | double | 1 | Left edge x-coordinate |
| y0 | double | 1 | Bottom edge y-coordinate |
| numPoints | int | 1 | Number of boundary points |
| **Flux q** | double | ngrid × numFluxes | All flux values |
| **Omega (vorticity)** | double | ngrid × (nx-1) × (ny-1) | Interior vorticity only |
| **Forces f** | double | numPoints × 2 | Boundary forces (fx, fy pairs) |
| timestep | int | 1 | Current timestep |
| time | double | 1 | Current simulation time |

where `numFluxes = (nx+1)×ny + nx×(ny+1)` (X-fluxes followed by Y-fluxes).

### Flux Data Layout

The Flux array is stored as a flat 1D vector for each grid level:

**X-fluxes** (indices 0 to (nx+1)×ny−1):
- Dimensions: (nx+1) × ny
- Layout: row-major, `index = i*ny + j` for point (i,j)
- Covers nx+1 x-edges × ny y-cells

**Y-fluxes** (indices (nx+1)×ny to end):
- Dimensions: nx × (ny+1)
- Layout: row-major, `index = (nx+1)*ny + i*(ny+1) + j` for point (i,j)
- Covers nx x-cells × ny+1 y-edges

### Omega (Vorticity) Data Layout

For each grid level, stores only **interior points** (no boundary):
- Iteration order: `for i in 1..nx-1: for j in 1..ny-1`
- Total elements per level: (nx−1) × (ny−1)
- Boundary values (i=0, i=nx−1, j=0, j=ny−1) are implicitly zero

### Forces Data Layout

Flat array of N boundary points with x and y components:
- Order: `[fx₁, fy₁, fx₂, fy₂, ..., fₙ, fyₙ]`
- Dimension: numPoints × 2

## Python Reader

See `validate_restart_reader.py` for the complete reader implementation.

### Basic Usage

```python
from validate_restart_reader import IBPMRestart

# Load a restart file
restart = IBPMRestart('cylinder_test00000.bin')

# Access data
print(f"Timestep: {restart.data['timestep']}")
print(f"Time: {restart.data['time']:.6f}")

# Get flux components
q_x, q_y = restart.get_flux_xy()  # Returns (ngrid, nx+1, ny) and (ngrid, nx, ny+1)

# Get vorticity (full array with zero boundary)
omega = restart.get_omega()  # Returns (ngrid, nx, ny)

# Get forces
forces = restart.get_forces()  # Returns (numPoints, 2)

# Print summary
restart.info()
```

### Comparing Two Restarts

```python
from validate_restart_reader import compare_restarts

result = compare_restarts('cylinder_test00000.bin', 'cylinder_test00001.bin', tol=1e-10)
print(f"Match: {result['match']}")
if result['differences']:
    for key, val in result['differences'].items():
        print(f"  {key}: {val}")
```

## Validation Checklist

Each ported module should be validated against these reference files:

- [ ] **Grid module**: Grid parameters match (nx, ny, ngrid, dx, offsets)
- [ ] **State module**: Can load/save restart files and recover all fields exactly
- [ ] **Flux module**: Flux reconstruction from state matches reference
- [ ] **Regularizer module**: Force spreading and velocity interpolation match reference
- [ ] **Elliptic solver**: Poisson/Helmholtz solutions match reference
- [ ] **Projection solver**: Constraint solve produces correct forces
- [ ] **IB solver**: Full timestep advances produce expected state changes
- [ ] **Full solver**: 10 timesteps produce forces/energy matching reference

## Testing Strategy

### Stage 1: Data Loading (Low Risk)
Start with the Python reader itself. Verify it:
1. Reads all file formats without errors
2. Produces identical arrays when re-saved and re-read
3. Comparison function correctly identifies differences

### Stage 2: Core Data Structures (Low Risk)
Port Grid, State, Scalar, Flux:
1. Load reference data into Python objects
2. Verify all indexing/slicing matches C++ semantics
3. Re-export to binary, compare byte-for-byte

### Stage 3: Linear Algebra (Medium Risk)
Port Regularizer, Elliptic solver, Projection solver:
1. Load reference flux/omega
2. Apply regularization, check results against reference
3. Solve Poisson equations, verify solutions

### Stage 4: Time Integration (High Risk)
Port IBSolver:
1. Load reference state at t=0
2. Advance one timestep, compare to reference at t=dt
3. Iteratively advance 10 steps, compare forces

## Example Validation Test

```python
import numpy as np
from validate_restart_reader import IBPMRestart
from ibpm_python import Grid, State, Regularizer  # (not yet implemented)

# Load reference
ref0 = IBPMRestart('cylinder_test00000.bin')
ref1 = IBPMRestart('cylinder_test00001.bin')

# Create Python grid/state
grid = Grid(nx=64, ny=64, ngrid=1, length=4.0, xOffset=-2.0, yOffset=-2.0)
state = State(grid, numPoints=160)

# Load reference data into Python state
state.omega = ref0.get_omega()[0]  # Level 0
q_x, q_y = ref0.get_flux_xy()
state.q = construct_flux_from_xy(q_x[0], q_y[0])  # (not implemented yet)
state.f = ref0.get_forces()

# Test regularizer
regularizer = Regularizer(grid, geometry)  # (not implemented)
flux_from_forces = regularizer.to_flux(state.f)
# Compare to reference...

# Test solver
solver = IBSolver(grid, model, dt=0.01)
solver.advance(state)
# Compare to reference at next timestep
assert np.allclose(state.omega, ref1.get_omega()[0], rtol=1e-10)
```

## Notes

- **Floating-point tolerance**: Use `rtol=1e-10` for comparing double-precision data
- **NaN handling**: If any restart file contains NaN (solver divergence), that is a hard error
- **Grid compatibility**: Reader validates nx, ny, ngrid, dx, x0, y0 on load
- **Offsets**: xShift, yShift are **not** saved in the binary format (see State.cc comments)

## Files Included

1. `validate_restart_reader.py` — Binary reader and comparison functions
2. `VALIDATION_SETUP.md` — This document
3. `examples/cylinder_test*.bin` — Reference restart files (generated)
4. `examples/cylinder_test.force` — Reference force data (generated)

