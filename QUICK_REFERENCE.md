# IBPM Python Port — Quick Reference Card

## Validation Workflow TL;DR

```bash
# Check reference data is valid
python3 test_validation_harness.py summary

# Read a file
python3 validate_restart_reader.py examples/ibpm00000.bin

# Compare two files
python3 -c "
from validate_restart_reader import compare_restarts
r = compare_restarts('a.bin', 'b.bin', tol=1e-10)
print('Match:', r['match'])
"
```

## Python API

```python
from validate_restart_reader import IBPMRestart, compare_restarts

# Load a file
r = IBPMRestart('examples/ibpm00000.bin')

# Get data
q_x, q_y = r.get_flux_xy()      # (ngrid, nx+1, ny), (ngrid, nx, ny+1)
omega = r.get_omega()           # (ngrid, nx, ny)
forces = r.get_forces()         # (numPoints, 2)

# Check properties
print(f"Grid: {r.data['nx']}x{r.data['ny']}")
print(f"Time: {r.data['time']:.6f}")
print(f"Step: {r.data['timestep']}")

# Compare files
match, diffs = compare_restarts('file1.bin', 'file2.bin')
```

## Binary File Format

| Component | Type | Count | Bytes |
|-----------|------|-------|-------|
| nx, ny, ngrid | int × 3 | 1 | 12 |
| dx, x0, y0 | double × 3 | 1 | 24 |
| numPoints | int | 1 | 4 |
| Flux q | double | ngrid×numFluxes | — |
| Omega ω | double | ngrid×(nx-1)×(ny-1) | — |
| Forces f | double | numPoints×2 | — |
| timestep, time | int, double | 1 | 12 |

where `numFluxes = (nx+1)×ny + nx×(ny+1)`

## Data Access Patterns

```python
# Extract at timestep
q_x, q_y = r.get_flux_xy()          # Separate X/Y components
omega = r.get_omega()               # Full array with zero boundary
forces = r.get_forces()             # (numPoints, 2)

# Access metadata
grid_params = {
    'nx': r.data['nx'],
    'ny': r.data['ny'],
    'dx': r.data['dx'],
    'x0': r.data['x0'],
    'y0': r.data['y0'],
}

# Access time info
t = r.data['time']
step = r.data['timestep']
```

## Porting Stages

| Stage | Files | Duration | Risk |
|-------|-------|----------|------|
| 1 | Grid, State, Scalar, Flux | 1 week | Low |
| 2 | Regularizer, Elliptic, Projection | 1 week | Medium |
| 3 | Geometry, Motion, Output | 3-4 days | Low |
| 4 | NavierStokes, BaseFlow | 1 week | Medium |
| 5 | IBSolver, Scheme, Integration | 1-2 weeks | High |

## Validation Checklist

- [ ] Reference data loads without error
- [ ] All timesteps have valid (non-NaN) data
- [ ] Grid parameters are consistent across files
- [ ] Time increases monotonically
- [ ] Data ranges are reasonable (check extreme values)

## Common Operations

### Load reference and create Python objects
```python
from validate_restart_reader import IBPMRestart

ref = IBPMRestart('examples/ibpm00000.bin')
grid = Grid(nx=ref.data['nx'], ny=ref.data['ny'], ...)
state = State(grid, numPoints=ref.data['numPoints'])
state.omega = ref.get_omega()[0]
```

### Test a ported component
```python
ref0 = IBPMRestart('examples/ibpm00000.bin')
ref1 = IBPMRestart('examples/ibpm00001.bin')

# Your ported code
solver = IBSolver(...)
state = State.from_reference(ref0)
solver.advance(state)

# Validate
assert np.allclose(state.omega, ref1.get_omega()[0], rtol=1e-8)
print("✓ Test passed")
```

### Compare outputs
```python
from validate_restart_reader import compare_restarts

result = compare_restarts(
    'my_output.bin',
    'examples/ibpm00001.bin',
    tol=1e-10
)

if not result['match']:
    print("Differences:")
    for key, val in result['differences'].items():
        print(f"  {key}: {val}")
else:
    print("✓ Output matches reference")
```

## Key Numbers (Reference Data)

- **Grid**: 200×200 cells, dx=0.02
- **Domain**: 4.0×4.0, origin at (-2,-2)
- **Timestep**: dt=0.02
- **Reynolds**: Re=100
- **Boundary points**: 160 (cylinder)
- **Available timesteps**: 0-200 (t=0 to t=4.0)
- **Grid levels**: 1 (single domain)

## Files to Know

| File | Purpose |
|------|---------|
| `validate_restart_reader.py` | Load/compare binary files |
| `test_validation_harness.py` | Pre-loaded reference data + tests |
| `VALIDATION_SETUP.md` | Binary format spec |
| `PORTING_GUIDE.md` | Detailed porting examples |
| `examples/ibpm*.bin` | Reference restart files |

## Floating-Point Tolerance

Use these guidelines:
- **Same algorithm**: `rtol=1e-14` (machine epsilon)
- **Same numerical method**: `rtol=1e-10` (default)
- **Different algorithm**: `rtol=1e-6` (order of magnitude different)
- **Always**: Check absolute difference too if near zero: `atol=1e-15`

```python
np.allclose(a, b, rtol=1e-10, atol=0)
```

## Debugging Failed Validation

1. **Check file validity**
   ```bash
   python3 validate_restart_reader.py examples/ibpm00000.bin
   # Should print data without errors
   ```

2. **Identify mismatch type**
   ```python
   from validate_restart_reader import compare_restarts
   r = compare_restarts('a.bin', 'b.bin')
   print(r['differences'].keys())
   # One of: grid_*, q_x, q_y, omega, forces, time, timestep
   ```

3. **Inspect specific component**
   ```python
   r1 = IBPMRestart('a.bin')
   r2 = IBPMRestart('b.bin')
   
   import numpy as np
   q_x1, _ = r1.get_flux_xy()
   q_x2, _ = r2.get_flux_xy()
   
   diff = np.abs(q_x1 - q_x2)
   print(f"Max diff: {diff.max()}")
   print(f"Location: {np.unravel_index(np.argmax(diff), diff.shape)}")
   ```

4. **Check indices**
   - C++ interior: 1..nx-1, 1..ny-1
   - Python interior: 0..nx-2, 0..ny-2 (if 0-indexed)
   - Boundaries are implicit zero

## Performance Notes

- **Early stage**: Correctness first, performance later
- **Array operations**: Vectorize with NumPy (avoid Python loops)
- **FFT**: Use scipy.fft.dst, check normalization carefully
- **Linear algebra**: Use scipy.linalg.cholesky, scipy.sparse.linalg.cg
- **Profile later**: Use cProfile or line_profiler after porting

## File Locations

```
ibpm-master/
├── validate_restart_reader.py       # Reader (use this)
├── test_validation_harness.py       # Test framework (use this)
├── VALIDATION_*.md                  # Documentation
├── QUICK_REFERENCE.md               # This file
└── examples/
    ├── ibpm*.bin                    # Reference data (11 files)
    └── cylinder.geom                # Geometry
```

## Emergency Reference

Can't remember what method to call? Check here:

```python
# Loading
r = IBPMRestart('file.bin')

# Accessing flux
q_x, q_y = r.get_flux_xy()
# q_x shape: (ngrid, nx+1, ny)
# q_y shape: (ngrid, nx, ny+1)

# Accessing vorticity
omega = r.get_omega()
# shape: (ngrid, nx, ny)

# Accessing forces
forces = r.get_forces()
# shape: (numPoints, 2)

# Accessing metadata
t = r.data['time']
step = r.data['timestep']
nx, ny = r.data['nx'], r.data['ny']

# Comparing
match = compare_restarts(file1, file2)[0]
```

---

**For detailed documentation**, see:
- `VALIDATION_README.md` — Full API and examples
- `PORTING_GUIDE.md` — Step-by-step guide with test templates
- `VALIDATION_SETUP.md` — Complete binary format specification
