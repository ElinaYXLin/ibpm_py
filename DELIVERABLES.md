# IBPM Python Port — Validation Workflow Setup

## Deliverables Summary

This package provides everything needed to validate a Python port of the IBPM C++ solver against reference data.

### 1. Python Validation Tools

**`validate_restart_reader.py`** (406 lines)
- `IBPMRestart` class: Read binary restart files
- Methods: `get_flux_xy()`, `get_omega()`, `get_forces()`, `info()`
- `compare_restarts()`: Compare two files with tolerance
- Reverse-engineered from `src/State.cc`
- Tested and working ✓

**`test_validation_harness.py`** (315 lines)
- `ValidationHarness` class: Pre-loaded reference data
- Methods: `test_grid_parameters()`, `test_data_validity()`, `test_monotonic_time()`
- `get_reference_state()`: Extract data for ported code validation
- Tested and working ✓

### 2. Reference Data

**`examples/ibpm00000.bin` through `ibpm00200.bin`** (11 files, 99 KB each)
- C++ solver output from complete simulation
- Grid: 200×200, dx=0.02, Re=100
- Timesteps: 0-200, dt=0.02, total time t=4.0
- Quality: Valid double-precision (no NaN)
- Verified with reader ✓

### 3. Documentation (6 files)

| File | Size | Purpose |
|------|------|---------|
| `START_HERE.md` | 5 KB | **Read this first** — Overview and quickstart |
| `QUICK_REFERENCE.md` | 8 KB | API summary, common operations, debugging |
| `VALIDATION_README.md` | 8.3 KB | Full quick-start and API reference |
| `VALIDATION_SETUP.md` | 6.9 KB | Binary format specification (complete) |
| `PORTING_GUIDE.md` | 12 KB | 5-stage porting plan with code examples |
| `VALIDATION_SUMMARY.txt` | 9.6 KB | Status report and next steps |

### 4. What It Enables

✓ **Load C++ restart files in Python**
```python
from validate_restart_reader import IBPMRestart
r = IBPMRestart('examples/ibpm00000.bin')
```

✓ **Compare outputs to reference**
```python
from validate_restart_reader import compare_restarts
match = compare_restarts('my_output.bin', 'examples/ibpm00001.bin')[0]
```

✓ **Validate each ported component**
- Load reference input
- Run Python code
- Compare to reference output
- Pass/fail with tolerance checking

✓ **5-stage porting plan**
1. Data structures (Grid, State, Flux)
2. Linear algebra (Regularizer, Elliptic solver)
3. Geometry (Geometry, Motion, Output)
4. Model (NavierStokes, BaseFlow)
5. Integration (IBSolver, Scheme, time stepping)

## Getting Started

### Step 1: Verify Setup (1 minute)
```bash
cd /Users/elina/Desktop/SURF2026/ibpm-master
python3 test_validation_harness.py summary
```

Expected: Grid info, 11 timesteps, no errors

### Step 2: Read Documentation (30 minutes)
1. `START_HERE.md` — Overview
2. `QUICK_REFERENCE.md` — API summary
3. Pick a component from `PORTING_GUIDE.md`

### Step 3: Try the Reader (2 minutes)
```python
from validate_restart_reader import IBPMRestart
r = IBPMRestart('examples/ibpm00000.bin')
r.info()
q_x, q_y = r.get_flux_xy()
print(f"Flux shape: {q_x.shape}")
```

Expected: q_x shape (1, 201, 200), all values valid

### Step 4: Start Porting
Pick a component, follow the template from `PORTING_GUIDE.md`, validate against reference.

## Binary File Format

Each `.bin` file contains:
- Grid info (nx, ny, ngrid, dx, x0, y0)
- Flux data (ngrid × [(nx+1)×ny + nx×(ny+1)] doubles)
- Vorticity (ngrid × (nx-1)×(ny-1) interior points only)
- Boundary forces (numPoints × 2 pairs)
- Time info (timestep, time)

See `VALIDATION_SETUP.md` for complete details.

## Reference Data Properties

- **Grid**: 200×200 cells, dx=0.02
- **Domain**: 4.0×4.0 ([-2, 2] × [-2, 2])
- **Simulation**: Re=100, RK3, nonlinear N-S
- **Boundary**: Cylinder, 160 points
- **Timesteps**: 11 available (0–200, steps 0, 1, 2, ..., 200)
- **Time range**: t=0 to t=4.0 (dt=0.02)
- **Quality**: All data valid (no NaN, proper ranges)

## Validation Strategy

For each ported component:

1. **Load reference input**
   ```python
   ref = IBPMRestart('examples/ibpm00000.bin')
   ```

2. **Create Python object**
   ```python
   from ibpm import YourClass
   obj = YourClass(ref.data['grid_params'])
   ```

3. **Run Python code**
   ```python
   output = obj.compute(ref.get_omega())
   ```

4. **Compare to reference**
   ```python
   from validate_restart_reader import compare_restarts
   assert compare_restarts('my_output.bin', 'reference.bin')[0]
   ```

See `PORTING_GUIDE.md` for detailed examples for each stage.

## What's NOT Included

- Python source code for solver components (that's your job!)
- Setup.py or pyproject.toml (minimal: numpy, scipy)
- CI/CD configuration
- Build system

## Testing Philosophy

- **Validation over speculation**: Test against C++ output
- **Correctness first**: Accuracy > speed initially
- **Incremental**: Port one component, validate, move on
- **Document decisions**: Record translation choices

## Files Checklist

- ✓ `validate_restart_reader.py` — Binary reader
- ✓ `test_validation_harness.py` — Test framework
- ✓ `examples/ibpm*.bin` — Reference data (11 files)
- ✓ `examples/cylinder.geom` — Geometry
- ✓ `START_HERE.md` — Entry point
- ✓ `QUICK_REFERENCE.md` — API summary
- ✓ `VALIDATION_README.md` — Full documentation
- ✓ `VALIDATION_SETUP.md` — Format specification
- ✓ `PORTING_GUIDE.md` — Implementation guide
- ✓ `VALIDATION_SUMMARY.txt` — Status report
- ✓ `DELIVERABLES.md` — This file

## Dependencies

- Python 3.6+
- numpy
- scipy

Install:
```bash
pip3 install numpy scipy
```

## Time Estimate

- Setup verification: 5 minutes
- Documentation review: 30 minutes
- Stage 1 (Data structures): 1 week
- Stage 2 (Linear algebra): 1 week
- Stage 3 (Geometry): 3-4 days
- Stage 4 (Model): 1 week
- Stage 5 (Integration): 1-2 weeks

**Total: 4-5 weeks for complete port**

## Success Criteria

✓ Data structures pass I/O validation (load/save exact match)
✓ Linear algebra produces reference results (rtol=1e-10)
✓ Full integration matches reference forces and energy
✓ 10 timesteps produce correct evolution
✓ All tests pass with default tolerance (rtol=1e-10, atol=0)

## Next Action

**Read `START_HERE.md` now.**

It will guide you through verification, understanding the framework, and choosing your first component to port.

---

**Contact**: Questions? Check the documentation files. Validation failing? See Troubleshooting in `VALIDATION_README.md`.
