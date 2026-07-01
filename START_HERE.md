# IBPM Python Port — Start Here

## What You Have

A complete **validation framework** for porting the IBPM C++ solver to Python. The framework uses **binary restart file comparison** — the C++ solver is ground truth, and every ported Python module is validated against its output.

## In 5 Minutes

```bash
# 1. Check setup
python3 test_validation_harness.py summary

# 2. Try the reader
python3 -c "
from validate_restart_reader import IBPMRestart
r = IBPMRestart('examples/ibpm00000.bin')
r.info()
print('✓ Validation framework is working')
"

# 3. You're ready to port!
```

## What's Inside

| Item | Purpose |
|------|---------|
| **validate_restart_reader.py** | Read IBPM binary files (with examples) |
| **test_validation_harness.py** | Pre-loaded reference data for testing |
| **Documentation** | 5 files explaining format, approach, and examples |
| **Reference data** | 11 timesteps of C++ output (200×200 grid, t=0–4.0) |

## How It Works

```
┌─────────────────────────────────────────────────────────┐
│  C++ Solver (Ground Truth)                              │
│  ↓ Produces restart files ↓                             │
│  examples/ibpm00000.bin ... ibpm00200.bin               │
└──────────────────┬──────────────────────────────────────┘
                   │
                   │ Load with Python
                   ↓
┌─────────────────────────────────────────────────────────┐
│  validate_restart_reader.py                             │
│  - IBPMRestart class                                    │
│  - compare_restarts() function                          │
│  - get_flux_xy(), get_omega(), get_forces()            │
└──────────────────┬──────────────────────────────────────┘
                   │
                   │ Use as test oracle
                   ↓
┌─────────────────────────────────────────────────────────┐
│  Your Python Port                                       │
│  - Port Grid → validate                                 │
│  - Port State → validate                                │
│  - Port IBSolver → validate                             │
│  Each component tested against C++ output               │
└─────────────────────────────────────────────────────────┘
```

## Quick Examples

### Load a reference file
```python
from validate_restart_reader import IBPMRestart

r = IBPMRestart('examples/ibpm00000.bin')
q_x, q_y = r.get_flux_xy()      # Flux components
omega = r.get_omega()           # Vorticity (with zero boundary)
forces = r.get_forces()         # Boundary forces

print(f"Grid: {r.data['nx']}×{r.data['ny']}")
print(f"Time: {r.data['time']:.6f}, Step: {r.data['timestep']}")
```

### Validate your ported code
```python
from validate_restart_reader import compare_restarts
from ibpm import IBSolver  # Your ported code

# Run solver
solver = IBSolver(grid, model, dt=0.02)
solver.advance(state)
state.save('my_output.bin')

# Compare to reference
match, diffs = compare_restarts('my_output.bin', 'examples/ibpm00001.bin')
if match:
    print("✓ Your code produces identical output!")
else:
    print("Differences:", diffs)
```

## The 5-Stage Plan

1. **Data Structures** (1 week) — Port Grid, State, Flux
   - Test: Load reference → Save → Re-load → Match

2. **Linear Algebra** (1 week) — Port Regularizer, Elliptic solver, Projection solver
   - Test: Apply operations, compare to reference

3. **Geometry** (3-4 days) — Port Geometry, Motion, Output
   - Test: Load geometry, verify point coordinates

4. **Model** (1 week) — Port NavierStokes, BaseFlow
   - Test: Verify B and C operators produce correct results

5. **Solver** (1-2 weeks) — Port IBSolver, Scheme, time integration
   - Test: Single step, then 10+ steps, compare forces

**Total: 4-5 weeks** (can parallelize some stages)

## Key Files to Read

1. **Quick Reference**: `QUICK_REFERENCE.md` (2 min read)
   - API summary, common operations, debugging checklist

2. **Validation README**: `VALIDATION_README.md` (10 min read)
   - Quick start, API reference, troubleshooting

3. **Binary Format**: `VALIDATION_SETUP.md` (20 min read)
   - Complete format specification, layout details

4. **Porting Guide**: `PORTING_GUIDE.md` (30 min read)
   - 5-stage plan with code examples for each stage

## The Reference Data

- **Grid**: 200×200 cells, dx=0.02 (from x=-2 to x=2, y=-2 to y=2)
- **Simulation**: Re=100, RK3 scheme, 200 timesteps, dt=0.02
- **Boundary**: Cylinder (160 points, diameter=1.0)
- **Files**: `examples/ibpm00000.bin` (t=0) to `ibpm00200.bin` (t=4.0)
- **Quality**: Valid double-precision data (no NaN or divergence)

## Getting Started

### Option 1: Quick Start (5 minutes)
```bash
# Verify everything works
python3 test_validation_harness.py summary

# You're ready — pick a file to port from PORTING_GUIDE.md
```

### Option 2: Understand First (1 hour)
1. Read `QUICK_REFERENCE.md`
2. Read `VALIDATION_README.md`
3. Try the examples above
4. Read `PORTING_GUIDE.md`

### Option 3: Full Deep Dive (3 hours)
1. Read all 4 documentation files
2. Study the binary format in detail
3. Review C++ source code for the component you'll port first
4. Set up your Python project structure
5. Start with Stage 1 (data structures)

## What's NOT Included

- Python implementation of any solver components
- CMake or build configuration
- Poetry/pip requirements.txt (minimal: just numpy, scipy)
- GitHub Actions CI/CD

These are **your** responsibility, but the validation framework will help you get them right.

## Validation Mindset

Before you write any solver code:

1. **Understand the C++ code**
   - Read the corresponding `.h` and `.cc` files
   - Understand data flow (inputs → outputs)
   - Identify any assumptions (index conventions, boundary handling)

2. **Load reference data**
   - Use `IBPMRestart` to load C++ output
   - Extract the specific data your component uses
   - Understand the shapes and ranges

3. **Plan the test**
   - Define inputs (state, geometry, parameters)
   - Define expected outputs
   - Calculate tolerance (usually rtol=1e-10 for doubles)

4. **Port the code**
   - Write as Python code (not C++ translated line-by-line)
   - Use NumPy/SciPy idioms
   - Avoid premature optimization

5. **Validate**
   - Load reference input
   - Run your Python code
   - Compare to reference output
   - If no match: debug, compare intermediate values, repeat

## Troubleshooting

**"ModuleNotFoundError: No module named 'numpy'"**
```bash
pip3 install numpy scipy
```

**"No reference files found"**
```bash
# Make sure you're in the project root
ls -la examples/ibpm*.bin   # Should show 11 files
```

**Validation fails**
1. Check that reference file is valid: `python3 validate_restart_reader.py examples/ibpm00000.bin`
2. Verify grid parameters (nx, ny, dx)
3. Check for NaN or infinity in your output
4. Review PORTING_GUIDE.md for debugging strategies

## Next Steps

1. **Read** `QUICK_REFERENCE.md` (2 minutes)
2. **Verify** `python3 test_validation_harness.py summary` (1 minute)
3. **Choose** a component to port from `PORTING_GUIDE.md`
4. **Port** that component (refer to documentation)
5. **Validate** using the harness
6. **Commit** and move to next component

## Files You'll Use Most Often

```
# Reading reference data
python3 validate_restart_reader.py examples/ibpm00100.bin

# Testing your code (after porting)
python3 -c "
from validate_restart_reader import IBPMRestart, compare_restarts
from ibpm import YourComponent  # Your ported code

ref = IBPMRestart('examples/ibpm00000.bin')
# ... run your code ...
result = compare_restarts('my_output.bin', 'examples/ibpm00001.bin')
assert result['match']
"

# Reading documentation
cat QUICK_REFERENCE.md          # 2-minute lookup
cat VALIDATION_README.md        # Full API reference
cat PORTING_GUIDE.md            # Detailed examples
```

## Philosophy

- **Validation over speculation**: Test against C++ output, not against other Python code
- **Correctness first**: Get it right, optimize later
- **Incremental progress**: Port one component, validate it, move on
- **Preserve behavior**: Don't "improve" algorithms during porting; match C++ exactly
- **Document decisions**: Record why you translated something a particular way

## You're Ready!

Everything is set up. Reference data is generated. Readers are written. Documentation is complete.

**Pick a component from PORTING_GUIDE.md and start porting.**

The validation framework will tell you when you've got it right.

---

Questions? Check:
- Stuck on API? → `QUICK_REFERENCE.md`
- Need examples? → `PORTING_GUIDE.md`
- Want format details? → `VALIDATION_SETUP.md`
- Validation failing? → `VALIDATION_README.md` (Troubleshooting section)

Good luck! 🚀
