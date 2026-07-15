# Cross-validation tests: py/vector_operations.py vs. C++ reference

These tools numerically cross-validate the Python port of
`src/VectorOperations.cc` against the original C++ implementation, at full
`float64` precision, on the exact scenarios used by the pre-existing (but
non-automated / interactive) standalone C++ programs in `test/CheckCurl.cc`,
`test/CheckLaplacian.cc`, `test/CheckFluxToX.cc`, and `test/CheckAdjoint.cc`.

They are separate from (and do not replace) the C++ gtest suite in `test/`
(which does not build on modern toolchains -- gtest 1.6.0 requires
`<tr1/tuple>`, removed from libc++ in C++11+) and from
`test_validation_harness.py` at the repo root (which only checks the
checked-in `examples/cylinder_test*.bin` reference files for internal
consistency; it doesn't exercise this Python port at all).

## Why raw-value dumps instead of diffing printed output

The original `test/Check*.cc` programs print grids with C++'s default
`cout` precision (6 significant digits) via `Scalar::print()` /
`Flux::print()`, and some end in an interactive `cin.get()` loop meant for
manual eyeballing, not automated comparison. Each `cpp/dump_*.cc` here is a
minimal rewrite of the corresponding `Check*.cc` scenario that:

- removes the interactive loop (replaced with a fixed spot-check in
  `dump_fluxtox.cc`'s case -- see its header comment),
- prints every value at full precision (`setprecision(17)`), one value per
  line, under a `DUMP_<name>` marker, so a bitwise/near-bitwise comparison
  is trivial.

`compare_dumps.py` parses two such dumps and reports, section by section,
whether every value matches to a given relative tolerance (default `1e-9`),
with NaN-aware comparison (a NaN in one file only matches a NaN at the same
position in the other -- see "About NaNs" below).

## What's here

```
cpp/                    C++ dump drivers (not integrated into test/Makefile;
                         build standalone against libibpm.a, see each file's
                         header comment for the exact g++ invocation)
  dump_curl.cc             Curl (Flux->Scalar and Scalar->Flux), + adjoint check
  dump_laplacian.cc        Laplacian + Scalar.coarsify() (checkLCEqualsCLC)
  dump_fluxtox.cc          FluxToXVelocity / XVelocityToFlux (round-trip + spot-check)
  dump_adjoint.cc          FluxToYVelocity / YVelocityToFlux adjoint check

python/                 Python-side equivalents (py package, run with -m)
  dump_curl.py
  dump_laplacian.py
  dump_fluxtox.py
  dump_adjoint.py
  cross_product_adjoint.py   Standalone (no C++ build needed): validates
                              CrossProduct's documented adjoint identity
                              <a, q1 x q2> = <q1, q2 x a> at Ngrid=1,2,3
  inner_products.py          Test-only InnerProduct(Scalar,Scalar) /
                              InnerProduct(Flux,Flux), transcribed from
                              VectorOperations.cc. NOT part of the production
                              port (py/vector_operations.py doesn't define
                              InnerProduct -- nothing there needs it yet).

compare_dumps.py        Parses and diffs two dump files
```

## How to run

From the repo root, with `build/libibpm.a` already built (`cd build && make`)
and `fftw3` available:

```bash
# 1. Build a C++ dump driver
g++ -Wall -O2 -I src -c py/tests/cross_validation/cpp/dump_curl.cc -o /tmp/dump_curl.o
g++ /tmp/dump_curl.o -L build -libpm -lfftw3 -lm -o /tmp/dump_curl
/tmp/dump_curl > /tmp/dump_curl_cpp.out

# 2. Run the Python equivalent
python3 -m py.tests.cross_validation.python.dump_curl > /tmp/dump_curl_py.out

# 3. Compare
python3 py/tests/cross_validation/compare_dumps.py /tmp/dump_curl_cpp.out /tmp/dump_curl_py.out
```

Repeat for `dump_laplacian`, `dump_fluxtox`, `dump_adjoint`.

`cross_product_adjoint.py` needs no C++ build:

```bash
python3 -m py.tests.cross_validation.python.cross_product_adjoint
```

## Results as of this writing

All four dump-based comparisons and the standalone adjoint-identity check
passed at (near-)machine precision:

| Check | Result |
|---|---|
| `dump_curl` (Curl both directions + adjoint, 8 sections) | bitwise-exact match (`max_rel_diff = 0`) |
| `dump_laplacian` (`checkLCEqualsCLC`, 98 points x 2 levels) | identical pass/fail pattern (65/98 fail in both -- a genuine property of the discretization, not a bug) |
| `dump_fluxtox` (round-trip + 7-point spot-check, 3 levels) | bitwise-exact match |
| `dump_adjoint` (FluxToYVelocity adjoint, per-basis-vector, 3 levels) | max relative diff `2.1e-16` (float64 epsilon) |
| `cross_product_adjoint` (Ngrid = 1, 2, 3) | relative diff `< 1e-10` (assertion in script) |

## About NaNs

Several checks legitimately produce NaN at specific grid points, identically
in C++ and Python: for `lev >= 1`, the multigrid `InnerProduct(Scalar,
Scalar)` formula (see `VectorOperations.cc` lines 249-310, transcribed in
`inner_products.py`) deliberately excludes the "interior duplicate" square of
each coarse level (points that coincide with the next-finer grid, to avoid
double-counting). A basis vector `e` set at one of those interior points
therefore has `InnerProduct(e, e) == 0`, so any adjoint computation that
divides by that norm (`AdjYVelToFlux`, `AdjointOfScalarCurl`) produces `0/0
= NaN` there in both languages. `compare_dumps.py` treats this as a match
when the NaN positions agree exactly, which they do everywhere in these
tests.
