# Computational cost analysis: C++ vs. Python

Measures and compares runtime, CPU, and RAM cost of the C++ (`src/`,
`build/ibpm`) and Python (`py/ibpm.py`) implementations, across a range of
grid resolutions, under two different core-usage conditions:

- **[`1-multi-core/`](1-multi-core/)** -- each backend run under its
  environment's *default* threading behavior (whatever numpy's BLAS
  backend / FFTW decide to do on their own). This is where the original
  cost analysis lived; see its own `README.md` for the full writeup,
  including a documented benchmark-methodology bug (a `.cholesky` cache
  file silently making reruns measure cache-load time instead of a fresh
  factorization) found and fixed while re-verifying it.
- **[`2-single_core/`](2-single_core/)** -- the same benchmark, same
  methodology, but with Python forced onto a single CPU core's worth of
  parallelism (via `OPENBLAS_NUM_THREADS=1` and friends, set as
  environment variables at subprocess-launch time -- no code in `py/` was
  touched). Added because C++ was already found to be inherently
  single-core (no OpenMP, no threaded FFTW anywhere in `src/`), but
  Python's numpy is linked against a multi-threaded OpenBLAS build, so a
  fair single-core comparison required explicitly ruling that out rather
  than assuming it away. **Headline finding**: multi-core Python spiked to
  610-787% peak CPU (6-8 cores) at its busiest moment, vs. C++'s ~112-115%
  ceiling in either condition -- but wall-clock time barely changed (at
  most ~0.03x), showing that parallelism wasn't buying real speed at these
  problem sizes, just extra (wasted) CPU-seconds. See that directory's
  `README.md` for the full numbers.
