#!/usr/bin/env python3
"""Compare two raw-value dump files produced by a cpp/dump_*.cc driver and
its python/dump_*.py counterpart (see README.md).

Dump file format: a "DUMP_<name>" marker line, followed by one repr()'d /
setprecision(17) float per line, until the next marker or EOF.

NaN values are treated specially: a NaN in one file only matches a NaN at
the same position in the other file (this arises legitimately at multigrid
points whose basis-vector norm is exactly zero -- see the README).

Usage:
    python3 compare_dumps.py <cpp_output> <python_output> [--tol 1e-9]
"""

import argparse
import sys

import numpy as np


def parse(path):
    sections = {}
    cur = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("DUMP_"):
                cur = line
                sections[cur] = []
            elif cur is not None and line:
                sections[cur].append(float(line))
    return sections


def compare(cpp_path, py_path, tol):
    cpp = parse(cpp_path)
    py = parse(py_path)

    all_ok = True
    for key in cpp:
        if key not in py:
            print(f"{key}: MISSING IN PYTHON OUTPUT")
            all_ok = False
            continue
        c = np.array(cpp[key])
        p = np.array(py[key])
        if c.shape != p.shape:
            print(f"{key}: SHAPE MISMATCH cpp={c.shape} python={p.shape}")
            all_ok = False
            continue
        nan_match = np.array_equal(~np.isfinite(c), ~np.isfinite(p))
        fin = np.isfinite(c) & np.isfinite(p)
        rel = np.abs(c[fin] - p[fin]) / np.maximum(np.abs(c[fin]), 1e-14)
        worst = rel.max() if fin.any() else 0.0
        ok = nan_match and worst < tol
        status = "MATCH" if ok else "MISMATCH"
        print(f"{key}: n={len(c)} nan_match={nan_match} max_rel_diff={worst:.3e} -> {status}")
        all_ok = all_ok and ok

    extra = set(py) - set(cpp)
    if extra:
        print(f"Sections only in Python output (ignored): {sorted(extra)}")

    print()
    print("ALL SECTIONS MATCH:" if all_ok else "SOME SECTIONS MISMATCH:", all_ok)
    return all_ok


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpp_output")
    parser.add_argument("python_output")
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()
    ok = compare(args.cpp_output, args.python_output, args.tol)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
