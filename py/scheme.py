# scheme.py
#
# Python port of src/Scheme.h
#
# Storage for time-integration scheme coefficients (Explicit Euler, Adams-
# Bashforth, and two 3rd-order Runge-Kutta variants). Each scheme is a table
# of per-substep coefficients (an, bn, cn).
#
# NOTE(port): Scheme.h is header-only and, unusually for this codebase, lives
# in the *global* namespace (not `ibpm`). That has no Python analogue; it is
# simply a module-level class here.
#
# ---------------------------------------------------------------------------
# JUDGMENT CALLS / PORT NOTES:
#
#   1. SchemeType enum.  C++ declares `enum SchemeType { EULER, AB2, RK3, RK3b }`
#      inside the class, referenced as `Scheme::EULER` etc. Ported as an
#      IntEnum `SchemeType` (module level, so it can also be used as a type
#      hint), with the four values additionally exposed as class attributes on
#      Scheme (Scheme.EULER, ...) to mirror the C++ `Scheme::EULER` spelling.
#
#   2. _coeff storage.  C++ uses Array::Array2<double> `_coeff`, allocated as
#      (nsteps, 3). Ported as a numpy array of shape (nsteps, 3). `nsteps()`
#      returns `_coeff.Nx()`, i.e. the number of rows (first dimension).
#
#   3. Unrecognized scheme.  The C++ prints an error and calls exit(1). Ported
#      as raising ValueError, the idiomatic Python equivalent of aborting on
#      bad input (rather than terminating the interpreter).
#
#   4. JAX-readiness.  Coefficients are held in a small numpy array and only
#      read element-wise; nothing here needs changing for a later JAX port.
# ---------------------------------------------------------------------------

from __future__ import annotations

from enum import IntEnum

import numpy as np


class SchemeType(IntEnum):
    EULER = 0
    AB2 = 1
    RK3 = 2
    RK3b = 3


class Scheme:
    """Time-integration scheme: a table of (an, bn, cn) coefficients, one row
    per substep."""

    # NOTE(port): expose enum values as class attributes so callers can write
    # Scheme.EULER etc., mirroring the C++ `Scheme::EULER`.
    SchemeType = SchemeType
    EULER = SchemeType.EULER
    AB2 = SchemeType.AB2
    RK3 = SchemeType.RK3
    RK3b = SchemeType.RK3b

    def __init__(self, scheme: SchemeType) -> None:
        if scheme == SchemeType.EULER:
            self._coeff = np.zeros((1, 3), dtype=np.float64)
            self._coeff[0, 0] = 1.0        # an
            self._coeff[0, 1] = 0.0        # bn
            self._coeff[0, 2] = 1.0        # cn
            self._name = "Explicit Euler"
        elif scheme == SchemeType.AB2:
            self._coeff = np.zeros((1, 3), dtype=np.float64)
            self._coeff[0, 0] = 3.0 / 2.0  # an
            self._coeff[0, 1] = -1.0 / 2.0  # bn
            self._coeff[0, 2] = 1.0        # cn
            self._name = "Adams Bashforth"
        elif scheme == SchemeType.RK3:
            self._coeff = np.zeros((3, 3), dtype=np.float64)
            self._coeff[0, 0] = 8.0 / 15.0  # an
            self._coeff[1, 0] = 5.0 / 12.0
            self._coeff[2, 0] = 3.0 / 4.0
            self._coeff[0, 1] = 0.0         # bn
            self._coeff[1, 1] = -17.0 / 60.0
            self._coeff[2, 1] = -5.0 / 12.0
            self._coeff[0, 2] = 8.0 / 15.0  # cn
            self._coeff[1, 2] = 2.0 / 3.0
            self._coeff[2, 2] = 1.0
            self._name = "3rd-order Runge Kutta (3-step)"
        elif scheme == SchemeType.RK3b:
            self._coeff = np.zeros((4, 3), dtype=np.float64)
            self._coeff[0, 0] = 8.0 / 17.0  # an
            self._coeff[1, 0] = 17.0 / 60.0
            self._coeff[2, 0] = 5.0 / 12.0
            self._coeff[3, 0] = 3.0 / 4.0
            self._coeff[0, 1] = 0.0         # bn
            self._coeff[1, 1] = -15.0 / 68.0
            self._coeff[2, 1] = -17.0 / 60.0
            self._coeff[3, 1] = -5.0 / 12.0
            self._coeff[0, 2] = 8.0 / 17.0  # cn
            self._coeff[1, 2] = 8.0 / 15.0
            self._coeff[2, 2] = 2.0 / 3.0
            self._coeff[3, 2] = 1.0
            self._name = "3rd-order Runge Kutta (4-step)"
        else:
            # NOTE(port): C++ prints an error and exit(1); ValueError is the
            # Python equivalent (port note 3).
            raise ValueError(f"ERROR: unrecognized solver: {scheme}")

    def an(self, i: int) -> float:
        return float(self._coeff[i, 0])

    def bn(self, i: int) -> float:
        return float(self._coeff[i, 1])

    def cn(self, i: int) -> float:
        return float(self._coeff[i, 2])

    def nsteps(self) -> int:
        return int(self._coeff.shape[0])

    def name(self) -> str:
        return self._name
