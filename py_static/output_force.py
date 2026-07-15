# output_force.py
#
# Python port of src/OutputForce.h / src/OutputForce.cc
#
# Output routine for writing a list of force coefficients.

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional, TextIO

from .output import Output

if TYPE_CHECKING:
    from .base_flow import BaseFlow
    from .state import State


class OutputForce(Output):
    """Output routine for writing a list of force coefficients."""

    def __init__(self, filename: str) -> None:
        self._filename: str = filename
        # NOTE(port): C++ `_fp` is an uninitialized `FILE*` data member
        # until `init()` is called (it is not set in the constructor body);
        # initialized here to None instead of leaving it undefined.
        self._fp: Optional[TextIO] = None

    def init(self) -> bool:
        """Open the file for writing. If a file with the same name is
        already present, it is overwritten. Returns true if successful."""
        try:
            self._fp = open(self._filename, "w")
        except OSError:
            return False
        return True

    def cleanup(self) -> bool:
        """Close the file. Returns true if successful."""
        status = True
        if self._fp is not None:
            self._fp.close()
        return status

    def doOutput(self, *args) -> bool:
        """Write data to the force file.

        NOTE(port): collapses the three C++ overloads
            bool doOutput(const double alpha, const double mag, const State& x);
            bool doOutput(const BaseFlow& q, const State& x);
            bool doOutput(const State& x);
        into one method dispatching on the number of positional arguments
        (3, 2, or 1 respectively), matching the overload-collapse
        convention used throughout this port.
        """
        if len(args) == 3:
            alpha, mag, x = args
            return self._doOutputForceCoeffs(alpha, mag, x)
        elif len(args) == 1:
            (x,) = args
            # If no other information is provided, assume zero angle of
            # attack, unity freestream velocity
            alpha = 0.0
            mag = 1.0
            return self._doOutputForceCoeffs(alpha, mag, x)
        elif len(args) == 2:
            q, x = args
            alpha = q.getAlpha()
            mag = q.getMag()
            return self._doOutputForceCoeffs(alpha, mag, x)
        else:
            raise TypeError("OutputForce.doOutput: unsupported arguments")

    def _doOutputForceCoeffs(self, alpha: float, mag: float, x: "State") -> bool:
        """Compute lift, drag from state (x), angle of attack (alpha), and
        freestream velocity (mag).

        NOTE(port): matches the C++ body of the 3-argument
        `doOutput(alpha, mag, x)` exactly, including that `mag` is accepted
        but never used in the computation below -- that is already true of
        the original C++ function body, not a translation artifact.
        """
        xF, yF = x.computeNetForce()
        drag = xF * math.cos(alpha) + yF * math.sin(alpha)
        lift = xF * -1.0 * math.sin(alpha) + yF * math.cos(alpha)

        # Convert forces to lift and drag coefficients:
        # If L_d is dimensional lift, then in the nondimensionalization of
        # the code (lengths by c, density by rho, velocity by U), we have
        #    L = L_d / (c rho U^2)
        # so
        #    C_L = L_d / (1/2 rho U^2)
        #        = 2 L
        drag *= 2.0
        lift *= 2.0

        if self._fp is None:
            return False
        self._fp.write("%5d %.5e %.5e %.5e\n" % (x.timestep, x.time, drag, lift))
        self._fp.flush()

        return True
