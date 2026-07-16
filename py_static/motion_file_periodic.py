# motion_file_periodic.py
#
# Python port of src/MotionFilePeriodic.h
#
# Subclass of Motion, for periodic motions defined piecewise-linearly by
# samples in a file.

from __future__ import annotations

import math
import sys
from typing import List, NamedTuple

import numpy as np

from .motion import Motion
from .tangent_se2 import TangentSE2


class TTSE2(NamedTuple):
    """A timestamped sample of (x, y, theta, dx, dy, dtheta), one period of
    which is repeated with period `_period`.

    (C++: `struct TTSE2` in MotionFilePeriodic.h.)
    """

    t: float
    x: float
    y: float
    theta: float
    dx: float
    dy: float
    dtheta: float


class MotionFilePeriodic(Motion):
    """A periodic Motion defined by piecewise-linear interpolation between
    samples read from a file, repeating with the given period."""

    def __init__(self, filename: str, period: float) -> None:
        self._filename = filename
        self._period = period
        self._data: List[TTSE2] = []
        # NOTE(port): see the identical NOTE(port) in motion_file.py --
        # the file-format and error-message handling here mirrors
        # MotionFile's constructor (same `>>`-style whitespace
        # tokenization, same non-fatal-in-C++-but-stop-early-in-Python
        # treatment of malformed input).
        try:
            with open(filename) as f:
                tokens = f.read().split()
        except OSError:
            print(f"ERROR:: MotionFilePeriodic: file {filename} formatted incorrectly! (err1)", file=sys.stderr)
            return
        if not tokens:
            print(f"ERROR:: MotionFilePeriodic: file {filename} formatted incorrectly! (err1)", file=sys.stderr)
            return
        try:
            n = int(tokens[0])
        except ValueError:
            print(f"ERROR:: MotionFilePeriodic: file {filename} formatted incorrectly! (err1)", file=sys.stderr)
            return

        tlast = -1.0e5
        pos = 1
        for _ in range(n):
            try:
                t, x, y, theta, dx, dy, dtheta = (float(v) for v in tokens[pos:pos + 7])
            except ValueError:
                print(f"ERROR:: MotionFilePeriodic: file {filename} formatted incorrectly! (err2)", file=sys.stderr)
                break
            pos += 7
            if t < tlast:
                print("ERROR:: MotionFilePeriodic: time must increase monotonically!", file=sys.stderr)
            tlast = t
            self.addTTSE2(t, x, y, theta, dx, dy, dtheta)

    def addTTSE2(self, t: float, x: float, y: float, theta: float, dx: float, dy: float, dtheta: float) -> None:
        """Adds an element of TTSE2 to the list of samples."""
        self._data.append(TTSE2(t, x, y, theta, dx, dy, dtheta))

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation for the piecewise-linear, periodic data
        in `_filename`."""
        time = math.fmod(time, self._period)
        n = len(self._data)
        # NOTE(port): see the identical NOTE(port) in motion_file.py's
        # getTransformation regarding "last match wins" vectorized search.
        index = -1
        if n >= 2:
            ts = np.array([d.t for d in self._data])
            mask = (time >= ts[:-1]) & (time < ts[1:])
            matches = np.nonzero(mask)[0]
            if matches.size > 0:
                index = int(matches[-1])

        xo = 0.0
        yo = 0.0
        thetao = 0.0
        dxo = 0.0
        dyo = 0.0
        dthetao = 0.0
        if index != -1:
            d0 = self._data[index]
            d1 = self._data[index + 1]
            tdiff = d1.t - d0.t
            alpha = (time - d0.t) / tdiff
            xo = alpha * d1.x + (1 - alpha) * d0.x
            yo = alpha * d1.y + (1 - alpha) * d0.y
            thetao = alpha * d1.theta + (1 - alpha) * d0.theta
            dxo = alpha * d1.dx + (1 - alpha) * d0.dx
            dyo = alpha * d1.dy + (1 - alpha) * d0.dy
            dthetao = alpha * d1.dtheta + (1 - alpha) * d0.dtheta
        return TangentSE2(xo, yo, thetao, dxo, dyo, dthetao)

    def clone(self) -> "MotionFilePeriodic":
        # NOTE(port): as in MotionFile.clone(), only the constructor
        # arguments are forwarded, so the clone re-reads and re-parses the
        # file from scratch.
        return MotionFilePeriodic(self._filename, self._period)
