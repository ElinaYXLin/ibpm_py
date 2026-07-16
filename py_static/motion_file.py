# motion_file.py
#
# Python port of src/MotionFile.h
#
# Subclass of Motion, for a (non-periodic) motion defined piecewise-
# linearly by samples in a file.

from __future__ import annotations

import sys
from typing import List, NamedTuple

import numpy as np

from .motion import Motion
from .tangent_se2 import TangentSE2


class TxTSE2(NamedTuple):
    """A timestamped sample of (x, y, theta, dx, dy, dtheta).

    (C++: `struct TxTSE2` in MotionFile.h.)
    """

    t: float
    x: float
    y: float
    theta: float
    dx: float
    dy: float
    dtheta: float


class MotionFile(Motion):
    """A Motion defined by piecewise-linear interpolation between samples
    read from a file: each line/whitespace-separated record is
    `t x y theta dx dy dtheta`, preceded by a count `n`."""

    def __init__(self, filename: str) -> None:
        self._filename = filename
        self._data: List[TxTSE2] = []
        # NOTE(port): C++ reads with `ifstream >> n` then `>> t >> x >> y
        # >> theta >> dx >> dy >> dtheta` repeated n times -- i.e. tokens
        # are whitespace-delimited, not necessarily one record per line.
        # On a missing file or a parse failure, C++ prints an error to
        # cerr via `in.fail()` checks but does *not* stop construction --
        # it continues attempting `n` reads of a stream already in a fail
        # state, which yields undefined/garbage values for t,x,y,... (a
        # real bug in the original). There is no faithful way to
        # reproduce reading garbage values from a failed C++ stream in
        # Python, so this instead stops after printing the same error
        # message, leaving `_data` empty -- a deliberate, flagged
        # deviation from the (buggy) exact C++ behavior on malformed
        # input, while preserving the original's behavior on well-formed
        # input and its error *message*.
        try:
            with open(filename) as f:
                tokens = f.read().split()
        except OSError:
            print(f"ERROR:: MotionFile: file {filename} formatted incorrectly! (err1)", file=sys.stderr)
            return
        if not tokens:
            print(f"ERROR:: MotionFile: file {filename} formatted incorrectly! (err1)", file=sys.stderr)
            return
        try:
            n = int(tokens[0])
        except ValueError:
            print(f"ERROR:: MotionFile: file {filename} formatted incorrectly! (err1)", file=sys.stderr)
            return

        tlast = -1.0e5
        pos = 1
        for _ in range(n):
            try:
                t, x, y, theta, dx, dy, dtheta = (float(v) for v in tokens[pos:pos + 7])
            except ValueError:
                print(f"ERROR:: MotionFile: file {filename} formatted incorrectly! (err2)", file=sys.stderr)
                break
            pos += 7
            if t < tlast:
                print("ERROR:: MotionFile: time must increase monotonically!", file=sys.stderr)
            tlast = t
            self.addTxTSE2(t, x, y, theta, dx, dy, dtheta)

    def addTxTSE2(self, t: float, x: float, y: float, theta: float, dx: float, dy: float, dtheta: float) -> None:
        """Adds an element of TxTSE2 to the list of samples."""
        self._data.append(TxTSE2(t, x, y, theta, dx, dy, dtheta))

    def getTransformation(self, time: float) -> TangentSE2:
        """Returns transformation for the piecewise-linear data in
        `_filename`."""
        n = len(self._data)
        # NOTE(port): the C++ loop `for (i=0; i<n-1; ++i) if (data[i].t <=
        # time < data[i+1].t) index = i;` scans every bracketing interval
        # and keeps the *last* match (later loop iterations overwrite
        # `index`). Since t is expected to increase monotonically, at most
        # one interval matches in practice, but the vectorized numpy
        # version below reproduces "last match wins" exactly (via [-1])
        # rather than assuming uniqueness, to stay faithful even if the
        # monotonicity precondition is violated.
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

    def clone(self) -> "MotionFile":
        # NOTE(port): C++ `clone()` only forwards `_filename`, so the
        # clone re-reads and re-parses the file from scratch (rather than
        # copying `_data`) -- reproduced exactly, even though copying
        # `_data` directly would be cheaper.
        return MotionFile(self._filename)
