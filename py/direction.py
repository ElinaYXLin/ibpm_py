# direction.py
#
# Python port of src/Direction.h
#
# Type for specifying directions of vectors (X or Y).

from __future__ import annotations

from enum import IntEnum


class Direction(IntEnum):
    """Direction of a vector component (X or Y), or XY for "both"."""

    X = 0
    Y = 1
    XY = 2

# NOTE(port): C++ defines postfix/prefix `operator++` on Direction, used
# nowhere in Grid/Field/Scalar/Flux (grep shows no call sites in the files
# ported so far). Not ported to keep this a faithful-but-minimal port;
# add `Direction(int(dir) + 1)` at the call site if/when a dependent file
# needs increment semantics.
