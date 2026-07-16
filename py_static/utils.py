# utils.py
#
# Python port of src/utils.h / src/utils.cc
#
# Utilities for working with strings.

from __future__ import annotations


def EatWhitespace(s: str) -> str:
    """Remove whitespace at the beginning of the string s.

    NOTE(port): C++ takes `string&` and mutates the caller's string in
    place (erasing leading characters one at a time while `isspace(s[0])`).
    Python strings are immutable, so this cannot mutate in place; it
    returns the trimmed string instead. Callers must be updated to use the
    return value (`s = EatWhitespace(s)`) rather than relying on a
    side-effecting `s.strip()`-in-place. Behaviorally equivalent to
    `s.lstrip()` using the same whitespace definition as C's `isspace`
    (`str.lstrip()` with no arguments strips the same set of ASCII
    whitespace characters for ordinary input).
    """
    return s.lstrip()


def MakeLowercase(s: str) -> str:
    """Convert the string s to lower case.

    NOTE(port): same in-place-vs-immutable caveat as EatWhitespace above;
    returns the lowercased string instead of mutating in place.
    """
    return s.lower()


def AddSlashToPath(s: str) -> str:
    """Add a slash to the end of the string s, if not already present.

    NOTE(port): same in-place-vs-immutable caveat as EatWhitespace above.
    """
    if len(s) > 0 and s[-1] != "/":
        return s + "/"
    return s
