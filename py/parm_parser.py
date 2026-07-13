# parm_parser.py
#
# Python port of src/ParmParser.h / src/ParmParser.cc
#
# Class for parsing command-line arguments.

from __future__ import annotations

import re
import sys
from typing import Callable, IO, List, Optional, TypeVar

_OPTION_WIDTH = 20

_T = TypeVar("_T")

# NOTE(port): mimic C++ `istringstream >> T` extraction: it succeeds as
# long as the token *begins* with a valid literal for T (any trailing,
# non-matching characters are simply left unconsumed rather than causing
# failure); it only fails if no valid prefix is present at all. Python's
# `int()`/`float()` instead require the *entire* string to match, so a
# leading-prefix regex scan is used below to reproduce the C++ "prefix
# succeeds" behavior for getInt/getDouble/getBool.
_INT_RE = re.compile(r"[+-]?\d+")
_FLOAT_RE = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")


def _extract_int(token: str) -> Optional[int]:
    """Mimic `istringstream >> int` on `token`; see module note above."""
    m = _INT_RE.match(token)
    if m is None:
        return None
    return int(m.group(0))


def _extract_float(token: str) -> Optional[float]:
    """Mimic `istringstream >> double` on `token`; see module note above."""
    m = _FLOAT_RE.match(token)
    if m is None:
        return None
    return float(m.group(0))


def _extract_bool(token: str) -> Optional[bool]:
    """Mimic `istringstream >> bool` (without `boolalpha`, the stream
    default) on `token`: C++ parses a bool the same way as an int (0 ->
    false, nonzero -> true); see module note above."""
    val = _extract_int(token)
    if val is None:
        return None
    return val != 0


def _extract_str(token: str) -> str:
    """Mimic `istringstream >> string` on `token`: always succeeds (the
    token itself, already whitespace-delimited)."""
    return token


def _format_default_double(val: float) -> str:
    """Format a double the way C++'s default `ostream <<` does (general
    format, 6 significant digits, no trailing zeros/decimal point) --
    used when building the usage-message default-value string."""
    return f"{val:.6g}"


class ParmParser:
    """Class for parsing command-line arguments."""

    def __init__(self, argc: int, argv: List[str]) -> None:
        """Constructor, taking command line arguments as input.

        NOTE(port): preserves the exact two-argument C++ signature
        `ParmParser(int argc, char* argv[])`, even though `argc` is
        redundant with `len(argv)` in Python; callers should pass e.g.
        `ParmParser(len(sys.argv), sys.argv)`.
        """
        self._argc: int = argc
        self._argv: List[str] = list(argv)
        self._used: List[bool] = [False] * argc

        # NOTE(port): the C++ constructor builds `_args` by concatenating
        # argv[1..argc-1] with a trailing space after each, then later
        # re-tokenizes `_args` via whitespace-delimited istream extraction
        # (`in >> s`) in getFlag/getParm below, rather than indexing
        # `_argv` directly. Reproduced faithfully (join then re-split via
        # `str.split()`, which -- like `istringstream >> string` -- splits
        # on runs of whitespace); this means an argument containing
        # embedded whitespace would be incorrectly split into multiple
        # tokens here, exactly as in the original C++ code.
        self._args: str = ""
        for i in range(1, argc):
            self._args += argv[i] + " "

        self._argOut: str = argv[0]

        # remove path from executable name, if specified on command line
        exec_name = argv[0]
        slash_position = exec_name.rfind("/")
        if slash_position >= 0:
            exec_name = exec_name[slash_position + 1:]

        self._usageMessage: str = (
            f"USAGE: {exec_name} [options]\n"
            "where [options] are as follows (defaults shown in [ ]):\n"
        )

    def _appendUsageMessageFlag(self, flag: str, description: str) -> None:
        """Append a usage-message line for a flag (no argument).

        NOTE(port): collapses one of C++'s two private
        `appendUsageMessage` overloads (flag+description) -- see
        `_appendUsageMessageParm` for the other.
        """
        self._usageMessage += "  " + ("-" + flag).ljust(_OPTION_WIDTH) + description + "\n"

    def _appendUsageMessageParm(self, parm: str, typ: str, description: str, defVal: str) -> None:
        """Append a usage-message line for a parm taking an argument.

        NOTE(port): the other of C++'s two `appendUsageMessage` overloads
        (parm+type+description+defaultVal).
        """
        label = "-" + parm + " " + typ
        self._usageMessage += "  " + label.ljust(_OPTION_WIDTH) + description + " [" + defVal + "]\n"

    def getFlag(self, flag: str, description: str) -> bool:
        """Search the parameter list for the given flag, that does not
        take an argument. Returns true if the flag is present."""
        self._appendUsageMessageFlag(flag, description)
        target = "-" + flag

        tokens = self._args.split()
        for i, s in enumerate(tokens, start=1):
            if s == target:
                self._used[i] = True
                self._argOut += " " + target
                return True
        return False

    def _getParm(
        self,
        parm: str,
        defaultVal: _T,
        extractor: Callable[[str], Optional[_T]],
        formatter: Callable[[_T], str] = str,
    ) -> _T:
        """Generic function to search for the given entry parm and return
        its argument or a default value.

        NOTE(port): C++ implements this as a class template
        `template<class T> T ParmParser::getParm(...)`; Python has no
        templates, so this single method is parameterized by an
        `extractor` callable (one of _extract_int/_extract_float/
        _extract_bool/_extract_str) instead, dispatched from
        getInt/getDouble/getString/getBool below (matching the
        overload-collapse convention used throughout this port).

        NOTE(port): C++ echoes the found value into `_argOut` with
        `_argOut << " " << parm << " " << val`, i.e. using `ostream::operator<<`
        for the value's own type T. The `formatter` callable reproduces that
        type-specific formatting (default `str`, matching `<< int`/`<< string`;
        overridden for double -> `%.6g` and bool -> "1"/"0" by the callers), so
        the saved/echoed parameter list matches C++ byte-for-byte rather than
        using Python's fuller `str(float)` / "True"/"False" renderings.
        """
        target = "-" + parm
        tokens = self._args.split()
        n = len(tokens)
        for i, s in enumerate(tokens, start=1):
            if s != target:
                continue
            # Try to get the next argument in the list
            val = extractor(tokens[i]) if i < n else None
            if val is not None:
                self._used[i] = True
                self._used[i + 1] = True
                self._argOut += f" {target} {formatter(val)}"
                return val
            sys.stderr.write(f"Warning: cannot parse argument {i}: {target}\n")
            return defaultVal
        return defaultVal

    def getInt(self, parm: str, description: str, defaultVal: int) -> int:
        """Search the parameter list for the given entry parm and a
        single integer argument, returning defaultVal if not specified.
        If argument is invalid, print a warning message and return
        defaultVal."""
        self._appendUsageMessageParm(parm, "<int>", description, str(defaultVal))
        return self._getParm(parm, defaultVal, _extract_int)

    def getDouble(self, parm: str, description: str, defaultVal: float) -> float:
        """Search the parameter list for description, and return the
        corresponding double value, or defaultVal if not specified."""
        self._appendUsageMessageParm(parm, "<real>", description, _format_default_double(defaultVal))
        return self._getParm(parm, defaultVal, _extract_float, _format_default_double)

    def getString(self, parm: str, description: str, defaultVal: str) -> str:
        """Search the parameter list for description, and return the
        corresponding string value, or defaultVal if not specified."""
        self._appendUsageMessageParm(parm, "<string>", description, defaultVal)
        return self._getParm(parm, defaultVal, _extract_str)

    def getBool(self, parm: str, description: str, defaultVal: bool) -> bool:
        """Search the parameter list for description, and return the
        corresponding boolean value, or defaultVal if not specified."""
        self._appendUsageMessageParm(parm, "[0 or 1]", description, "1" if defaultVal else "0")
        return self._getParm(parm, defaultVal, _extract_bool, lambda v: "1" if v else "0")

    def inputIsValid(self) -> bool:
        """Check if any input parameters were invalid, or unused.

        Returns true if everything is valid and all parameters used. If
        any were not used, prints them to standard error and returns
        false.
        """
        allused = True
        for i in range(1, self._argc):
            allused = allused and self._used[i]

        if allused:
            return True
        sys.stderr.write("Warning: the following parameters were not used:\n  ")
        for i in range(1, self._argc):
            if not self._used[i]:
                sys.stderr.write(self._argv[i] + " ")
        sys.stderr.write("\n")
        return False

    def printUsage(self, out: IO[str]) -> None:
        """Print a usage message."""
        out.write(self._usageMessage)

    def getParameters(self) -> str:
        """Return argument list in a string."""
        return self._argOut

    def saveParameters(self, fname: str) -> None:
        """Save argument list to a file."""
        with open(fname, "w") as out:
            out.write(self._argOut + "\n")
