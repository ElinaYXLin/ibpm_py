# _fftw_native.py
#
# ctypes wrapper around _fftw_dst_shim.c -- gives elliptic_solver_2d.py a way
# to call the real FFTW3 C library, with the real FFTW_EXHAUSTIVE planning
# flag, exactly as src/EllipticSolver2d.cc does, instead of scipy.fft's
# pocketfft backend (which produces the same numbers but never performs
# FFTW's own timed search over its list of candidate sine-transform
# algorithms).
#
# Why this exists: an earlier version of this port dropped the FFTW plan
# machinery entirely (see the port-notes header this file's caller,
# elliptic_solver_2d.py, used to carry) on the grounds that scipy.fft is
# "functional" -- true for correctness, but not faithful to what the
# original C++ actually spends its startup time doing. This module restores
# that behavior for authenticity: same library, same flag, same per-instance
# plan lifecycle (built once at solver-construction time, reused for every
# solve for that solver's lifetime, exactly like C++'s _FFTWPlan member).
#
# Trade-off (intentional, per the mentor's request): FFTW_EXHAUSTIVE times
# every candidate algorithm it knows for the transform's exact size before
# picking one, which is slow -- for a short run this can dominate the whole
# run's wall-clock time (the "list of fast-sine-transform algorithms" is
# genuinely searched, not just picked heuristically as scipy does). This
# module deliberately does NOT fall back silently to scipy on failure: if
# the shim can't be built/loaded, it raises with a clear explanation, so a
# missing C compiler or FFTW3 dev headers doesn't quietly turn into a loss
# of this authenticity.

from __future__ import annotations

import ctypes
import pathlib
import subprocess
import sys
import sysconfig

_HERE = pathlib.Path(__file__).resolve().parent
_REPO = _HERE.parent
_SRC = _HERE / "_fftw_dst_shim.c"
_BUILD_DIR = _REPO / "build"
_DYLIB = _BUILD_DIR / ("_fftw_dst_shim" + sysconfig.get_config_var("SHLIB_SUFFIX") or ".so")

_INCLUDE_CANDIDATES = ["/usr/local/include", "/opt/homebrew/include", "/usr/include"]
_LIB_CANDIDATES = ["/usr/local/lib", "/opt/homebrew/lib", "/usr/lib"]


class FFTWUnavailableError(RuntimeError):
    """Raised when the native FFTW3 shim can't be built or loaded. Deliberately
    NOT caught anywhere to silently fall back to scipy -- see module docstring."""


def _build_shim() -> pathlib.Path:
    _BUILD_DIR.mkdir(parents=True, exist_ok=True)
    include_dir = next((d for d in _INCLUDE_CANDIDATES if (pathlib.Path(d) / "fftw3.h").exists()), None)
    lib_dir = next((d for d in _LIB_CANDIDATES
                     if list(pathlib.Path(d).glob("libfftw3.*"))), None)
    if include_dir is None or lib_dir is None:
        raise FFTWUnavailableError(
            "Could not find fftw3.h / libfftw3 on this machine (checked "
            f"{_INCLUDE_CANDIDATES} and {_LIB_CANDIDATES}). This port now calls "
            "the real FFTW3 library (via a small C shim) for its sine transform, "
            "for authenticity with src/EllipticSolver2d.cc's FFTW_EXHAUSTIVE "
            "planning -- install FFTW3 (the same dependency build/ibpm already "
            "needs; see README's 'Installation' section) to use it."
        )
    cc = "clang" if sys.platform == "darwin" else "cc"
    cmd = [cc, "-shared", "-O2", "-fPIC" if sys.platform != "darwin" else "-dynamiclib",
           "-I", include_dir, str(_SRC), "-L", lib_dir, "-lfftw3", "-o", str(_DYLIB)]
    # -dynamiclib already implies shared-object output on macOS; -shared is
    # redundant there but harmless to pass alongside it. On Linux, -shared
    # -fPIC is the standard pair; drop the (macOS-only) -dynamiclib flag.
    if sys.platform != "darwin":
        cmd = [c for c in cmd if c != "-dynamiclib"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not _DYLIB.exists():
        raise FFTWUnavailableError(
            f"Failed to compile {_SRC.name} into a shared library:\n"
            f"  command: {' '.join(cmd)}\n  stdout: {proc.stdout}\n  stderr: {proc.stderr}"
        )
    return _DYLIB


_lib = None


def _get_lib():
    global _lib
    if _lib is not None:
        return _lib
    if not _DYLIB.exists():
        _build_shim()
    try:
        lib = ctypes.CDLL(str(_DYLIB))
    except OSError:
        # stale/incompatible cached build (e.g. left over from a different
        # machine) -- rebuild once before giving up
        _build_shim()
        lib = ctypes.CDLL(str(_DYLIB))
    lib.dst2d_create.restype = ctypes.c_void_p
    lib.dst2d_create.argtypes = [ctypes.c_int, ctypes.c_int]
    lib.dst2d_execute.restype = None
    lib.dst2d_execute.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]
    lib.dst2d_destroy.restype = None
    lib.dst2d_destroy.argtypes = [ctypes.c_void_p]
    _lib = lib
    return lib


class NativeDST2D:
    """One FFTW_EXHAUSTIVE-planned 2D DST-I (FFTW_RODFT00) transform, sized
    n0 x n1. Mirrors EllipticSolver2d's _FFTWPlan/_fft members: the
    (expensive) plan search happens once, in __init__, and every call to
    execute() reuses it -- exactly the same lifecycle as the C++ object that
    owns it."""

    def __init__(self, n0: int, n1: int) -> None:
        import numpy as np  # local import: keep ctypes wiring numpy-free above

        self._np = np
        self._n0, self._n1 = n0, n1
        self._lib = _get_lib()
        self._handle = self._lib.dst2d_create(n0, n1)
        if not self._handle:
            raise FFTWUnavailableError("fftw_plan_r2r_2d returned NULL (out of memory?)")
        # reusable contiguous scratch buffer for the ctypes call
        self._scratch = np.empty((n0, n1), dtype=np.float64, order="C")

    def execute(self, u):
        """Return the DST-I of 2D array u (shape (n0, n1)), unnormalized --
        same convention as FFTW_RODFT00 / scipy.fft.dstn(type=1)."""
        np = self._np
        self._scratch[...] = u
        ptr = self._scratch.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        self._lib.dst2d_execute(self._handle, ptr)
        return self._scratch.copy()

    def close(self) -> None:
        if getattr(self, "_handle", None):
            self._lib.dst2d_destroy(self._handle)
            self._handle = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
