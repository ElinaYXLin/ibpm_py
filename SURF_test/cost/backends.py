# backends.py
#
# Registry of "backends" (implementations of the IBPM solver) that the
# computational cost analysis in this directory can benchmark. Each
# backend just needs to supply a subprocess command line; run_benchmark.py
# is otherwise implementation-agnostic.
#
# This exists as its own module (rather than being inlined in
# run_benchmark.py) specifically so that adding a future JAX port is a
# matter of writing one `build_cmd` function and appending one `Backend`
# entry below -- nothing else in this directory needs to change. The
# `implemented=False` JAX stub is deliberately left in place (per current
# instructions: build the architecture now, wire up JAX later).

from __future__ import annotations

import sys
import pathlib
from dataclasses import dataclass
from typing import Callable, List, Optional

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
COST_DIR = pathlib.Path(__file__).resolve().parent

DOMAIN = dict(length=4.0, xoffset=-2.0, yoffset=-2.0)
RE = 100.0
SCHEME = "rk3"


def _common_grid_args(nx: int, ny: int, dt: float, nsteps: int) -> List[str]:
    return [
        "-nx", str(nx), "-ny", str(ny), "-ngrid", "1",
        "-length", str(DOMAIN["length"]), "-xoffset", str(DOMAIN["xoffset"]),
        "-yoffset", str(DOMAIN["yoffset"]),
        "-Re", str(RE), "-dt", str(dt), "-nsteps", str(nsteps), "-scheme", SCHEME,
        "-tecplot", "0", "-restart", "0", "-force", "1",
    ]


def cpp_cmd(geom: pathlib.Path, outdir: pathlib.Path, name: str,
            nx: int, ny: int, dt: float, nsteps: int) -> List[str]:
    binpath = REPO_ROOT / "build" / "ibpm"
    return [str(binpath), "-geom", str(geom), "-outdir", str(outdir), "-name", name,
            *_common_grid_args(nx, ny, dt, nsteps)]


def python_cmd(geom: pathlib.Path, outdir: pathlib.Path, name: str,
               nx: int, ny: int, dt: float, nsteps: int) -> List[str]:
    runner = COST_DIR / "run_ibpm_case.py"
    return [sys.executable, "-u", str(runner), "-geom", str(geom), "-outdir", str(outdir),
            "-name", name, *_common_grid_args(nx, ny, dt, nsteps)]


def jax_cmd(geom: pathlib.Path, outdir: pathlib.Path, name: str,
            nx: int, ny: int, dt: float, nsteps: int) -> List[str]:
    raise NotImplementedError(
        "JAX backend not implemented yet -- this stub exists only to reserve "
        "the architecture slot. When a JAX port exists, point this at its "
        "entry point (e.g. `python jax_ibpm.py -geom ... -nx ... `) with the "
        "same CLI surface as cpp_cmd/python_cmd above and flip `implemented=True` "
        "on the JAX Backend entry below -- nothing else in run_benchmark.py or "
        "gen_cost_report.py needs to change."
    )


@dataclass
class Backend:
    name: str            # short id, used in filenames/JSON keys, e.g. "cpp"
    label: str           # display label for tables/figures
    implemented: bool
    build_cmd: Optional[Callable[..., List[str]]] = None
    note: str = ""


BACKENDS: List[Backend] = [
    Backend("cpp", "C++ (src/, build/ibpm)", True, cpp_cmd),
    Backend("python", "Python (py/ibpm.py)", True, python_cmd),
    Backend(
        "jax", "JAX (planned)", False, jax_cmd,
        note="Architecture reserved for a future JAX port; not implemented "
             "or benchmarked in this pass (see module docstring).",
    ),
]


def get_backend(name: str) -> Backend:
    for b in BACKENDS:
        if b.name == name:
            return b
    raise KeyError(f"no such backend: {name!r}")


def make_circle_geom(dx: float, out_path: pathlib.Path, diameter: float = 1.0) -> int:
    """Write a `circle_n`-based cylinder .geom file whose boundary-point
    spacing is matched to grid spacing `dx` (spacing ~ dx avoids both the
    'leaky boundary' problem of too few points and the 'over-resolved
    boundary -> singular projection matrix' problem of too many -- see
    SURF_test/high_re/SD7003/README.md for a worked example of the latter failure
    mode on a non-circular body)."""
    import math
    n = max(int(round(math.pi * diameter / dx)), 8)
    out_path.write_text(
        f"# Cylinder, diameter {diameter}, {n} points (spacing ~ dx={dx:g})\n\n"
        f"body Cylinder\n  circle_n 0 0 {diameter / 2} {n}\nend\n"
    )
    return n
