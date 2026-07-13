# run_vortall.py
#
# Runs the C++ reference and the Python port independently out to
# `TARGET_STEPS`, each continuing only from its OWN previously-written
# restart file.
#
# This script exists to fix a real bug in how this dataset was previously
# produced: the Python continuation was once launched with
# `-ic results/vortall/_run_data_cpp/vortall13900.bin` -- i.e. seeded from
# the C++ run's own output -- so the last 13900/14000 steps (99.3%) of
# what was reported as the "Python result" were actually C++-computed,
# and Python only independently integrated the final 100 steps. That made
# the reported Python/C++ agreement close to tautological (two runs
# sharing 99.3% identical history will of course agree at the end) and
# silently discarded the real, fully-independent Python trajectory that
# a correct multi-chunk run had already produced.
#
# The fix is structural, not just "don't do that": `_last_checkpoint`
# below only ever globs inside the implementation's own output
# directory, so it is not possible for this script to pick up the other
# implementation's restart file, even by accident. Each implementation's
# run is independently resumable and always continues from its own state.
#
# Usage: python3 SURF_test/vortall/run_vortall.py
# (run from the repository root; requires build/ibpm -- see build/Makefile)

from __future__ import annotations

import pathlib
import subprocess
import sys
from typing import List, Optional, Tuple

REPO = pathlib.Path(__file__).resolve().parents[3]
GEOM = REPO / "examples" / "cylinder.geom"
OUTDIR_CPP = REPO / "SURF_test" / "vortall" / "1-baseline" / "_run_data_cpp"
OUTDIR_PY = REPO / "SURF_test" / "vortall" / "1-baseline" / "_run_data"

# Grid/physics parameters pinned by VORTALL.mat's own shape (89351 x 151
# = 449 x 199 interior vorticity nodes -> nx=450, ny=200); see README.
GRID = dict(nx=450, ny=200, ngrid=1, length=9, xoffset=-1, yoffset=-2,
            Re=100, dt=0.02, restart=100, force=1, tecplot=0)

TARGET_STEPS = 14000  # t = 280, saturated periodic vortex shedding


def _grid_args() -> List[str]:
    return [
        "-nx", str(GRID["nx"]), "-ny", str(GRID["ny"]), "-ngrid", str(GRID["ngrid"]),
        "-length", str(GRID["length"]), "-xoffset", str(GRID["xoffset"]),
        "-yoffset", str(GRID["yoffset"]), "-Re", str(GRID["Re"]), "-dt", str(GRID["dt"]),
        "-restart", str(GRID["restart"]), "-force", str(GRID["force"]),
        "-tecplot", str(GRID["tecplot"]),
    ]


def _last_checkpoint(outdir: pathlib.Path, name: str) -> Tuple[int, Optional[pathlib.Path]]:
    """Highest-numbered restart file already written by THIS
    implementation's own run in `outdir`, or (0, None) if none exists.
    Deliberately never looks outside `outdir` -- see module docstring."""
    bins = sorted(outdir.glob(f"{name}?????.bin"))
    if not bins:
        return 0, None
    last = bins[-1]
    return int(last.stem[len(name):]), last


def run_to(cmd_prefix: List[str], outdir: pathlib.Path, name: str,
           target_steps: int, log_path: pathlib.Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    done_step, ic_path = _last_checkpoint(outdir, name)
    if done_step >= target_steps:
        print(f"[{name}] already at step {done_step} >= target {target_steps}; nothing to do")
        return
    nsteps = target_steps - done_step
    cmd = [*cmd_prefix, "-geom", str(GEOM), "-name", name, "-outdir", str(outdir),
           *_grid_args(), "-nsteps", str(nsteps)]
    if ic_path is not None:
        cmd += ["-ic", str(ic_path)]
    print(f"[{name}] continuing from step {done_step} -> {target_steps} ({nsteps} steps)")
    print("  " + " ".join(cmd))
    with open(log_path, "a") as logf:
        logf.write(f"\n=== {' '.join(cmd)} ===\n")
        logf.flush()
        proc = subprocess.run(cmd, cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"[{name}] run failed (exit {proc.returncode}); see {log_path}")


def main() -> None:
    cpp_bin = REPO / "build" / "ibpm"
    if not cpp_bin.exists():
        print(f"ERROR: {cpp_bin} not found. Build it first: cd build && make", file=sys.stderr)
        sys.exit(1)

    run_to([str(cpp_bin)], OUTDIR_CPP, "vortall", TARGET_STEPS, OUTDIR_CPP / "run_log.txt")

    run_ibpm_runner = REPO / "SURF_test" / "cost" / "1-multi-core" / "run_ibpm_case.py"
    run_to([sys.executable, "-u", str(run_ibpm_runner)], OUTDIR_PY, "vortall", TARGET_STEPS,
           OUTDIR_PY / "run_log.txt")

    print("done")


if __name__ == "__main__":
    main()
