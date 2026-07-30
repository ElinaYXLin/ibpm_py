"""
test5_regrow.py

Follow-up test #5 from ../7-stripe_investigation/README.md's Proposals
section: is the upstream noise a static artifact that got locked in
early and just sits there, or is it continuously re-seeded by the
steady boundary force every step? Takes the converged baseline
snapshot, zeroes the upstream region of omega by hand, saves that as a
new initial-condition file (State.save(), the same restart format the
solver's own `-ic` flag loads via State.load() -- confirmed compatible
since ibpm.py calls model.refreshState(x) right after loading an IC,
which recomputes the flux q from omega, so editing omega alone is
self-consistent), and restarts the simulation from that edited field.
If the noise regrows fast, it's continuously re-seeded (a method
change is needed to remove it for good); if it doesn't come back, it
was a one-time startup transient (a one-time cleanup filter would
suffice). One short, cheap new run (dx=0.02, 300x150 grid).

Usage:
  python3 test5_regrow.py run       # build the edited IC + launch the run
  python3 test5_regrow.py analyze   # analyze + plot (default)
Output: figures/test5_regrow.png, data/test5_regrow.csv
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

DX = 0.02
BUFFER_DX = 2.0
YLIM = (-0.5, 0.5)
BASELINE_RUN = c.KURT1 / "runs" / "dx0.020" / "steady_py_a00"
BASELINE_NSTEPS = 3000
IC_PATH = c.RUNS / "regrow_ic.bin"
REGROW_OUTDIR = c.RUNS / "regrow"
REGROW_NSTEPS = 600
REGROW_RESTART = 20
REGROW_STEPS = list(range(0, REGROW_NSTEPS + 1, REGROW_RESTART))

FINE_OUTDIR = c.RUNS / "regrow_fine"
FINE_NSTEPS = 40
FINE_RESTART = 2
FINE_STEPS = list(range(0, FINE_NSTEPS + 1, FINE_RESTART))


def build_edited_ic():
    if IC_PATH.exists():
        print(f"edited IC already exists at {IC_PATH}")
        return
    state = c.load_state(BASELINE_RUN, BASELINE_NSTEPS)
    X, Y = c.grid_xy(DX)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]
    mask = c.upstream_mask(X, x_le, BUFFER_DX, DX)
    ys = Y[0, :]
    row_sel = (ys >= YLIM[0]) & (ys <= YLIM[1])
    before = float(np.abs(state.omega._data[0][mask, :][:, row_sel]).max())
    state.omega._data[0][np.ix_(mask, row_sel)] = 0.0
    after = float(np.abs(state.omega._data[0][mask, :][:, row_sel]).max())
    print(f"zeroed upstream omega: peak before={before:.4f}, after={after:.4f}")
    ok = state.save(str(IC_PATH))
    print(f"saved edited IC to {IC_PATH}: {'OK' if ok else 'FAILED'}")


def _launch(outdir, nsteps, restart):
    import subprocess
    import time
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-u", str(c.PY_RUNNER),
        "-geom", str(c.BASE_GEOM_DX002), "-name", "flow", "-outdir", str(outdir),
        "-nx", "300", "-ny", "150", "-ngrid", "1",
        "-length", str(c.DOMAIN["length"]), "-xoffset", str(c.DOMAIN["xoffset"]),
        "-yoffset", str(c.DOMAIN["yoffset"]), "-alpha", "0.0", "-Re", str(c.RE),
        "-dt", "0.01", "-nsteps", str(nsteps), "-tecplot", "0",
        "-restart", str(restart), "-force", "1",
        "-ic", str(IC_PATH), "-resettime", "1",
    ]
    t0 = time.time()
    with open(outdir / "run_log.txt", "w") as logf:
        proc = subprocess.run(cmd, cwd=c.REPO, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    ok = proc.returncode == 0 and c.is_done(outdir, nsteps)
    print(f"  {'OK' if ok else 'FAILED'} in {elapsed:.1f}s")


def do_run():
    build_edited_ic()
    if c.is_done(REGROW_OUTDIR, REGROW_NSTEPS):
        print("regrow run: already done")
    else:
        print("regrow run: launching (600 steps from edited IC, every 20)...", flush=True)
        _launch(REGROW_OUTDIR, REGROW_NSTEPS, REGROW_RESTART)

    if c.is_done(FINE_OUTDIR, FINE_NSTEPS):
        print("fine regrow run: already done")
    else:
        print("fine regrow run: launching (40 steps from edited IC, every 2)...", flush=True)
        _launch(FINE_OUTDIR, FINE_NSTEPS, FINE_RESTART)


def main():
    X, Y = c.grid_xy(DX)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]

    converged_om = c.load_omega(BASELINE_RUN, BASELINE_NSTEPS)
    converged_m = c.upstream_scalar_metrics(X, Y, converged_om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
    converged_ens = converged_m["enstrophy"]

    if not c.is_done(REGROW_OUTDIR, REGROW_NSTEPS):
        print("regrow run not found -- run `python3 test5_regrow.py run` first")
        return

    rows = []
    for step in REGROW_STEPS:
        om = c.load_omega(REGROW_OUTDIR, step)
        m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
        rows.append(dict(step=step, enstrophy=m["enstrophy"], peak=m["peak"]))
        print(f"step={step}: enstrophy={m['enstrophy']:.6f} "
              f"({m['enstrophy']/converged_ens*100:.1f}% of pre-edit converged value)")

    fine_rows = []
    if c.is_done(FINE_OUTDIR, FINE_NSTEPS):
        for step in FINE_STEPS:
            om = c.load_omega(FINE_OUTDIR, step)
            m = c.upstream_scalar_metrics(X, Y, om, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
            fine_rows.append(dict(step=step, enstrophy=m["enstrophy"], peak=m["peak"]))
        print(f"\nfine regrow: {len(fine_rows)} snapshots loaded (0-{FINE_NSTEPS}, every {FINE_RESTART})")
    else:
        print("\nfine regrow run not found -- run `python3 test5_regrow.py run` first")

    with open(c.DATA / "test5_regrow.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["step", "enstrophy", "peak"])
        w.writeheader(); w.writerows(rows)
    if fine_rows:
        with open(c.DATA / "test5_regrow_fine.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["step", "enstrophy", "peak"])
            w.writeheader(); w.writerows(fine_rows)
    print("wrote test5_regrow.csv" + (" and _fine.csv" if fine_rows else ""))

    steps = [r["step"] for r in rows]
    ens = [r["enstrophy"] for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axes[0]
    ax.plot(steps, ens, "o-", color="#c0392b", lw=1.4, ms=6, alpha=0.5, label="coarse (every 20 steps)")
    if fine_rows:
        ax.plot([r["step"] for r in fine_rows], [r["enstrophy"] for r in fine_rows], "o-",
                 color="#c0392b", lw=1.8, ms=5, label="fine (every 2 steps)")
    ax.axhline(converged_ens, color="gray", ls="--", lw=1.2, label="original converged value")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("timestep after restart (t = step x 0.01)"); ax.set_ylabel("upstream enstrophy")
    ax.set_title("Full range (0-600 steps)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)

    ax = axes[1]
    if fine_rows:
        ax.plot([r["step"] for r in fine_rows], [r["enstrophy"] for r in fine_rows], "o-",
                 color="#c0392b", lw=1.8, ms=6)
        ax.axhline(converged_ens, color="gray", ls="--", lw=1.2, label="original converged value")
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlim(0, FINE_NSTEPS)
    ax.set_xlabel("timestep after restart (t = step x 0.01)"); ax.set_ylabel("upstream enstrophy")
    ax.set_title(f"Zoomed to 0-{FINE_NSTEPS} steps (fine cadence)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)

    fig.suptitle("Test 5: does the upstream noise regrow after being zeroed once?\n"
                 "NACA0012, dx=0.02, alpha=0, steady, Re=1000", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = c.FIGS / "test5_regrow.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")

    def step_to_pct(step_list, vals, converged, pct):
        thresh = pct * converged
        for s, v in zip(step_list, vals):
            if v >= thresh:
                return s
        return None
    print(f"\n[coarse] regrowth reaches 50% of converged value by step "
          f"{step_to_pct(steps, ens, converged_ens, 0.5)}")
    print(f"[coarse] regrowth reaches 90% of converged value by step "
          f"{step_to_pct(steps, ens, converged_ens, 0.9)}")
    if fine_rows:
        fsteps = [r["step"] for r in fine_rows]
        fens = [r["enstrophy"] for r in fine_rows]
        print(f"[fine] regrowth reaches 50% of converged value by step "
              f"{step_to_pct(fsteps, fens, converged_ens, 0.5)}")
        print(f"[fine] regrowth reaches 90% of converged value by step "
              f"{step_to_pct(fsteps, fens, converged_ens, 0.9)}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        do_run()
    else:
        main()
