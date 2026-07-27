"""
test8_mechanism.py

Test 8: mechanism isolation. Two cheap probes into WHICH stage of the
solver produces the upstream noise.

Probe A (zero new runs) -- overlay the upstream vorticity pattern
against the curl of the spread boundary force: `model.B(f, omega_out)`
computes exactly `omega_out = Curl(regularizer.toFlux(f))` (see
py_static/navier_stokes_model.py's B()), i.e. the vorticity-space
footprint the discrete delta function alone would produce from the
converged boundary force f. Loaded directly from an existing steady
snapshot's saved state (State.f is part of the restart format). If this
spread-force pattern correlates strongly with the real upstream omega,
the regularization/projection step is the source, not the elliptic
(streamfunction) solve.

Probe B (one new run, nsteps=1) -- every run in this repo starts from a
zero-vorticity (uniform flow) initial condition unless a restart file is
explicitly loaded; running exactly one timestep therefore isolates
whatever the FIRST application of the projection step produces, before
any advection has had a chance to accumulate or transport anything. If
the upstream pattern is already present after step 1, it is produced
directly by the projection/regularization, not built up over many
steps of nonlinear advection.

Usage: python3 test8_mechanism.py [run|analyze]  (default: run then analyze)
Output: figures/test8_mechanism.png, data/test8_mechanism.csv
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c

YLIM = (-0.5, 0.5)
BUFFER_DX = 2.0
DX = 0.02
NX, NY = 300, 150
NSTEPS_BASELINE = 3000
BASELINE_RUN = c.KURT1 / "runs" / "dx0.020" / "steady_py_a00"
ONE_STEP_OUTDIR = c.RUNS / "one_step"


def do_run():
    if c.is_done(ONE_STEP_OUTDIR, 1):
        print("one-step run: already done")
        return
    print("one-step run: launching (nsteps=1, from zero IC)...", flush=True)
    ok, elapsed = c.run_case("py", c.BASE_GEOM_DX002, ONE_STEP_OUTDIR, NX, NY, 0.01, 1,
                               alpha=0.0, restart=1)
    print(f"  {'OK' if ok else 'FAILED'} in {elapsed:.2f}s")


def probe_a():
    """Spread-force footprint vs real omega, on the existing steady baseline."""
    state = c.load_state(BASELINE_RUN, NSTEPS_BASELINE)
    grid = c.Grid(NX, NY, 1, c.DOMAIN["length"], c.DOMAIN["xoffset"], c.DOMAIN["yoffset"])
    geom = c.Geometry(str(c.BASE_GEOM_DX002))
    model = c.NavierStokesModel(grid, geom, c.RE)
    model.init()

    spread = c.Scalar(grid)
    model.B(state.f, spread)
    spread_arr = spread._data[0].copy()
    real_arr = state.omega._data[0].copy()

    X, Y = c.grid_xy(DX)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]
    xs = X[:, 0]
    ys = Y[0, :]
    iy = np.where((ys >= YLIM[0]) & (ys <= YLIM[1]))[0]

    # (1) compact-support extent: how far upstream does the RAW spread
    # footprint have any nonzero value at all, before any buffer exclusion?
    support_cells = 0
    for buf in range(0, 20):
        m = xs <= x_le - buf * DX
        if np.abs(spread_arr[m]).max() == 0.0:
            support_cells = buf - 1
            break
    else:
        support_cells = 19

    # (2) far-upstream comparison, using this folder's standard buffer: the
    # real field is expected to still show noise there (Tests 1-3); if the
    # raw spread footprint is EXACTLY zero there too, that alone proves the
    # far-upstream signal cannot be a direct footprint of B() and must be
    # solve-mediated -- reported as a flag rather than a (meaningless, NaN)
    # correlation against an all-zero array.
    mask = c.upstream_mask(X, x_le, BUFFER_DX, DX)
    real_up = real_arr[np.ix_(np.where(mask)[0], iy)]
    spread_up = spread_arr[np.ix_(np.where(mask)[0], iy)]
    spread_is_zero_far_upstream = bool(np.abs(spread_up).max() == 0.0)
    corr_far = (float(np.corrcoef(real_up.ravel(), spread_up.ravel())[0, 1])
                if not spread_is_zero_far_upstream else float("nan"))

    # (3) near-LE comparison, restricted to where the spread footprint is
    # actually nonzero (0 to support_cells+1 cells upstream) -- the
    # meaningful version of the correlation check.
    near_mask = (xs <= x_le) & (xs >= x_le - (support_cells + 1) * DX)
    real_near = real_arr[np.ix_(np.where(near_mask)[0], iy)]
    spread_near = spread_arr[np.ix_(np.where(near_mask)[0], iy)]
    corr_near = float(np.corrcoef(real_near.ravel(), spread_near.ravel())[0, 1])

    return dict(X=X, Y=Y, real_arr=real_arr, spread_arr=spread_arr, mask=mask, iy=iy,
                support_cells=support_cells, spread_is_zero_far_upstream=spread_is_zero_far_upstream,
                corr_far=corr_far, corr_near=corr_near)


def probe_b():
    om1 = c.load_omega(ONE_STEP_OUTDIR, 1)
    X, Y = c.grid_xy(DX)
    g = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le = g["x"][g["i_le"]]
    m = c.upstream_scalar_metrics(X, Y, om1, x_le, DX, buffer_dx=BUFFER_DX, ylim=YLIM)
    return X, Y, om1, m


def main():
    a = probe_a()
    X, Y, real_arr, spread_arr = a["X"], a["Y"], a["real_arr"], a["spread_arr"]
    print(f"Probe A: raw spread-force footprint has nonzero support out to "
          f"{a['support_cells']} cells upstream of the LE, then is EXACTLY zero "
          f"(spread_is_zero_far_upstream={a['spread_is_zero_far_upstream']})")
    print(f"Probe A: correlation(real, spread) in the near-LE support region = {a['corr_near']:.4f}")
    print(f"Probe A: correlation(real, spread) in the far-upstream (buffer={BUFFER_DX}) window = "
          f"{a['corr_far']}")

    Xb, Yb, om1, m1 = probe_b()
    print(f"Probe B: after nsteps=1 from zero IC, upstream enstrophy={m1['enstrophy']:.6e}, "
          f"peak={m1['peak']:.4f}")

    with open(c.DATA / "test8_mechanism.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["probe", "quantity", "value"])
        w.writeheader()
        w.writerow(dict(probe="A", quantity="spread_support_cells", value=a["support_cells"]))
        w.writerow(dict(probe="A", quantity="spread_is_zero_far_upstream", value=a["spread_is_zero_far_upstream"]))
        w.writerow(dict(probe="A", quantity="correlation_near_LE_support", value=a["corr_near"]))
        w.writerow(dict(probe="A", quantity="correlation_far_upstream", value=a["corr_far"]))
        w.writerow(dict(probe="B", quantity="one_step_upstream_enstrophy", value=m1["enstrophy"]))
        w.writerow(dict(probe="B", quantity="one_step_upstream_peak", value=m1["peak"]))
    print("wrote test8_mechanism.csv")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    xw, yw, real_w, _, _ = c.window(X, Y, real_arr, (-0.35, 0.35), YLIM)
    _, _, spread_w, _, _ = c.window(X, Y, spread_arr, (-0.35, 0.35), YLIM)
    V = 8.0
    axes[0, 0].pcolormesh(xw, yw, np.clip(real_w[:-1, :-1], -V, V), shading="flat", cmap="jet", vmin=-V, vmax=V)
    axes[0, 0].set_title("Real omega (converged steady state)", fontsize=10); axes[0, 0].set_aspect("equal")
    axes[0, 1].pcolormesh(xw, yw, np.clip(spread_w[:-1, :-1], -V, V), shading="flat", cmap="jet", vmin=-V, vmax=V)
    axes[0, 1].set_title(f"Curl(spread boundary force) alone\nzero beyond {a['support_cells']} cells upstream "
                          f"(near-LE corr={a['corr_near']:.2f})", fontsize=10)
    axes[0, 1].set_aspect("equal")

    xw1, yw1, om1_w, _, _ = c.window(Xb, Yb, om1, (-0.35, 0.35), YLIM)
    V1 = max(np.abs(om1_w).max(), 1e-6)
    axes[1, 0].pcolormesh(xw1, yw1, np.clip(om1_w[:-1, :-1], -V1, V1), shading="flat", cmap="jet", vmin=-V1, vmax=V1)
    axes[1, 0].set_title("omega after nsteps=1 from zero IC\n(is upstream pattern already present?)", fontsize=10)
    axes[1, 0].set_aspect("equal")

    ax = axes[1, 1]
    ys_row = Yb[0, :]
    iy0 = c.nearest_index(ys_row, 0.0)
    g2 = c.load_geom_points(c.BASE_GEOM_DX002)
    x_le2 = g2["x"][g2["i_le"]]
    m_up = c.upstream_mask(Xb, x_le2, BUFFER_DX, DX)
    xs_row = Xb[m_up, 0]
    vals_row = om1[np.ix_(np.where(m_up)[0], [iy0])].flatten()
    order = np.argsort(xs_row)
    ax.plot(xs_row[order], vals_row[order], "-o", ms=3, color="#c0392b")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("x"); ax.set_ylabel("omega at y=0")
    ax.set_title("Upstream y=0 row after nsteps=1 (raw values)", fontsize=10)
    ax.grid(alpha=0.3)

    fig.suptitle("Test 8: mechanism isolation -- regularization/projection vs. advection", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = c.FIGS / "test8_mechanism.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        do_run()
    main()
