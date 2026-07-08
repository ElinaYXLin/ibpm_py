"""
generate_module_graph.py

Statically parses every module under py/ with Python's `ast` module and
draws a "which file imports which file" dependency graph for the whole
ibpm_py package.

This is a MECHANICAL, fully-automated diagram: every edge comes directly
from a `from .module import ...` (or `from . import module`) statement
found by `ast.parse`. Nothing here is hand-traced or guessed, so it can be
regenerated at any time as the code changes, and every edge can be checked
by opening the source file and looking at its import block.

Usage:
    python3 flowchart/generate_module_graph.py

Output:
    flowchart/output/module_dependency_graph.png
    flowchart/output/module_dependency_edges.csv   (raw edge list used to draw the figure)
"""

from __future__ import annotations

import ast
import csv
import glob
import os
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY_DIR = os.path.join(REPO_ROOT, "py")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

# --------------------------------------------------------------------------
# 1. Extract edges by parsing every `from .x import y` / `from . import x`
#    statement in py/*.py. This is the only source of truth for this figure.
# --------------------------------------------------------------------------


def extract_edges() -> tuple[set[str], list[tuple[str, str]]]:
    nodes: set[str] = set()
    edges: list[tuple[str, str]] = []
    for path in sorted(glob.glob(os.path.join(PY_DIR, "*.py"))):
        mod = os.path.basename(path)[:-3]
        if mod == "__init__":
            continue
        nodes.add(mod)
        tree = ast.parse(open(path, "r", encoding="utf-8").read(), filename=path)
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 1:
                if node.module:
                    edges.append((mod, node.module.split(".")[0]))
                else:
                    for alias in node.names:
                        edges.append((mod, alias.name))
    return nodes, edges


# --------------------------------------------------------------------------
# 2. Layer nodes automatically by longest dependency path (a topological
#    layering): a module with no in-package imports is layer 0; everything
#    else is 1 + max(layer of everything it imports). This makes the entry
#    point (ibpm.py, which transitively imports the most) float to the top.
# --------------------------------------------------------------------------


def compute_layers(nodes: set[str], edges: list[tuple[str, str]]) -> dict[str, int]:
    deps = defaultdict(list)
    for src, dst in edges:
        deps[src].append(dst)

    layer: dict[str, int] = {}

    def layer_of(n: str, stack: tuple[str, ...] = ()) -> int:
        if n in layer:
            return layer[n]
        if n in stack:  # defensive: no cycles expected in this codebase
            return 0
        if not deps[n]:
            layer[n] = 0
        else:
            layer[n] = 1 + max(layer_of(d, stack + (n,)) for d in deps[n] if d in nodes)
        return layer[n]

    for n in nodes:
        layer_of(n)
    return layer


# --------------------------------------------------------------------------
# 3. Purely cosmetic role grouping (color only) so related files cluster
#    visually. This grouping does NOT affect which edges are drawn -- it is
#    a readability aid on top of the mechanically-extracted graph above.
# --------------------------------------------------------------------------

ROLE_GROUPS = {
    "entry": {"ibpm"},
    "solvers": {
        "ib_solver", "projection_solver", "cholesky_solver",
        "conjugate_gradient_solver", "elliptic_solver", "elliptic_solver_2d",
        "regularizer",
    },
    "model": {"navier_stokes_model", "scheme"},
    "core_data": {"state", "scalar", "flux", "boundary_vector", "field", "grid", "bc"},
    "geometry_motion": {
        "geometry", "rigid_body", "motion", "tangent_se2", "direction",
        "eldredge1", "eldredge2", "eldredge_combined2", "eldredge_maneuver",
        "fixed_position", "fixed_velocity", "lag_step1", "lag_step2",
        "motion_file", "motion_file_periodic", "pitch_plunge",
        "sigmoidal_step", "base_flow",
    },
    "io_utils": {
        "logger", "output", "output_energy", "output_force", "output_probes",
        "output_restart", "output_tecplot", "scalar_to_tecplot",
        "parm_parser", "checkgeom", "utils", "vector_operations",
    },
}

ROLE_COLOR = {
    "entry": "#c0392b",
    "solvers": "#8e44ad",
    "model": "#2980b9",
    "core_data": "#16a085",
    "geometry_motion": "#d35400",
    "io_utils": "#7f8c8d",
}


def role_of(mod: str) -> str:
    for role, members in ROLE_GROUPS.items():
        if mod in members:
            return role
    return "io_utils"


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    nodes, edges = extract_edges()
    layers = compute_layers(nodes, edges)

    # dump the raw edge list so the picture can be checked line-by-line
    edge_csv = os.path.join(OUT_DIR, "module_dependency_edges.csv")
    with open(edge_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["importer", "imports_from"])
        for src, dst in sorted(set(edges)):
            w.writerow([src, dst])

    # group nodes by layer, order alphabetically within a layer for stability
    by_layer: dict[int, list[str]] = defaultdict(list)
    for n in nodes:
        by_layer[layers[n]].append(n)
    for layer_nodes in by_layer.values():
        layer_nodes.sort()

    max_layer = max(by_layer)
    max_width = max(len(v) for v in by_layer.values())

    # --- layout: x by position within layer, y by layer (top = entry point)
    # col_w must comfortably fit the longest label in the widest layer, or
    # neighboring boxes/text overlap -- sized empirically from label length.
    pos: dict[str, tuple[float, float]] = {}
    col_w = 3.6
    row_h = 2.3
    for layer_idx, layer_nodes in by_layer.items():
        n = len(layer_nodes)
        total_w = n * col_w
        x0 = -total_w / 2 + col_w / 2
        y = layer_idx * row_h
        for i, node in enumerate(layer_nodes):
            pos[node] = (x0 + i * col_w, y)

    fig_w = max(22, max_width * col_w * 0.62)
    fig_h = max(13, (max_layer + 1) * row_h * 0.95)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # edges first (so nodes draw on top)
    for src, dst in sorted(set(edges)):
        if src == dst:
            continue
        x1, y1 = pos[src]
        x2, y2 = pos[dst]
        arrow = FancyArrowPatch(
            (x1, y1 - 0.32), (x2, y2 + 0.32),
            connectionstyle="arc3,rad=0.08",
            arrowstyle="-|>", mutation_scale=9,
            color=ROLE_COLOR[role_of(src)], alpha=0.35, linewidth=0.9,
            zorder=1,
        )
        ax.add_patch(arrow)

    box_w, box_h = 3.2, 0.6
    for node, (x, y) in pos.items():
        color = ROLE_COLOR[role_of(node)]
        box = FancyBboxPatch(
            (x - box_w / 2, y - box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            linewidth=1.2, edgecolor=color, facecolor="white",
            zorder=2,
        )
        ax.add_patch(box)
        ax.text(x, y, f"{node}.py", ha="center", va="center", fontsize=7.6,
                 fontweight="bold", color=color, zorder=3)

    # legend
    handles = [
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=c,
                    markeredgecolor=c, markersize=12, label=role.replace("_", " "))
        for role, c in ROLE_COLOR.items()
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.0),
               fontsize=10, title="role (cosmetic grouping only)", frameon=True)

    ax.set_title(
        "ibpm_py/py -- module import dependency graph\n"
        "(auto-generated from `ast`-parsed `from .x import y` statements; "
        "arrow A -> B means \"A.py imports B.py\")",
        fontsize=13,
    )
    ax.set_xlim(-fig_w * 1.05, fig_w * 1.05)
    ax.set_ylim(-1.2, (max_layer + 1) * row_h)
    ax.axis("off")

    out_png = os.path.join(OUT_DIR, "module_dependency_graph.png")
    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    print(f"wrote {out_png}")
    print(f"wrote {edge_csv}")
    print(f"{len(nodes)} nodes, {len(set(edges))} unique edges")


if __name__ == "__main__":
    main()
