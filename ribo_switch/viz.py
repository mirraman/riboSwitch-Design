"""
Server-side rendering for constraint-graph (arc/circos via nxviz) and
RNA structure (VARNA) visualizations.  All public functions return PNG
bytes; callers base64-encode them for JSON transport.

VARNA_JAR env var must point to VARNAv3-93.jar (or compatible).
"""
from __future__ import annotations

import io
import os
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as _cm
# nxviz 0.9+ calls matplotlib.cm.get_cmap which was removed in matplotlib 3.9
if not hasattr(_cm, "get_cmap"):
    _cm.get_cmap = lambda name, lut=None: matplotlib.colormaps[name]
import networkx as nx
import nxviz

from ribo_switch.graph import build_constraint_graph
from ribo_switch.structure import parse_dot_bracket

_VARNA_JAR = os.environ.get("VARNA_JAR", "")
_VARNA_CLASS = "fr.orsay.lri.varna.applications.VARNAcmd"


# ── constraint-graph helpers ───────────────────────────────────────────────

def _build_nx_graph(s_on: str, s_off: str, mode: str) -> nx.Graph:
    struct_on = parse_dot_bracket(s_on)
    struct_off = parse_dot_bracket(s_off)
    cgraph = build_constraint_graph(struct_on, struct_off)
    n = cgraph.n

    edge_sources: dict[tuple[int, int], set[str]] = {}
    for comp in cgraph.components:
        for edge in comp.edges:
            if mode == "on" and edge.source != "on":
                continue
            if mode == "off" and edge.source != "off":
                continue
            i, j = (min(edge.i, edge.j), max(edge.i, edge.j))
            edge_sources.setdefault((i, j), set()).add(edge.source)

    graph = nx.Graph()
    for node in range(n):
        graph.add_node(node, index=node + 1)
    for (i, j), sources in edge_sources.items():
        source = "both" if (mode == "combined" and len(sources) > 1) else next(iter(sources))
        graph.add_edge(i, j, source=source)
    return graph


def _plot_to_bytes(plot_fn, graph: nx.Graph, dpi: int = 150) -> bytes:
    fig = plt.figure(figsize=(8, 8))
    try:
        plot_fn(graph, sort_by="index", edge_color_by="source", edge_palette="Set2")
    except TypeError:
        try:
            plot_fn(graph, sort_by="index")
        except TypeError:
            plot_fn(graph)
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.read()


# ── public rendering functions ─────────────────────────────────────────────

def render_arc(s_on: str, s_off: str, mode: str = "combined") -> bytes:
    return _plot_to_bytes(nxviz.arc, _build_nx_graph(s_on, s_off, mode))


def render_circos(s_on: str, s_off: str, mode: str = "combined") -> bytes:
    return _plot_to_bytes(nxviz.circos, _build_nx_graph(s_on, s_off, mode))


def render_varna(
    sequence: str,
    structure: str,
    *,
    title: str = "",
    algorithm: str = "radiate",
    resolution: float = 2.0,
) -> bytes:
    """Render an RNA secondary structure to PNG bytes using VARNAcmd.

    Requires VARNA_JAR env var to point at VARNAv3-93.jar and java on PATH.
    """
    jar = _VARNA_JAR
    if not jar or not os.path.isfile(jar):
        raise RuntimeError(
            f"VARNA_JAR env var is not set or file not found: {jar!r}. "
            "Set VARNA_JAR=/path/to/VARNAv3-93.jar"
        )

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        out_path = f.name

    try:
        cmd = [
            "java", "-Djava.awt.headless=true", "-cp", jar, _VARNA_CLASS,
            "-sequenceDBN", sequence,
            "-structureDBN", structure,
            "-algorithm", algorithm,
            "-resolution", str(resolution),
            "-o", out_path,
        ]
        if title:
            cmd += ["-title", title]

        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            timeout=60,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"VARNAcmd failed (exit {result.returncode}): "
                f"{result.stderr.decode(errors='replace')}"
            )

        with open(out_path, "rb") as f:
            return f.read()
    finally:
        try:
            os.unlink(out_path)
        except FileNotFoundError:
            pass
