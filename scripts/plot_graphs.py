from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import nxviz

from ribo_switch.graph import build_constraint_graph
from ribo_switch.structure import parse_dot_bracket


DOT_BRACKET_CHARS = {".", "(", ")"}


def _is_dot_bracket(line: str) -> bool:
    return bool(line) and set(line) <= DOT_BRACKET_CHARS


def _read_structure_pair(input_path: Path) -> tuple[str, str]:
    lines = [line.strip() for line in input_path.read_text(encoding="utf-8").splitlines()]
    structures = [line for line in lines if _is_dot_bracket(line)]
    if len(structures) < 2:
        raise ValueError("Input file must contain at least two dot-bracket lines")
    return (structures[0], structures[1])


def _collect_edge_sources(s_on: str, s_off: str, mode: str) -> tuple[int, dict[tuple[int, int], set[str]]]:
    struct_on = parse_dot_bracket(s_on)
    struct_off = parse_dot_bracket(s_off)
    cgraph = build_constraint_graph(struct_on, struct_off)
    edge_sources: dict[tuple[int, int], set[str]] = {}
    for comp in cgraph.components:
        for edge in comp.edges:
            if mode == "on" and edge.source != "on":
                continue
            if mode == "off" and edge.source != "off":
                continue
            i, j = (edge.i, edge.j)
            if i > j:
                i, j = (j, i)
            edge_sources.setdefault((i, j), set()).add(edge.source)
    return (cgraph.n, edge_sources)


def _build_nx_graph(s_on: str, s_off: str, mode: str) -> nx.Graph:
    n, edge_sources = _collect_edge_sources(s_on, s_off, mode)
    graph = nx.Graph()
    for node in range(n):
        graph.add_node(node, index=node + 1)
    for (i, j), sources in edge_sources.items():
        if mode == "combined" and len(sources) > 1:
            source = "both"
        else:
            source = next(iter(sources))
        graph.add_edge(i, j, source=source)
    return graph


def _create_plot(plot_fn, graph: nx.Graph):
    try:
        return plot_fn(
            graph,
            sort_by="index",
            edge_color_by="source",
            edge_palette="Set2",
        )
    except TypeError:
        try:
            return plot_fn(graph, sort_by="index")
        except TypeError:
            return plot_fn(graph)


def _save_plot(plot_fn, graph: nx.Graph, out_path: Path, dpi: int) -> None:
    plt.figure(figsize=(8, 8))
    _create_plot(plot_fn, graph)
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate arc and circos plots from a riboswitch constraint graph.")
    parser.add_argument("--on", dest="s_on", help="S_ON dot-bracket string")
    parser.add_argument("--off", dest="s_off", help="S_OFF dot-bracket string")
    parser.add_argument("--input", "-i", dest="input_file", help="File with two dot-bracket lines")
    parser.add_argument("--out-dir", default=".", help="Output directory")
    parser.add_argument("--prefix", default="graph", help="Output file prefix")
    parser.add_argument(
        "--graph",
        choices=["combined", "on", "off"],
        default="combined",
        help="Which base-pair graph to render",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output DPI")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.input_file:
        s_on, s_off = _read_structure_pair(Path(args.input_file))
    elif args.s_on and args.s_off:
        s_on, s_off = (args.s_on, args.s_off)
    else:
        raise SystemExit("Provide --on/--off or --input")

    if len(s_on) != len(s_off):
        raise SystemExit("S_ON and S_OFF must have the same length")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    graph = _build_nx_graph(s_on, s_off, args.graph)
    arc_path = out_dir / f"{args.prefix}_arc.png"
    circos_path = out_dir / f"{args.prefix}_circos.png"

    _save_plot(nxviz.arc, graph, arc_path, args.dpi)
    _save_plot(nxviz.circos, graph, circos_path, args.dpi)

    print(f"Wrote: {arc_path}")
    print(f"Wrote: {circos_path}")


if __name__ == "__main__":
    main()
