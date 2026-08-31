"""Consistent visualizations for method-specific protein graphs."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D


def shared_graph_layout(
    graphs: Iterable[nx.Graph],
    *,
    seed: int = 42,
    iterations: int = 100,
) -> dict[str, tuple[float, float]]:
    """Lay out the union of several graphs for fair side-by-side comparison."""

    graph_list = list(graphs)
    if not graph_list:
        raise ValueError("at least one graph is required")
    union = nx.Graph()
    for graph in graph_list:
        union.add_nodes_from(graph.nodes)
        union.add_edges_from(graph.edges)
    if not union:
        raise ValueError("graphs must contain at least one node")
    return nx.spring_layout(
        union,
        seed=seed,
        iterations=iterations,
        weight=None,
    )


def draw_similarity_graph(
    graph: nx.Graph,
    *,
    positions: Mapping[str, tuple[float, float]] | None = None,
    color_attribute: str | None = None,
    category_colors: Mapping[str, str] | None = None,
    title: str | None = None,
    ax: Axes | None = None,
    node_size: float = 45,
    edge_alpha: float = 0.18,
    with_labels: bool = False,
) -> tuple[Figure, Axes]:
    """Draw a readable protein graph with optional categorical node colours."""

    if graph.number_of_nodes() == 0:
        raise ValueError("graph must contain at least one node")
    if positions is None:
        positions = nx.spring_layout(graph, seed=42, weight=None)
    missing_positions = set(graph) - set(positions)
    if missing_positions:
        raise ValueError(
            "positions are missing graph node(s): "
            + ", ".join(sorted(missing_positions))
        )
    if ax is None:
        figure, ax = plt.subplots(figsize=(8, 7))
    else:
        figure = ax.figure

    colors: list[str] | str = "#4C78A8"
    legend_handles: list[Line2D] = []
    if color_attribute is not None:
        categories = [
            str(graph.nodes[node].get(color_attribute, "Unknown")) for node in graph
        ]
        unique_categories = sorted(set(categories))
        if category_colors is None:
            color_map = plt.get_cmap("tab20", max(1, len(unique_categories)))
            category_colors = {
                category: color_map(index)
                for index, category in enumerate(unique_categories)
            }
        missing_colors = set(unique_categories) - set(category_colors)
        if missing_colors:
            raise ValueError(
                "category_colors is missing value(s): "
                + ", ".join(sorted(missing_colors))
            )
        colors = [category_colors[category] for category in categories]
        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=category_colors[category],
                markeredgecolor="none",
                markersize=7,
                label=category,
            )
            for category in unique_categories
        ]

    nx.draw_networkx_edges(
        graph,
        positions,
        ax=ax,
        alpha=edge_alpha,
        width=0.7,
        edge_color="#65737E",
    )
    nx.draw_networkx_nodes(
        graph,
        positions,
        ax=ax,
        node_color=colors,
        node_size=node_size,
        linewidths=0.25,
        edgecolors="white",
    )
    if with_labels:
        nx.draw_networkx_labels(graph, positions, ax=ax, font_size=6)
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            title=color_attribute.replace("_", " ").title(),
            loc="best",
            frameon=False,
            fontsize=7,
        )
    ax.set_title(title or str(graph.graph.get("method", "Protein graph")))
    ax.set_axis_off()
    return figure, ax
