import matplotlib
import networkx as nx
import pytest

matplotlib.use("Agg")

from protein_alignment_networks.visualization import (
    draw_similarity_graph,
    shared_graph_layout,
)


def test_shared_layout_contains_union_of_graph_nodes():
    first = nx.Graph([("a", "b")])
    second = nx.Graph([("b", "c")])
    positions = shared_graph_layout([first, second], seed=7)
    assert set(positions) == {"a", "b", "c"}


def test_draw_similarity_graph_uses_metadata_categories():
    graph = nx.Graph([("a", "b")])
    graph.nodes["a"]["family"] = "x"
    graph.nodes["b"]["family"] = "y"
    positions = {"a": (0.0, 0.0), "b": (1.0, 0.0)}
    figure, axis = draw_similarity_graph(
        graph,
        positions=positions,
        color_attribute="family",
        category_colors={"x": "red", "y": "blue"},
    )
    assert figure is axis.figure
    assert axis.axison is False


def test_draw_similarity_graph_rejects_incomplete_positions():
    graph = nx.Graph([("a", "b")])
    with pytest.raises(ValueError, match="missing graph node"):
        draw_similarity_graph(graph, positions={"a": (0.0, 0.0)})
