import networkx as nx
import pandas as pd
import pytest

from protein_alignment_networks.evaluation import (
    louvain_partition,
    partition_agreement,
)


def test_partition_agreement_is_one_for_the_reference_partition():
    metadata = pd.DataFrame(
        {
            "node_id": ["a", "b", "c", "d"],
            "family_accession": ["F1", "F1", "F2", "F2"],
            "clan_accession": ["C1", "C1", "C2", "C2"],
        }
    )
    agreement = partition_agreement(
        {"a": 0, "b": 0, "c": 1, "d": 1}, metadata
    ).set_index("reference_label")

    assert agreement.loc["family_accession", "nmi"] == pytest.approx(1.0)
    assert agreement.loc["family_accession", "ami"] == pytest.approx(1.0)
    assert agreement.loc["clan_accession", "nmi"] == pytest.approx(1.0)


def test_louvain_partition_covers_isolates_and_is_reproducible():
    graph = nx.Graph()
    graph.add_edge("a", "b", score=3.0)
    graph.add_edge("c", "d", score=3.0)
    graph.add_node("e")

    first = louvain_partition(graph, seed=9)
    second = louvain_partition(graph, seed=9)

    assert first == second
    assert set(first) == set(graph)
    assert first["a"] == first["b"]
    assert first["c"] == first["d"]


def test_edgeless_graph_assigns_each_node_to_a_singleton():
    graph = nx.empty_graph(["b", "a"])
    assert louvain_partition(graph) == {"a": 0, "b": 1}
