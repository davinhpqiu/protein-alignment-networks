from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from protein_alignment_networks.graphs import (
    build_similarity_graph,
    graph_summary,
    save_graph_bundle,
)


@pytest.fixture
def pair_scores():
    return pd.DataFrame(
        {
            "protein_a": ["a", "a", "a", "b", "b", "c"],
            "protein_b": ["b", "c", "d", "c", "d", "d"],
            "biopython_score": [9.0, 7.0, 2.0, 6.0, 3.0, 1.0],
            "blast_bit_score": [90.0, np.nan, 20.0, 60.0, np.nan, 10.0],
            "dedal_sw_score": [8.0, 5.0, 4.0, 7.0, 2.0, 1.0],
        }
    )


def test_absolute_threshold_builds_separate_method_graphs_and_keeps_isolates(
    pair_scores,
):
    biopython = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="biopython",
        threshold=7.0,
    )
    blast = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="blast",
        threshold=50.0,
    )
    assert set(biopython.edges) == {("a", "b"), ("a", "c")}
    assert set(blast.edges) == {("a", "b"), ("b", "c")}
    assert "d" in biopython
    assert biopython["a"]["b"]["score"] == 9.0


def test_missing_blast_hits_are_ineligible_not_zero(pair_scores):
    graph = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="blast",
        threshold=-1.0,
    )
    assert set(graph.edges) == {("a", "b"), ("a", "d"), ("b", "c"), ("c", "d")}
    assert ("a", "c") not in graph.edges


def test_percentile_target_density_top_n_and_top_k_rules_are_toggleable(pair_scores):
    percentile = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="dedal",
        threshold=0.75,
        threshold_rule="percentile",
    )
    top_n = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="dedal",
        threshold=2,
        threshold_rule="top_n",
    )
    top_k = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="dedal",
        threshold=1,
        threshold_rule="top_k",
    )
    target_density = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="dedal",
        threshold=0.5,
        threshold_rule="target_density",
    )
    assert percentile.number_of_edges() == 2
    assert set(top_n.edges) == {("a", "b"), ("b", "c")}
    assert set(top_k.edges) == {("a", "b"), ("a", "d"), ("b", "c")}
    assert set(target_density.edges) == {("a", "b"), ("a", "c"), ("b", "c")}
    assert graph_summary(target_density)["target_edges"] == 3
    assert not list(nx.isolates(top_k))


def test_target_density_records_missing_score_shortfall(pair_scores):
    graph = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="blast",
        threshold=1.0,
        threshold_rule="target_density",
    )
    summary = graph_summary(graph)
    assert graph.number_of_edges() == 4
    assert summary["target_edges"] == 6
    assert summary["selection_shortfall"] == 2
    assert summary["available_pair_fraction"] == pytest.approx(4 / 6)


def test_node_metadata_is_attached(pair_scores):
    metadata = pd.DataFrame(
        {
            "accession": ["a", "b", "c", "d"],
            "subgroup": ["x", "x", "y", "y"],
            "length": [100, 101, 102, 103],
        }
    )
    graph = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="biopython",
        threshold=2,
        node_metadata=metadata,
        node_id_column="accession",
    )
    assert graph.nodes["c"] == {"subgroup": "y", "length": 102}


def test_graph_bundle_round_trip(tmp_path: Path, pair_scores):
    graph = build_similarity_graph(
        pair_scores,
        ["a", "b", "c", "d"],
        method="biopython",
        threshold=2,
    )
    paths = save_graph_bundle(graph, tmp_path / "example")
    restored = nx.read_graphml(paths["graphml"])
    assert set(restored.nodes) == set(graph.nodes)
    assert set(restored.edges) == set(graph.edges)
    assert pd.read_csv(paths["edges"], sep="\t").shape[0] == graph.number_of_edges()
    assert graph_summary(graph)["connected_components"] == 1


def test_invalid_graph_options_are_rejected(pair_scores):
    with pytest.raises(ValueError, match="between 0 and 1"):
        build_similarity_graph(
            pair_scores,
            ["a", "b", "c", "d"],
            method="dedal",
            threshold=2,
            threshold_rule="percentile",
        )
