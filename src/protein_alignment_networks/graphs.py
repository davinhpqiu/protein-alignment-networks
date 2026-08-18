"""Convert alignment-score matrices into protein networks."""

from __future__ import annotations

import networkx as nx
import pandas as pd


def score_matrix_to_graph(
    score_matrix: pd.DataFrame,
    *,
    threshold: float,
    include_self_loops: bool = False,
) -> nx.Graph:
    """Create an undirected graph containing edges with score >= threshold.

    All proteins become nodes, including isolates. Each retained edge receives
    a numeric ``score`` attribute.
    """

    if score_matrix.empty:
        raise ValueError("score_matrix cannot be empty")
    if list(score_matrix.index) != list(score_matrix.columns):
        raise ValueError("score_matrix rows and columns must have identical labels")
    if not score_matrix.equals(score_matrix.T):
        raise ValueError("score_matrix must be symmetric")

    graph = nx.Graph(threshold=float(threshold))
    labels = list(score_matrix.index)
    graph.add_nodes_from(labels)
    for i, label_a in enumerate(labels):
        start = i if include_self_loops else i + 1
        for j in range(start, len(labels)):
            score = float(score_matrix.iat[i, j])
            if score >= threshold:
                graph.add_edge(label_a, labels[j], score=score)
    return graph

