"""Unsupervised community recovery and agreement with external labels."""

from __future__ import annotations

from collections.abc import Sequence

import networkx as nx
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, normalized_mutual_info_score


def louvain_partition(
    graph: nx.Graph,
    *,
    seed: int = 42,
    resolution: float = 1.0,
    weight: str | None = "score",
) -> dict[str, int]:
    """Infer a complete node partition using NetworkX's Louvain algorithm."""

    if resolution <= 0:
        raise ValueError("resolution must be positive")
    if graph.number_of_nodes() == 0:
        raise ValueError("graph must contain at least one node")
    if graph.number_of_edges() == 0:
        return {str(node): index for index, node in enumerate(sorted(graph.nodes))}

    communities = nx.community.louvain_communities(
        graph,
        weight=weight,
        resolution=resolution,
        seed=seed,
    )
    ordered = sorted(
        (sorted(str(node) for node in community) for community in communities),
        key=lambda nodes: nodes[0],
    )
    partition = {
        node: community_index
        for community_index, nodes in enumerate(ordered)
        for node in nodes
    }
    if set(partition) != {str(node) for node in graph.nodes}:
        raise RuntimeError("inferred communities do not cover every graph node")
    return partition


def partition_agreement(
    partition: dict[str, int],
    node_metadata: pd.DataFrame,
    *,
    label_columns: Sequence[str] = ("family_accession", "clan_accession"),
) -> pd.DataFrame:
    """Compare an inferred partition with one or more unused reference labels."""

    required = {"node_id", *label_columns}
    missing = required - set(node_metadata.columns)
    if missing:
        raise ValueError(
            f"node metadata is missing columns: {', '.join(sorted(missing))}"
        )
    if node_metadata["node_id"].duplicated().any():
        raise ValueError("node metadata identifiers must be unique")
    indexed = node_metadata.set_index("node_id")
    node_order = sorted(partition)
    missing_nodes = set(node_order) - set(indexed.index)
    if missing_nodes:
        raise ValueError(
            f"node metadata is missing nodes: {', '.join(sorted(missing_nodes))}"
        )
    inferred = [partition[node] for node in node_order]
    rows = []
    for label_column in label_columns:
        reference = indexed.loc[node_order, label_column]
        if reference.isna().any():
            raise ValueError(f"reference labels contain missing {label_column} values")
        rows.append(
            {
                "reference_label": label_column,
                "reference_group_count": int(reference.nunique()),
                "inferred_community_count": len(set(inferred)),
                "nmi": normalized_mutual_info_score(
                    reference, inferred, average_method="arithmetic"
                ),
                "ami": adjusted_mutual_info_score(
                    reference, inferred, average_method="arithmetic"
                ),
            }
        )
    return pd.DataFrame(rows)
