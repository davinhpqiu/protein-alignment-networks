"""Convert alignment-score matrices into protein networks."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

METHOD_SCORE_COLUMNS = {
    "biopython": "biopython_score",
    "blast": "blast_bit_score",
    "dedal": "dedal_sw_score",
}
THRESHOLD_RULES = {"absolute", "percentile", "target_density", "top_k", "top_n"}


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


def build_similarity_graph(
    pair_scores: pd.DataFrame,
    protein_ids: Iterable[str],
    *,
    method: str,
    threshold: float,
    threshold_rule: str = "absolute",
    score_column: str | None = None,
    node_metadata: pd.DataFrame | None = None,
    node_id_column: str | None = None,
) -> nx.Graph:
    """Build one method-specific graph from a canonical pair-score table.

    ``threshold_rule`` may be ``absolute`` (native score at least the supplied
    threshold), ``percentile`` (within-method percentile at least a value in
    ``[0, 1]``), ``target_density`` (the strongest pairs up to an exact fraction
    of all possible node pairs), ``top_n`` (the strongest exact number of
    available pairs), or ``top_k`` (the union of each protein's strongest ``k``
    available pairs). Missing method outputs are ineligible edges and are never
    replaced by zero. A target-density graph records any edge-count shortfall
    caused by insufficient available method scores.
    """

    identifiers = list(protein_ids)
    if not identifiers:
        raise ValueError("at least one protein identifier is required")
    if any(not isinstance(identifier, str) or not identifier for identifier in identifiers):
        raise ValueError("protein identifiers must be non-empty strings")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("protein identifiers must be unique")
    if threshold_rule not in THRESHOLD_RULES:
        raise ValueError(
            "threshold_rule must be one of: " + ", ".join(sorted(THRESHOLD_RULES))
        )
    if score_column is None:
        try:
            score_column = METHOD_SCORE_COLUMNS[method]
        except KeyError as error:
            raise ValueError(
                "unknown method; provide score_column or choose one of: "
                + ", ".join(sorted(METHOD_SCORE_COLUMNS))
            ) from error
    required = {"protein_a", "protein_b", score_column}
    missing = required - set(pair_scores.columns)
    if missing:
        raise ValueError(f"missing pair-score column(s): {', '.join(sorted(missing))}")

    pairs = _canonical_pair_scores(pair_scores)
    unknown = (set(pairs["protein_a"]) | set(pairs["protein_b"])) - set(identifiers)
    if unknown:
        raise ValueError(f"pair table contains unknown protein(s): {', '.join(sorted(unknown))}")
    if not pd.api.types.is_numeric_dtype(pairs[score_column]):
        raise TypeError(f"score column must be numeric: {score_column}")

    available = pairs.loc[pairs[score_column].notna()].copy()
    possible_pair_count = len(identifiers) * (len(identifiers) - 1) // 2
    selected = _select_edges(
        available,
        score_column=score_column,
        threshold=threshold,
        threshold_rule=threshold_rule,
        possible_pair_count=possible_pair_count,
    )

    graph_attributes = {
        "method": method,
        "score_column": score_column,
        "threshold_rule": threshold_rule,
        "threshold": float(threshold),
        "available_pair_count": len(available),
        "possible_pair_count": possible_pair_count,
    }
    if threshold_rule == "target_density":
        target_edge_count = int(np.ceil(float(threshold) * possible_pair_count))
        graph_attributes["target_edge_count"] = target_edge_count
        graph_attributes["selection_shortfall"] = max(
            0, target_edge_count - len(selected)
        )
    graph = nx.Graph(**graph_attributes)
    graph.add_nodes_from(identifiers)
    _attach_node_metadata(
        graph,
        node_metadata=node_metadata,
        node_id_column=node_id_column,
    )

    method_prefix = f"{method}_"
    method_columns = [
        column
        for column in selected.columns
        if column.startswith(method_prefix) and column not in {"protein_a", "protein_b"}
    ]
    for row in selected.itertuples(index=False):
        values = row._asdict()
        attributes = {
            "method": method,
            "score_column": score_column,
            "score": float(values[score_column]),
        }
        for column in method_columns:
            value = values[column]
            if not pd.isna(value):
                attributes[column] = _python_scalar(value)
        graph.add_edge(values["protein_a"], values["protein_b"], **attributes)

    graph.graph["selected_edge_count"] = graph.number_of_edges()
    return graph


def graph_summary(graph: nx.Graph) -> dict[str, float | int | str]:
    """Return small, interpretable structural summaries for one graph."""

    if graph.is_directed():
        raise ValueError("graph_summary expects an undirected graph")
    possible_pair_count = int(graph.graph.get("possible_pair_count", 0))
    available_pair_count = int(graph.graph.get("available_pair_count", 0))
    return {
        "method": str(graph.graph.get("method", "")),
        "threshold_rule": str(graph.graph.get("threshold_rule", "")),
        "threshold": float(graph.graph.get("threshold", float("nan"))),
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "isolates": nx.number_of_isolates(graph),
        "connected_components": nx.number_connected_components(graph),
        "density": nx.density(graph),
        "average_clustering": nx.average_clustering(graph),
        "possible_pairs": possible_pair_count,
        "available_pairs": available_pair_count,
        "available_pair_fraction": (
            available_pair_count / possible_pair_count
            if possible_pair_count
            else float("nan")
        ),
        "target_edges": int(graph.graph.get("target_edge_count", -1)),
        "selection_shortfall": int(graph.graph.get("selection_shortfall", 0)),
    }


def save_graph_bundle(graph: nx.Graph, output_stem: str | Path) -> dict[str, Path]:
    """Save GraphML, edge/node tables, and a JSON structural summary."""

    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = {
        "graphml": stem.with_suffix(".graphml"),
        "edges": stem.with_suffix(".edges.tsv"),
        "nodes": stem.with_suffix(".nodes.tsv"),
        "summary": stem.with_suffix(".summary.json"),
    }
    nx.write_graphml(graph, paths["graphml"])

    edge_rows = [
        {"protein_a": a, "protein_b": b, **attributes}
        for a, b, attributes in sorted(graph.edges(data=True))
    ]
    edge_columns = ["protein_a", "protein_b", "method", "score_column", "score"]
    pd.DataFrame(edge_rows, columns=edge_columns if not edge_rows else None).to_csv(
        paths["edges"], sep="\t", index=False
    )
    node_rows = [
        {"protein_id": identifier, **attributes}
        for identifier, attributes in sorted(graph.nodes(data=True))
    ]
    pd.DataFrame(node_rows).to_csv(paths["nodes"], sep="\t", index=False)
    paths["summary"].write_text(
        json.dumps(graph_summary(graph), indent=2) + "\n", encoding="utf-8"
    )
    return paths


def _canonical_pair_scores(pair_scores: pd.DataFrame) -> pd.DataFrame:
    pairs = pair_scores.copy()
    endpoints = pairs[["protein_a", "protein_b"]].astype(str)
    pairs["protein_a"] = endpoints.min(axis=1)
    pairs["protein_b"] = endpoints.max(axis=1)
    if (pairs["protein_a"] == pairs["protein_b"]).any():
        raise ValueError("pair table cannot contain self-pairs")
    if pairs.duplicated(["protein_a", "protein_b"]).any():
        raise ValueError("pair table contains duplicate protein pairs")
    return pairs.sort_values(["protein_a", "protein_b"], ignore_index=True)


def _select_edges(
    available: pd.DataFrame,
    *,
    score_column: str,
    threshold: float,
    threshold_rule: str,
    possible_pair_count: int,
) -> pd.DataFrame:
    if threshold_rule == "absolute":
        if not isinstance(threshold, (int, float)) or not np.isfinite(threshold):
            raise ValueError("absolute threshold must be a finite number")
        return available.loc[available[score_column] >= float(threshold)]
    if threshold_rule == "percentile":
        if not isinstance(threshold, (int, float)) or not 0.0 <= threshold <= 1.0:
            raise ValueError("percentile threshold must be between 0 and 1")
        percentiles = available[score_column].rank(
            method="average", ascending=True, pct=True
        )
        return available.loc[percentiles >= float(threshold)]

    if threshold_rule == "target_density":
        if not isinstance(threshold, (int, float)) or not 0.0 <= threshold <= 1.0:
            raise ValueError("target_density threshold must be between 0 and 1")
        target_edge_count = int(np.ceil(float(threshold) * possible_pair_count))
        return available.sort_values(
            [score_column, "protein_a", "protein_b"],
            ascending=[False, True, True],
            kind="stable",
        ).head(target_edge_count)

    if threshold_rule == "top_k":
        if not isinstance(threshold, int) or isinstance(threshold, bool):
            raise TypeError("top_k threshold must be an integer")
        if threshold < 0:
            raise ValueError("top_k threshold must be non-negative")
        ordered = available.sort_values(
            [score_column, "protein_a", "protein_b"],
            ascending=[False, True, True],
            kind="stable",
        )
        if threshold == 0 or ordered.empty:
            return ordered.iloc[0:0]
        endpoint_a = ordered.assign(_protein=ordered["protein_a"])
        endpoint_b = ordered.assign(_protein=ordered["protein_b"])
        endpoint_rows = pd.concat([endpoint_a, endpoint_b]).sort_values(
            [score_column, "protein_a", "protein_b"],
            ascending=[False, True, True],
            kind="stable",
        )
        selected_indices = endpoint_rows.groupby("_protein", sort=True).head(
            threshold
        ).index.unique()
        return ordered.loc[ordered.index.isin(selected_indices)]

    if not isinstance(threshold, int) or isinstance(threshold, bool):
        raise TypeError("top_n threshold must be an integer")
    if not 0 <= threshold <= len(available):
        raise ValueError("top_n threshold must be between zero and available pairs")
    return available.sort_values(
        [score_column, "protein_a", "protein_b"],
        ascending=[False, True, True],
        kind="stable",
    ).head(threshold)


def _attach_node_metadata(
    graph: nx.Graph,
    *,
    node_metadata: pd.DataFrame | None,
    node_id_column: str | None,
) -> None:
    if node_metadata is None:
        return
    if node_id_column is None:
        raise ValueError("node_id_column is required when node_metadata is supplied")
    if node_id_column not in node_metadata:
        raise ValueError(f"node metadata is missing ID column: {node_id_column}")
    if node_metadata[node_id_column].duplicated().any():
        raise ValueError("node metadata identifiers must be unique")

    indexed = node_metadata.set_index(node_id_column, drop=True)
    unknown = set(indexed.index.astype(str)) - set(graph.nodes)
    if unknown:
        raise ValueError(f"node metadata contains unknown protein(s): {', '.join(sorted(unknown))}")
    for identifier, row in indexed.iterrows():
        attributes = {
            column: _python_scalar(value)
            for column, value in row.items()
            if not pd.isna(value)
        }
        graph.nodes[str(identifier)].update(attributes)


def _python_scalar(value):
    return value.item() if hasattr(value, "item") else value
