"""Run a deterministic score-permutation negative control on sampled collections.

For each selected collection, finite pair scores are shuffled among pair
identities. This preserves the nodes, score distribution, missing-value count,
and graph rule while destroying the relationship between sequence pairs and
scores. Raw graph/community evaluations and a checksum-rich manifest are saved
for Notebook 03; the notebook loads these artifacts rather than copying results.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from protein_alignment_networks import (
    build_similarity_graph,
    graph_summary,
    louvain_partition,
    partition_agreement,
)

METHOD_SCORE_COLUMNS = {
    "biopython": "biopython_score",
    "blast": "blast_bit_score",
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_seed(base_seed: int, *parts: object) -> int:
    payload = "|".join(map(str, (base_seed, *parts))).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def permute_finite_scores(
    pairs: pd.DataFrame,
    score_column: str,
    *,
    seed: int,
) -> pd.DataFrame:
    """Shuffle finite scores while retaining endpoints and missing positions."""

    if score_column not in pairs:
        raise ValueError(f"pair table is missing score column: {score_column}")
    shuffled = pairs.copy()
    finite = shuffled[score_column].notna()
    values = shuffled.loc[finite, score_column].to_numpy(copy=True)
    shuffled.loc[finite, score_column] = np.random.default_rng(seed).permutation(values)
    return shuffled


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("collection_manifest", type=Path)
    parser.add_argument("collection_membership", type=Path)
    parser.add_argument("node_metadata", type=Path)
    parser.add_argument("pair_scores", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--method", choices=sorted(METHOD_SCORE_COLUMNS), default="biopython")
    parser.add_argument("--top-k", nargs="*", type=int, default=[5])
    parser.add_argument("--target-densities", nargs="*", type=float, default=[0.02])
    parser.add_argument("--clans", type=int, default=4)
    parser.add_argument("--families-per-clan", type=int, default=2)
    parser.add_argument("--domains-per-family", type=int, default=20)
    parser.add_argument("--max-collections", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--community-resolution", type=float, default=1.0)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if arguments.max_collections < 1:
        raise ValueError("--max-collections must be positive")
    if arguments.community_resolution <= 0:
        raise ValueError("--community-resolution must be positive")
    if any(value < 0 for value in arguments.top_k):
        raise ValueError("--top-k values must be non-negative")
    if any(not 0 <= value <= 1 for value in arguments.target_densities):
        raise ValueError("--target-densities values must be between zero and one")

    manifest = pd.read_csv(arguments.collection_manifest, sep="\t")
    manifest = manifest.loc[
        (manifest["clans_per_collection"] == arguments.clans)
        & (manifest["families_per_clan"] == arguments.families_per_clan)
        & (manifest["domains_per_family"] == arguments.domains_per_family)
    ].sort_values("collection_id", ignore_index=True)
    manifest = manifest.head(arguments.max_collections)
    if manifest.empty:
        raise ValueError("no collections match the requested composition")

    membership = pd.read_csv(
        arguments.collection_membership, sep="\t", dtype={"node_id": str}
    )
    metadata = pd.read_csv(
        arguments.node_metadata, sep="\t", dtype={"node_id": str}
    )
    pair_scores = pd.read_csv(
        arguments.pair_scores,
        sep="\t",
        dtype={"protein_a": str, "protein_b": str},
    )
    score_column = METHOD_SCORE_COLUMNS[arguments.method]
    if score_column not in pair_scores:
        raise ValueError(f"pair table is missing score column: {score_column}")

    members_by_collection = (
        membership.groupby("collection_id", sort=False)["node_id"].apply(list).to_dict()
    )
    metadata_index = metadata.set_index("node_id", drop=False)
    settings = [("top_k", value) for value in arguments.top_k] + [
        ("target_density", value) for value in arguments.target_densities
    ]
    rows = []
    for collection in manifest.itertuples(index=False):
        identifiers = members_by_collection[collection.collection_id]
        identifier_set = set(identifiers)
        collection_metadata = metadata_index.loc[identifiers].reset_index(drop=True)
        collection_pairs = pair_scores.loc[
            pair_scores["protein_a"].isin(identifier_set)
            & pair_scores["protein_b"].isin(identifier_set)
        ].copy()
        permutation_seed = stable_seed(
            arguments.seed, collection.collection_id, arguments.method, "score_permutation"
        )
        null_pairs = permute_finite_scores(
            collection_pairs, score_column, seed=permutation_seed
        )
        for rule, threshold in settings:
            graph = build_similarity_graph(
                null_pairs,
                identifiers,
                method=arguments.method,
                score_column=score_column,
                threshold=threshold,
                threshold_rule=rule,
                node_metadata=collection_metadata,
                node_id_column="node_id",
            )
            community_seed = stable_seed(
                arguments.seed,
                collection.collection_id,
                arguments.method,
                rule,
                threshold,
                "louvain",
            )
            partition = louvain_partition(
                graph,
                seed=community_seed,
                resolution=arguments.community_resolution,
                weight="score",
            )
            agreement = partition_agreement(partition, collection_metadata)
            summary = graph_summary(graph)
            for result in agreement.to_dict("records"):
                rows.append(
                    {
                        "graph_id": (
                            f"{collection.collection_id}__{arguments.method}__"
                            f"{rule}_{threshold}__score_permutation"
                        ),
                        **collection._asdict(),
                        "null_model": "within_collection_score_permutation",
                        "null_base_seed": arguments.seed,
                        "permutation_seed": permutation_seed,
                        "community_seed": community_seed,
                        **summary,
                        **result,
                        "community_method": "louvain",
                        "community_resolution": arguments.community_resolution,
                        "community_weight": "score",
                    }
                )

    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    results_path = output_directory / "graph_recovery_results.tsv"
    results = pd.DataFrame(rows).sort_values(
        ["graph_id", "reference_label"], ignore_index=True
    )
    results.to_csv(results_path, sep="\t", index=False)
    run_manifest = {
        "experiment": "within-collection score-permutation negative control",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "null_model": (
            "shuffle finite method scores among pair identities independently for "
            "each collection; preserve endpoints, score distribution, and missingness"
        ),
        "base_seed": arguments.seed,
        "method": arguments.method,
        "score_column": score_column,
        "composition": {
            "clans_per_collection": arguments.clans,
            "families_per_clan": arguments.families_per_clan,
            "domains_per_family": arguments.domains_per_family,
        },
        "collection_count": int(results["collection_id"].nunique()),
        "graph_count": int(results["graph_id"].nunique()),
        "graph_settings": [
            {"threshold_rule": rule, "threshold": threshold}
            for rule, threshold in settings
        ],
        "community_method": "NetworkX Louvain",
        "community_resolution": arguments.community_resolution,
        "community_weight": "score",
        "evaluation": {
            "primary": "adjusted mutual information (AMI)",
            "secondary": "normalized mutual information (NMI)",
            "average_method": "arithmetic",
            "reference_labels": ["family_accession", "clan_accession"],
        },
        "inputs": {
            str(path): sha256_file(path)
            for path in (
                arguments.collection_manifest,
                arguments.collection_membership,
                arguments.node_metadata,
                arguments.pair_scores,
            )
        },
        "results_sha256": sha256_file(results_path),
        "generator": "scripts/run_score_permutation_null.py",
    }
    (output_directory / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2) + "\n"
    )
    print(
        f"Evaluated {run_manifest['graph_count']} null graphs from "
        f"{run_manifest['collection_count']} collections"
    )
    print(f"Results: {results_path}")


if __name__ == "__main__":
    main()
