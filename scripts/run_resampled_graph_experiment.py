"""Build and evaluate graph series for reproducibly sampled Pfam collections."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

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
    "dedal": "dedal_sw_score",
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def graph_seed(sampling_seed: int, method: str, rule: str, threshold: float) -> int:
    payload = f"{sampling_seed}|{method}|{rule}|{threshold}".encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("collection_manifest", type=Path)
    parser.add_argument("collection_membership", type=Path)
    parser.add_argument("node_metadata", type=Path)
    parser.add_argument("pair_scores", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=sorted(METHOD_SCORE_COLUMNS),
        default=["biopython"],
    )
    parser.add_argument("--top-k", nargs="*", type=int, default=[1, 2, 5, 10])
    parser.add_argument(
        "--percentiles", nargs="*", type=float, default=[]
    )
    parser.add_argument(
        "--target-densities", nargs="*", type=float, default=[]
    )
    parser.add_argument("--community-resolution", type=float, default=1.0)
    parser.add_argument(
        "--unweighted-communities",
        action="store_true",
        help="ignore retained edge scores during Louvain community detection",
    )
    parser.add_argument("--clans", nargs="*", type=int)
    parser.add_argument("--families-per-clan", nargs="*", type=int)
    parser.add_argument("--domains-per-family", nargs="*", type=int)
    parser.add_argument("--replicates", nargs="*", type=int)
    parser.add_argument(
        "--max-collections",
        type=int,
        help="deterministic smoke-test limit after all other filters",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=25,
        help="write a recoverable partial results table after this many collections",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse complete graph IDs already present in the output table",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if arguments.community_resolution <= 0:
        raise ValueError("--community-resolution must be positive")
    if any(value < 0 for value in arguments.top_k):
        raise ValueError("--top-k values must be non-negative")
    if any(not 0 <= value <= 1 for value in arguments.percentiles):
        raise ValueError("--percentiles values must be between zero and one")
    if any(not 0 <= value <= 1 for value in arguments.target_densities):
        raise ValueError("--target-densities values must be between zero and one")
    if (
        not arguments.top_k
        and not arguments.percentiles
        and not arguments.target_densities
    ):
        raise ValueError("at least one graph rule/threshold is required")
    if arguments.checkpoint_every < 1:
        raise ValueError("--checkpoint-every must be positive")

    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    results_path = output_directory / "graph_recovery_results.tsv"
    run_manifest_path = output_directory / "run_manifest.json"

    manifest = pd.read_csv(arguments.collection_manifest, sep="\t")
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
    for column, values in (
        ("clans_per_collection", arguments.clans),
        ("families_per_clan", arguments.families_per_clan),
        ("domains_per_family", arguments.domains_per_family),
        ("replicate", arguments.replicates),
    ):
        if values:
            manifest = manifest.loc[manifest[column].isin(values)]
    manifest = manifest.sort_values("collection_id", ignore_index=True)
    if arguments.max_collections is not None:
        if arguments.max_collections < 1:
            raise ValueError("--max-collections must be positive")
        manifest = manifest.head(arguments.max_collections)
    if manifest.empty:
        raise ValueError("no collections remain after filtering")

    missing_scores = [
        METHOD_SCORE_COLUMNS[method]
        for method in arguments.methods
        if METHOD_SCORE_COLUMNS[method] not in pair_scores.columns
    ]
    if missing_scores:
        raise ValueError(f"pair table is missing scores: {', '.join(missing_scores)}")

    graph_settings = [
        ("top_k", value) for value in arguments.top_k
    ] + [
        ("percentile", value) for value in arguments.percentiles
    ] + [
        ("target_density", value) for value in arguments.target_densities
    ]
    metadata_index = metadata.set_index("node_id", drop=False)
    members_by_collection = (
        membership.groupby("collection_id", sort=False)["node_id"].apply(list).to_dict()
    )
    if arguments.resume and results_path.exists():
        previous_results = pd.read_csv(results_path, sep="\t")
        rows = previous_results.to_dict("records")
        completed_counts = previous_results.groupby("graph_id")[
            "reference_label"
        ].nunique()
        completed_graph_ids = set(completed_counts.loc[completed_counts == 2].index)
        print(f"Resuming with {len(completed_graph_ids)} complete graphs")
    else:
        rows: list[dict] = []
        completed_graph_ids: set[str] = set()

    def checkpoint() -> None:
        results = pd.DataFrame(rows).sort_values(
            ["graph_id", "reference_label"], ignore_index=True
        )
        results.to_csv(results_path, sep="\t", index=False)

    collection_count = len(manifest)
    for collection_index, collection in enumerate(
        manifest.itertuples(index=False), start=1
    ):
        identifiers = members_by_collection.get(collection.collection_id, [])
        identifier_set = set(identifiers)
        if len(identifier_set) != collection.node_count:
            raise ValueError(
                f"membership count does not match {collection.collection_id} manifest"
            )
        collection_metadata = metadata_index.loc[identifiers].reset_index(drop=True)
        collection_pairs = pair_scores.loc[
            pair_scores["protein_a"].isin(identifier_set)
            & pair_scores["protein_b"].isin(identifier_set)
        ]
        for method in arguments.methods:
            score_column = METHOD_SCORE_COLUMNS[method]
            for rule, threshold in graph_settings:
                graph_id = (
                    f"{collection.collection_id}__{method}__{rule}_{threshold}"
                )
                if graph_id in completed_graph_ids:
                    continue
                graph = build_similarity_graph(
                    collection_pairs,
                    identifiers,
                    method=method,
                    score_column=score_column,
                    threshold=threshold,
                    threshold_rule=rule,
                    node_metadata=collection_metadata,
                    node_id_column="node_id",
                )
                partition = louvain_partition(
                    graph,
                    seed=graph_seed(
                        collection.sampling_seed, method, rule, threshold
                    ),
                    resolution=arguments.community_resolution,
                    weight=None if arguments.unweighted_communities else "score",
                )
                agreement = partition_agreement(partition, collection_metadata)
                summary = graph_summary(graph)
                for result in agreement.to_dict("records"):
                    rows.append(
                        {
                            "graph_id": graph_id,
                            **collection._asdict(),
                            **summary,
                            **result,
                            "community_method": "louvain",
                            "community_resolution": arguments.community_resolution,
                            "community_weight": (
                                "none"
                                if arguments.unweighted_communities
                                else "score"
                            ),
                        }
                    )
                completed_graph_ids.add(graph_id)
        if (
            collection_index % arguments.checkpoint_every == 0
            or collection_index == collection_count
        ):
            checkpoint()
            print(
                f"Processed {collection_index}/{collection_count} collections; "
                f"{len(completed_graph_ids)} graphs complete",
                flush=True,
            )

    results = pd.read_csv(results_path, sep="\t")
    run_manifest = {
        "experiment": "repeated Pfam-domain graph label recovery",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "collection_count": int(manifest["collection_id"].nunique()),
        "graph_count": int(results["graph_id"].nunique()),
        "result_row_count": len(results),
        "methods": arguments.methods,
        "graph_settings": [
            {"threshold_rule": rule, "threshold": threshold}
            for rule, threshold in graph_settings
        ],
        "community_method": "NetworkX Louvain",
        "community_resolution": arguments.community_resolution,
        "community_weight": (
            "none" if arguments.unweighted_communities else "score"
        ),
        "evaluation": {
            "primary": "normalized mutual information (NMI)",
            "secondary": "adjusted mutual information (AMI)",
            "reference_labels": ["family_accession", "clan_accession"],
            "labels_used_during_graph_or_community_construction": False,
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
        "generator": "scripts/run_resampled_graph_experiment.py",
    }
    run_manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n")
    print(
        f"Evaluated {run_manifest['graph_count']} graphs from "
        f"{run_manifest['collection_count']} sampled collections"
    )
    print(f"Results: {results_path}")
    print(f"Manifest: {run_manifest_path}")


if __name__ == "__main__":
    main()
