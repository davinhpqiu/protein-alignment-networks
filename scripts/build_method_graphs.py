"""Build separate method-specific graphs from a saved protein-pair score table."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from protein_alignment_networks import (
    assess_score_correlations,
    build_similarity_graph,
    graph_summary,
    read_fasta,
    save_graph_bundle,
)

METHODS = ("biopython", "blast", "dedal")
RULES = ("absolute", "percentile", "top_n")
INDEPENDENT_SCORE_COLUMNS = {
    "biopython_score",
    "blast_bit_score",
    "dedal_sw_score",
}


def parse_graph_spec(value: str) -> tuple[str, str, float]:
    """Parse ``METHOD:RULE:THRESHOLD`` from the command line."""

    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "graph specification must be METHOD:RULE:THRESHOLD"
        )
    method, rule, raw_threshold = parts
    if method not in METHODS:
        raise argparse.ArgumentTypeError(f"unknown method: {method}")
    if rule not in RULES:
        raise argparse.ArgumentTypeError(f"unknown threshold rule: {rule}")
    try:
        threshold = float(raw_threshold)
    except ValueError as error:
        raise argparse.ArgumentTypeError("graph threshold must be numeric") from error
    if rule == "top_n":
        if not threshold.is_integer():
            raise argparse.ArgumentTypeError("top_n threshold must be an integer")
        threshold = int(threshold)
    return method, rule, threshold


def sha256_file(path: Path) -> str:
    """Return the SHA-256 checksum of a file."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--node-id-column", default="uniprot_accession")
    parser.add_argument("--correlations", type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument(
        "--graph",
        action="append",
        required=True,
        type=parse_graph_spec,
        help="repeat METHOD:RULE:THRESHOLD for each separate graph",
    )
    parser.add_argument("--minimum-rho", type=float)
    parser.add_argument("--minimum-overlap-fraction", type=float)
    parser.add_argument("--minimum-pairs", type=int)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    sequences = read_fasta(arguments.fasta)
    pairs = pd.read_csv(arguments.pairs, sep="\t")
    metadata = (
        pd.read_csv(arguments.metadata, sep="\t") if arguments.metadata else None
    )
    output = arguments.output_directory
    output.mkdir(parents=True, exist_ok=True)

    quality_path: Path | None = None
    quality_rows = 0
    if arguments.correlations:
        correlations = pd.read_csv(arguments.correlations, sep="\t")
        correlations = correlations.loc[
            correlations["score_a"].isin(INDEPENDENT_SCORE_COLUMNS)
            & correlations["score_b"].isin(INDEPENDENT_SCORE_COLUMNS)
        ].reset_index(drop=True)
        quality = assess_score_correlations(
            correlations,
            total_pairs=len(pairs),
            minimum_rho=arguments.minimum_rho,
            minimum_overlap_fraction=arguments.minimum_overlap_fraction,
            minimum_pairs=arguments.minimum_pairs,
        )
        quality_path = output / "quality_flags.tsv"
        quality.to_csv(quality_path, sep="\t", index=False)
        quality_rows = int((quality["status"] == "flag").sum())

    summaries = []
    generated_files: list[Path] = []
    graph_specs = []
    for method, rule, threshold in arguments.graph:
        graph = build_similarity_graph(
            pairs,
            sequences,
            method=method,
            threshold=threshold,
            threshold_rule=rule,
            node_metadata=metadata,
            node_id_column=arguments.node_id_column if metadata is not None else None,
        )
        label = str(threshold).replace(".", "p")
        stem = output / f"{method}_{rule}_{label}"
        bundle = save_graph_bundle(graph, stem)
        generated_files.extend(bundle.values())
        summaries.append(graph_summary(graph))
        graph_specs.append(
            {"method": method, "threshold_rule": rule, "threshold": threshold}
        )

    summaries_path = output / "graph_summaries.tsv"
    pd.DataFrame(summaries).to_csv(summaries_path, sep="\t", index=False)
    generated_files.append(summaries_path)
    if quality_path is not None:
        generated_files.append(quality_path)

    manifest = {
        "experiment": "protein set to method-specific graphs",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "fasta": str(arguments.fasta),
        "fasta_sha256": sha256_file(arguments.fasta),
        "pair_scores": str(arguments.pairs),
        "pair_scores_sha256": sha256_file(arguments.pairs),
        "protein_count": len(sequences),
        "pair_count": len(pairs),
        "graphs": graph_specs,
        "quality_thresholds": {
            "minimum_rho": arguments.minimum_rho,
            "minimum_overlap_fraction": arguments.minimum_overlap_fraction,
            "minimum_pairs": arguments.minimum_pairs,
        },
        "quality_flag_count": quality_rows,
        "files": {str(path): sha256_file(path) for path in generated_files},
    }
    manifest_path = output / "graph_run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Built {len(graph_specs)} separate graph(s) from {len(sequences)} proteins")
    print(f"Quality thresholds configured: {any(value is not None for value in manifest['quality_thresholds'].values())}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
