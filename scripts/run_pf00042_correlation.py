"""Run pair-level alignment scoring for a prepared protein collection."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from datetime import UTC, datetime
from pathlib import Path

import Bio
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

from protein_alignment_networks import (
    add_descending_ranks,
    blastp_all_vs_all,
    blastp_version,
    canonical_pair_index,
    dedal_all_vs_all,
    merge_pair_tables,
    pairwise_score_matrix,
    read_fasta,
    score_matrix_to_pairs,
    spearman_correlations,
    summarize_blast_pairs,
)

METHODS = ("biopython", "blast", "dedal")
PRIMARY_SCORE_COLUMNS = (
    "biopython_score",
    "blast_bit_score",
    "dedal_sw_score",
    "dedal_homology_logit",
)


def sha256_file(path: Path) -> str:
    """Return the SHA-256 checksum of a file."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def recorded_path(path: Path, project_root: Path) -> str:
    """Prefer a portable project-relative path in the run manifest."""

    try:
        return str(path.resolve().relative_to(project_root.resolve()))
    except ValueError:
        return str(path.resolve())


def configure_biopython_aligner() -> Align.PairwiseAligner:
    """Return the project's recorded local protein-alignment baseline."""

    return Align.PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,
        extend_gap_score=-1,
    )


def run_biopython(
    sequences: dict[str, str],
    output_directory: Path,
) -> tuple[pd.DataFrame, float]:
    """Run and save the Biopython score matrix and unique-pair table."""

    started = time.perf_counter()
    matrix = pairwise_score_matrix(
        sequences,
        aligner=configure_biopython_aligner(),
    )
    pairs = score_matrix_to_pairs(matrix, score_column="biopython_score")
    elapsed = time.perf_counter() - started
    matrix.to_csv(output_directory / "biopython_score_matrix.tsv", sep="\t")
    pairs.to_csv(output_directory / "biopython_pairs.tsv", sep="\t", index=False)
    return pairs, elapsed


def run_blast(
    sequences: dict[str, str],
    output_directory: Path,
    *,
    threads: int,
) -> tuple[pd.DataFrame, float]:
    """Run and save BLAST HSPs plus the strongest result for each pair."""

    started = time.perf_counter()
    hits = blastp_all_vs_all(
        sequences,
        output_path=output_directory / "blast_hsps.tsv",
        threads=threads,
        evalue=10.0,
    )
    pairs = summarize_blast_pairs(hits)
    elapsed = time.perf_counter() - started
    pairs.to_csv(output_directory / "blast_pairs.tsv", sep="\t", index=False)
    return pairs, elapsed


def run_dedal(
    sequences: dict[str, str],
    output_directory: Path,
) -> tuple[pd.DataFrame, float]:
    """Run and save pretrained DEDAL results for every unique pair."""

    started = time.perf_counter()
    pairs = dedal_all_vs_all(
        sequences,
        output_path=output_directory / "dedal_pairs.tsv",
    )
    return pairs, time.perf_counter() - started


def blast_scores(pairs: pd.DataFrame) -> pd.DataFrame:
    """Select and prefix BLAST fields used in the combined pair table."""

    columns = [
        "protein_a",
        "protein_b",
        "raw_score",
        "bit_score",
        "evalue",
        "percent_identity",
        "alignment_length",
        "query_coverage_hsp",
    ]
    return pairs[columns].rename(
        columns={column: f"blast_{column}" for column in columns[2:]}
    )


def dedal_scores(pairs: pd.DataFrame) -> pd.DataFrame:
    """Select and prefix DEDAL fields used in the combined pair table."""

    columns = [
        "protein_a",
        "protein_b",
        "sw_score",
        "homology_logit",
        "homology_probability",
        "identity_count",
        "similarity_count",
        "gap_count",
        "alignment_length",
    ]
    return pairs[columns].rename(
        columns={column: f"dedal_{column}" for column in columns[2:]}
    )


def load_available_tables(output_directory: Path) -> list[pd.DataFrame]:
    """Load every completed method table for the merged analysis."""

    tables: list[pd.DataFrame] = []
    biopython_path = output_directory / "biopython_pairs.tsv"
    blast_path = output_directory / "blast_pairs.tsv"
    dedal_path = output_directory / "dedal_pairs.tsv"
    if biopython_path.exists():
        tables.append(pd.read_csv(biopython_path, sep="\t"))
    if blast_path.exists():
        tables.append(blast_scores(pd.read_csv(blast_path, sep="\t")))
    if dedal_path.exists():
        tables.append(dedal_scores(pd.read_csv(dedal_path, sep="\t")))
    return tables


def write_combined_outputs(
    identifiers: list[str],
    output_directory: Path,
) -> tuple[pd.DataFrame, list[str]]:
    """Merge completed methods, rank scores, and save correlations."""

    pair_index = canonical_pair_index(identifiers)
    combined = merge_pair_tables(pair_index, *load_available_tables(output_directory))
    score_columns = [
        column for column in PRIMARY_SCORE_COLUMNS if column in combined.columns
    ]
    combined = (
        add_descending_ranks(combined, score_columns)
        if len(score_columns) >= 2
        else combined
    )
    combined.to_csv(output_directory / "paired_scores.tsv", sep="\t", index=False)
    if len(score_columns) >= 2:
        correlations = spearman_correlations(combined, score_columns)
    else:
        correlations = pd.DataFrame(
            columns=["score_a", "score_b", "spearman_rho", "n_pairs"]
        )
    correlations.to_csv(
        output_directory / "spearman_correlations.tsv", sep="\t", index=False
    )
    return combined, score_columns


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--fasta",
        type=Path,
        default=root
        / "data/processed/pfam_pf00042_globin_pilot/PF00042_globin_pilot.fasta",
    )
    parser.add_argument(
        "--dataset-manifest",
        type=Path,
        default=root
        / "data/processed/pfam_pf00042_globin_pilot/"
        "PF00042_globin_pilot_manifest.json",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=root / "outputs/tables/pfam_pf00042_correlation",
    )
    parser.add_argument(
        "--experiment-name",
        help="descriptive name recorded in the run manifest",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=METHODS,
        default=list(METHODS),
        help="methods to run; existing outputs from other methods are still merged",
    )
    parser.add_argument("--threads", type=int, default=1, help="BLAST worker threads")
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace outputs for the requested methods",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if arguments.threads < 1:
        raise ValueError("threads must be positive")
    sequences = read_fasta(arguments.fasta)
    project_root = Path(__file__).resolve().parents[1]
    output_directory = arguments.output_directory
    output_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = output_directory / "run_manifest.json"
    previous_manifest = (
        json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    )
    runtimes: dict[str, float] = {}

    runners = {
        "biopython": (
            output_directory / "biopython_pairs.tsv",
            lambda: run_biopython(sequences, output_directory),
        ),
        "blast": (
            output_directory / "blast_pairs.tsv",
            lambda: run_blast(sequences, output_directory, threads=arguments.threads),
        ),
        "dedal": (
            output_directory / "dedal_pairs.tsv",
            lambda: run_dedal(sequences, output_directory),
        ),
    }
    for method in arguments.methods:
        destination, runner = runners[method]
        if destination.exists() and not arguments.force:
            print(f"Reusing {method}: {destination}")
            continue
        print(f"Running {method}...")
        _result, runtimes[method] = runner()
        print(f"Completed {method} in {runtimes[method]:.3f} seconds")

    combined, score_columns = write_combined_outputs(
        list(sequences), output_directory
    )
    completed_methods = [
        method for method, (path, _runner) in runners.items() if path.exists()
    ]
    output_files = sorted(
        path for path in output_directory.iterdir() if path.name != "run_manifest.json"
    )
    cumulative_runtimes = dict(
        previous_manifest.get(
            "runtimes_seconds",
            previous_manifest.get("runtimes_seconds_this_run", {}),
        )
    )
    cumulative_runtimes.update(runtimes)
    dataset_record = json.loads(arguments.dataset_manifest.read_text())
    dataset_name = dataset_record.get("dataset", arguments.dataset_manifest.stem)
    manifest = {
        "experiment": arguments.experiment_name
        or f"{dataset_name} pair-level alignment scores",
        "updated_at_utc": datetime.now(UTC).isoformat(),
        "dataset_manifest": recorded_path(arguments.dataset_manifest, project_root),
        "dataset_manifest_sha256": sha256_file(arguments.dataset_manifest),
        "fasta": recorded_path(arguments.fasta, project_root),
        "fasta_sha256": sha256_file(arguments.fasta),
        "sequence_count": len(sequences),
        "pair_count": len(combined),
        "methods_requested_this_run": arguments.methods,
        "methods_completed": completed_methods,
        "score_columns": score_columns,
        "runtimes_seconds_this_run": runtimes,
        "runtimes_seconds": cumulative_runtimes,
        "parameters": {
            "biopython": {
                "mode": "local",
                "substitution_matrix": "BLOSUM62",
                "open_gap_score": -11,
                "extend_gap_score": -1,
            },
            "blast": {
                "program": "blastp",
                "matrix": "BLOSUM62",
                "gap_open": 11,
                "gap_extend": 1,
                "composition_based_statistics": 2,
                "evalue_reporting_threshold": 10.0,
                "pair_summary": "highest bit score HSP",
            },
            "dedal": {
                "model": "https://tfhub.dev/google/dedal/3",
                "primary_alignment_score": "sw_score",
                "additional_classifier_output": "homology_logit",
            },
            "correlation": {
                "statistic": "Spearman rank correlation",
                "missing_values": "pairwise deletion; never replace absent hits with zero",
                "ties": "average ranks",
            },
        },
        "software": {
            "python": platform.python_version(),
            "biopython": Bio.__version__,
            "pandas": pd.__version__,
            "blastp": blastp_version() if "blast" in completed_methods else None,
        },
        "files": {
            recorded_path(path, project_root): sha256_file(path)
            for path in output_files
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Merged {len(combined)} canonical pairs")
    print(f"Completed methods: {', '.join(completed_methods) or 'none'}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
