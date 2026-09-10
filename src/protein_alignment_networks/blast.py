"""Local NCBI BLAST+ integration for protein pair comparisons."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from .io import validate_protein_sequence, write_fasta

BLAST_COLUMNS = [
    "query_id",
    "subject_id",
    "raw_score",
    "bit_score",
    "evalue",
    "alignment_length",
    "percent_identity",
    "identical_residues",
    "mismatches",
    "gap_openings",
    "gaps",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "query_length",
    "subject_length",
    "query_coverage_hsp",
]

_INTEGER_COLUMNS = [
    "raw_score",
    "alignment_length",
    "identical_residues",
    "mismatches",
    "gap_openings",
    "gaps",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "query_length",
    "subject_length",
]

_OUTPUT_FORMAT = (
    "6 qseqid sseqid score bitscore evalue length pident nident mismatch "
    "gapopen gaps qstart qend sstart send qlen slen qcovhsp"
)


def require_blast_executable(name: str) -> str:
    """Return an executable path or explain how to install BLAST+."""

    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(
            f"NCBI BLAST+ executable '{name}' was not found. "
            "Install it on macOS with 'brew install blast'."
        )
    return executable


def blastp_version(blastp: str | None = None) -> str:
    """Return the first line of ``blastp -version``."""

    executable = blastp or require_blast_executable("blastp")
    completed = _run([executable, "-version"])
    return completed.stdout.splitlines()[0].strip()


def parse_blast_tabular(
    path: str | Path,
    *,
    identifiers: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Read this project's explicit BLAST tabular output format."""

    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return _empty_blast_frame()

    hits = pd.read_csv(path, sep="\t", names=BLAST_COLUMNS, header=None)
    for column in _INTEGER_COLUMNS:
        hits[column] = pd.to_numeric(hits[column], errors="raise").astype("int64")
    for column in ["bit_score", "evalue", "percent_identity", "query_coverage_hsp"]:
        hits[column] = pd.to_numeric(hits[column], errors="raise").astype("float64")

    if identifiers is not None:
        lookup = {f"sequence_{index:06d}": identifier for index, identifier in enumerate(identifiers)}
        for column in ["query_id", "subject_id"]:
            unknown = sorted(set(hits[column]) - set(lookup))
            if unknown:
                raise ValueError(f"unknown internal BLAST identifier(s): {', '.join(unknown)}")
            hits[column] = hits[column].map(lookup)
    return hits


def summarize_blast_pairs(hits: pd.DataFrame) -> pd.DataFrame:
    """Select the strongest HSP for each unique, non-self protein pair.

    BLAST searches are directional and may report several high-scoring segment
    pairs (HSPs). For an undirected graph, this function canonicalises each
    pair and retains the row with the highest bit score, using the lowest
    E-value and longest alignment as deterministic tie-breakers.
    """

    missing = set(BLAST_COLUMNS) - set(hits.columns)
    if missing:
        raise ValueError(f"missing BLAST column(s): {', '.join(sorted(missing))}")
    pairs = hits.loc[hits["query_id"] != hits["subject_id"]].copy()
    if pairs.empty:
        return pd.DataFrame(columns=["protein_a", "protein_b", *BLAST_COLUMNS])

    endpoints = pairs.apply(
        lambda row: sorted((str(row["query_id"]), str(row["subject_id"]))),
        axis=1,
        result_type="expand",
    )
    pairs[["protein_a", "protein_b"]] = endpoints
    pairs = pairs.sort_values(
        ["protein_a", "protein_b", "bit_score", "evalue", "alignment_length"],
        ascending=[True, True, False, True, False],
        kind="stable",
    )
    pairs = pairs.drop_duplicates(["protein_a", "protein_b"], keep="first")
    return pairs[["protein_a", "protein_b", *BLAST_COLUMNS]].reset_index(drop=True)


def summarize_blast_directionality(hits: pd.DataFrame) -> dict[str, int | float | None]:
    """Audit whether an all-vs-all BLAST pair was reported in both directions.

    BLAST searches are directional even though the project ultimately builds
    undirected graphs.  The relative score spread for bidirectional pairs is
    ``abs(forward - reverse) / max(forward, reverse)`` after retaining the
    strongest bit score in each ordered direction.
    """

    required = {"query_id", "subject_id", "bit_score", "evalue", "alignment_length"}
    missing = required - set(hits.columns)
    if missing:
        raise ValueError(f"missing BLAST column(s): {', '.join(sorted(missing))}")
    nonself = hits.loc[hits["query_id"] != hits["subject_id"]].copy()
    if nonself.empty:
        return {
            "hsp_row_count": 0,
            "undirected_pair_count": 0,
            "bidirectional_pair_count": 0,
            "single_direction_pair_count": 0,
            "bidirectional_relative_spread_median": None,
            "bidirectional_relative_spread_p90": None,
            "bidirectional_relative_spread_max": None,
        }
    ordered = (
        nonself.sort_values(
            ["query_id", "subject_id", "bit_score", "evalue", "alignment_length"],
            ascending=[True, True, False, True, False],
            kind="stable",
        )
        .drop_duplicates(["query_id", "subject_id"], keep="first")
        .copy()
    )
    endpoints = ordered.apply(
        lambda row: sorted((str(row["query_id"]), str(row["subject_id"]))),
        axis=1,
        result_type="expand",
    )
    ordered[["protein_a", "protein_b"]] = endpoints
    grouped = ordered.groupby(["protein_a", "protein_b"])["bit_score"]
    direction_counts = grouped.size()
    bidirectional_scores = grouped.apply(list).loc[direction_counts == 2]
    spreads = bidirectional_scores.map(
        lambda values: abs(values[0] - values[1]) / max(values)
    )
    return {
        "hsp_row_count": len(hits),
        "undirected_pair_count": len(direction_counts),
        "bidirectional_pair_count": int((direction_counts == 2).sum()),
        "single_direction_pair_count": int((direction_counts == 1).sum()),
        "bidirectional_relative_spread_median": (
            float(spreads.median()) if not spreads.empty else None
        ),
        "bidirectional_relative_spread_p90": (
            float(spreads.quantile(0.9)) if not spreads.empty else None
        ),
        "bidirectional_relative_spread_max": (
            float(spreads.max()) if not spreads.empty else None
        ),
    }


def blastp_all_vs_all(
    sequences: Mapping[str, str],
    *,
    output_path: str | Path | None = None,
    threads: int = 1,
    evalue: float = 10.0,
    blastp: str | None = None,
    makeblastdb: str | None = None,
) -> pd.DataFrame:
    """Run standard local BLASTP for every query against the same sequence set.

    The explicit scoring settings match standard BLASTP defaults: BLOSUM62,
    gap-open cost 11, gap-extension cost 1, and composition-based statistics.
    The returned table contains all reported HSPs. An absent pair means BLAST
    reported no hit at the chosen E-value threshold; it is not a zero score.
    """

    if not sequences:
        raise ValueError("at least one sequence is required")
    if threads < 1:
        raise ValueError("threads must be positive")
    if evalue <= 0:
        raise ValueError("evalue must be positive")

    identifiers = list(sequences)
    clean_sequences = {
        f"sequence_{index:06d}": validate_protein_sequence(sequence, allow_empty=False)
        for index, sequence in enumerate(sequences.values())
    }
    blastp_executable = blastp or require_blast_executable("blastp")
    makeblastdb_executable = makeblastdb or require_blast_executable("makeblastdb")

    with tempfile.TemporaryDirectory(prefix="protein_alignment_blast_") as temporary:
        temporary_path = Path(temporary)
        fasta_path = temporary_path / "sequences.fasta"
        database_path = temporary_path / "database"
        temporary_output = temporary_path / "hits.tsv"
        write_fasta(clean_sequences, fasta_path)

        _run(
            [
                makeblastdb_executable,
                "-in",
                str(fasta_path),
                "-dbtype",
                "prot",
                "-parse_seqids",
                "-out",
                str(database_path),
            ]
        )
        _run(
            [
                blastp_executable,
                "-query",
                str(fasta_path),
                "-db",
                str(database_path),
                "-matrix",
                "BLOSUM62",
                "-gapopen",
                "11",
                "-gapextend",
                "1",
                "-comp_based_stats",
                "2",
                "-evalue",
                str(evalue),
                "-max_target_seqs",
                str(len(clean_sequences)),
                "-num_threads",
                str(threads),
                "-outfmt",
                _OUTPUT_FORMAT,
                "-out",
                str(temporary_output),
            ]
        )
        hits = parse_blast_tabular(temporary_output, identifiers=identifiers)
        if output_path is not None:
            destination = Path(output_path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            hits.to_csv(destination, sep="\t", index=False)
        return hits


def _empty_blast_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=BLAST_COLUMNS)


def _run(command: list[str]) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(f"executable not found: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip() or "no diagnostic output"
        raise RuntimeError(f"BLAST+ command failed: {detail}") from error
