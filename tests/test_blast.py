import shutil
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import protein_alignment_networks.blast as blast_module
from protein_alignment_networks.blast import (
    BLAST_COLUMNS,
    blastp_all_vs_all,
    parse_blast_tabular,
    summarize_blast_directionality,
    summarize_blast_pairs,
)


def test_parse_and_summarize_blast_hits(tmp_path: Path):
    output = tmp_path / "hits.tsv"
    output.write_text(
        "sequence_000000\tsequence_000001\t20\t15.0\t0.01\t3\t66.667\t2\t1\t0\t0\t1\t3\t1\t3\t3\t3\t100\n"
        "sequence_000001\tsequence_000000\t22\t17.0\t0.005\t3\t66.667\t2\t1\t0\t0\t1\t3\t1\t3\t3\t3\t100\n",
        encoding="utf-8",
    )
    hits = parse_blast_tabular(output, identifiers=["alpha", "beta"])
    assert list(hits.columns) == BLAST_COLUMNS
    assert set(hits["query_id"]) == {"alpha", "beta"}

    pairs = summarize_blast_pairs(hits)
    assert len(pairs) == 1
    assert pairs.loc[0, "protein_a"] == "alpha"
    assert pairs.loc[0, "protein_b"] == "beta"
    assert pairs.loc[0, "bit_score"] == 17.0


def test_summary_removes_self_hits():
    hits = pd.DataFrame(
        [["a", "a", 10, 10.0, 0.1, 3, 100.0, 3, 0, 0, 0, 1, 3, 1, 3, 3, 3, 100.0]],
        columns=BLAST_COLUMNS,
    )
    assert summarize_blast_pairs(hits).empty


def test_directionality_summary_counts_pairs_and_defines_relative_spread():
    rows = [
        ["a", "b", 10, 20.0, 0.1, 3, 100.0, 3, 0, 0, 0, 1, 3, 1, 3, 3, 3, 100.0],
        ["b", "a", 10, 18.0, 0.1, 3, 100.0, 3, 0, 0, 0, 1, 3, 1, 3, 3, 3, 100.0],
        ["a", "c", 10, 12.0, 0.1, 3, 100.0, 3, 0, 0, 0, 1, 3, 1, 3, 3, 3, 100.0],
        ["a", "a", 10, 30.0, 0.1, 3, 100.0, 3, 0, 0, 0, 1, 3, 1, 3, 3, 3, 100.0],
    ]
    summary = summarize_blast_directionality(pd.DataFrame(rows, columns=BLAST_COLUMNS))

    assert summary["hsp_row_count"] == 4
    assert summary["undirected_pair_count"] == 2
    assert summary["bidirectional_pair_count"] == 1
    assert summary["single_direction_pair_count"] == 1
    assert summary["bidirectional_relative_spread_median"] == pytest.approx(0.1)


def test_all_vs_all_sets_target_cap_to_full_sequence_count(monkeypatch):
    commands = []

    def fake_run(command):
        commands.append(command)
        if "-outfmt" in command:
            Path(command[command.index("-out") + 1]).write_text("")
        return SimpleNamespace(stdout="")

    monkeypatch.setattr(blast_module, "_run", fake_run)
    hits = blastp_all_vs_all(
        {"a": "AAAA", "b": "AAAT", "c": "AATT"},
        blastp="blastp",
        makeblastdb="makeblastdb",
    )

    blast_command = commands[1]
    assert blast_command[blast_command.index("-max_target_seqs") + 1] == "3"
    assert hits.empty


@pytest.mark.skipif(
    shutil.which("blastp") is None or shutil.which("makeblastdb") is None,
    reason="NCBI BLAST+ is not installed",
)
def test_local_blast_all_vs_all_smoke_test():
    sequences = {
        "protein_a": "MKTAYIAKQRQISFVKSHFSRQ",
        "protein_b": "MKTAYIAKQRQISFVKSHFSRN",
        "protein_c": "GGGGGGGGGGGGGGGGGGGGGG",
    }
    hits = blastp_all_vs_all(sequences, evalue=1000)
    assert not hits.empty
    assert set(hits["query_id"]).issubset(sequences)
    assert set(hits["subject_id"]).issubset(sequences)
    assert {"bit_score", "evalue", "percent_identity"}.issubset(hits.columns)
    pairs = summarize_blast_pairs(hits)
    assert ((pairs["protein_a"] == "protein_a") & (pairs["protein_b"] == "protein_b")).any()
