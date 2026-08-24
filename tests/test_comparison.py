import pandas as pd
import pytest

from protein_alignment_networks.comparison import (
    add_descending_ranks,
    canonical_pair_index,
    merge_pair_tables,
    pair_rank_disagreements,
    score_matrix_to_pairs,
    spearman_correlations,
)


def test_canonical_pair_index_is_complete_and_sorted():
    pairs = canonical_pair_index(["c", "a", "b"])
    assert pairs.to_dict("records") == [
        {"protein_a": "a", "protein_b": "b"},
        {"protein_a": "a", "protein_b": "c"},
        {"protein_a": "b", "protein_b": "c"},
    ]


def test_score_matrix_to_pairs_uses_labels_and_excludes_diagonal():
    matrix = pd.DataFrame(
        [[10.0, 3.0, 4.0], [3.0, 11.0, 5.0], [4.0, 5.0, 12.0]],
        index=["b", "a", "c"],
        columns=["b", "a", "c"],
    )
    pairs = score_matrix_to_pairs(matrix, score_column="score")
    assert pairs.to_dict("records") == [
        {"protein_a": "a", "protein_b": "b", "score": 3.0},
        {"protein_a": "a", "protein_b": "c", "score": 5.0},
        {"protein_a": "b", "protein_b": "c", "score": 4.0},
    ]


def test_merge_pair_tables_canonicalises_direction_and_preserves_missing_pairs():
    universe = canonical_pair_index(["a", "b", "c"])
    method = pd.DataFrame(
        {
            "protein_a": ["b", "c"],
            "protein_b": ["a", "a"],
            "method_score": [7.0, 8.0],
        }
    )
    merged = merge_pair_tables(universe, method)
    assert merged["method_score"].tolist()[:2] == [7.0, 8.0]
    assert pd.isna(merged.loc[2, "method_score"])


def test_ranks_and_spearman_correlation_handle_different_scales_and_missing():
    scores = pd.DataFrame(
        {
            "a": [1.0, 2.0, 3.0, 4.0],
            "b": [10.0, 20.0, 30.0, 40.0],
            "c": [4.0, 3.0, 2.0, None],
        }
    )
    ranked = add_descending_ranks(scores, ["a", "b"])
    assert ranked["a_rank"].tolist() == [4.0, 3.0, 2.0, 1.0]
    assert ranked["a_percentile"].tolist() == [0.25, 0.5, 0.75, 1.0]

    correlations = spearman_correlations(scores, ["a", "b", "c"])
    ab = correlations.query("score_a == 'a' and score_b == 'b'").iloc[0]
    ac = correlations.query("score_a == 'a' and score_b == 'c'").iloc[0]
    assert ab["spearman_rho"] == pytest.approx(1.0)
    assert ab["n_pairs"] == 4
    assert ac["spearman_rho"] == pytest.approx(-1.0)
    assert ac["n_pairs"] == 3


def test_pair_rank_disagreements_uses_the_same_complete_pairs_for_both_scores():
    scores = pd.DataFrame(
        {
            "protein_a": ["a", "a", "b", "b"],
            "protein_b": ["b", "c", "c", "d"],
            "score_a": [1.0, 2.0, 3.0, 100.0],
            "score_b": [10.0, 30.0, 20.0, None],
        }
    )
    disagreements = pair_rank_disagreements(scores, "score_a", "score_b")
    assert len(disagreements) == 3
    assert disagreements.loc[0, ["protein_a", "protein_b"]].tolist() == ["a", "c"]
    assert disagreements.loc[0, "absolute_percentile_difference"] == pytest.approx(
        1 / 3
    )
