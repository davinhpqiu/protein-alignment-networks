"""Utilities for comparing method scores on the same protein pairs."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from itertools import combinations

import numpy as np
import pandas as pd

PAIR_COLUMNS = ["protein_a", "protein_b"]


def canonical_pair_index(identifiers: Iterable[str]) -> pd.DataFrame:
    """Return one lexicographically canonical row for every unique pair."""

    names = list(identifiers)
    if not names:
        raise ValueError("at least one protein identifier is required")
    if any(not isinstance(name, str) or not name for name in names):
        raise ValueError("protein identifiers must be non-empty strings")
    if len(names) != len(set(names)):
        raise ValueError("protein identifiers must be unique")
    pairs = [sorted(pair) for pair in combinations(names, 2)]
    return pd.DataFrame(pairs, columns=PAIR_COLUMNS).sort_values(
        PAIR_COLUMNS, ignore_index=True
    )


def score_matrix_to_pairs(
    matrix: pd.DataFrame,
    *,
    score_column: str,
) -> pd.DataFrame:
    """Convert a symmetric labelled score matrix to a canonical pair table."""

    if not score_column or score_column in PAIR_COLUMNS:
        raise ValueError("score_column must be a non-key column name")
    if matrix.empty or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("score matrix must be non-empty and square")
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError("score matrix labels must be unique")
    if set(matrix.index) != set(matrix.columns):
        raise ValueError("score matrix row and column labels must match")
    ordered = matrix.loc[matrix.index, matrix.index]
    values = ordered.to_numpy(dtype=float)
    if not np.allclose(values, values.T, equal_nan=True):
        raise ValueError("score matrix must be symmetric")

    lookup = {
        tuple(sorted((str(ordered.index[i]), str(ordered.index[j])))): values[i, j]
        for i in range(len(ordered.index))
        for j in range(i + 1, len(ordered.index))
    }
    rows = [(*pair, score) for pair, score in sorted(lookup.items())]
    return pd.DataFrame(rows, columns=[*PAIR_COLUMNS, score_column])


def merge_pair_tables(
    pair_index: pd.DataFrame,
    *score_tables: pd.DataFrame,
) -> pd.DataFrame:
    """Left-join method outputs onto one complete canonical pair index."""

    merged = _canonicalise_pair_table(pair_index)
    if list(merged.columns) != PAIR_COLUMNS:
        raise ValueError("pair_index must contain only protein_a and protein_b")
    for table in score_tables:
        canonical = _canonicalise_pair_table(table)
        overlapping = (set(merged.columns) & set(canonical.columns)) - set(PAIR_COLUMNS)
        if overlapping:
            raise ValueError(
                "duplicate score column(s): " + ", ".join(sorted(overlapping))
            )
        merged = merged.merge(
            canonical,
            on=PAIR_COLUMNS,
            how="left",
            validate="one_to_one",
            sort=False,
        )
    return merged


def add_descending_ranks(
    scores: pd.DataFrame,
    score_columns: Sequence[str],
) -> pd.DataFrame:
    """Add average ranks and high-is-one percentiles for each score."""

    _require_score_columns(scores, score_columns)
    ranked = scores.copy()
    for column in score_columns:
        ranked[f"{column}_rank"] = ranked[column].rank(
            method="average", ascending=False, na_option="keep"
        )
        ranked[f"{column}_percentile"] = ranked[column].rank(
            method="average", ascending=True, na_option="keep", pct=True
        )
    return ranked


def spearman_correlations(
    scores: pd.DataFrame,
    score_columns: Sequence[str],
) -> pd.DataFrame:
    """Calculate pairwise Spearman rank correlations without rescaling scores."""

    _require_score_columns(scores, score_columns)
    rows = []
    for score_a, score_b in combinations(score_columns, 2):
        complete = scores[[score_a, score_b]].dropna()
        ranks = complete.rank(method="average")
        rows.append(
            {
                "score_a": score_a,
                "score_b": score_b,
                "spearman_rho": ranks[score_a].corr(ranks[score_b]),
                "n_pairs": len(complete),
            }
        )
    return pd.DataFrame(
        rows,
        columns=["score_a", "score_b", "spearman_rho", "n_pairs"],
    )


def pair_rank_disagreements(
    scores: pd.DataFrame,
    score_a: str,
    score_b: str,
) -> pd.DataFrame:
    """Order complete protein pairs by disagreement in score percentile rank."""

    _require_score_columns(scores, [score_a, score_b])
    missing_keys = set(PAIR_COLUMNS) - set(scores.columns)
    if missing_keys:
        raise ValueError(f"missing pair column(s): {', '.join(sorted(missing_keys))}")
    complete = scores[[*PAIR_COLUMNS, score_a, score_b]].dropna().copy()
    rank_a = f"{score_a}_complete_percentile"
    rank_b = f"{score_b}_complete_percentile"
    complete[rank_a] = complete[score_a].rank(method="average", pct=True)
    complete[rank_b] = complete[score_b].rank(method="average", pct=True)
    complete["absolute_percentile_difference"] = (
        complete[rank_a] - complete[rank_b]
    ).abs()
    return complete.sort_values(
        ["absolute_percentile_difference", *PAIR_COLUMNS],
        ascending=[False, True, True],
        ignore_index=True,
    )


def _canonicalise_pair_table(table: pd.DataFrame) -> pd.DataFrame:
    missing = set(PAIR_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f"missing pair column(s): {', '.join(sorted(missing))}")
    canonical = table.copy()
    endpoints = canonical[PAIR_COLUMNS].astype(str)
    canonical["protein_a"] = endpoints.min(axis=1)
    canonical["protein_b"] = endpoints.max(axis=1)
    if (canonical["protein_a"] == canonical["protein_b"]).any():
        raise ValueError("pair tables cannot contain self-pairs")
    if canonical.duplicated(PAIR_COLUMNS).any():
        raise ValueError("pair table contains duplicate protein pairs")
    return canonical.sort_values(PAIR_COLUMNS, ignore_index=True)


def _require_score_columns(scores: pd.DataFrame, score_columns: Sequence[str]) -> None:
    if len(score_columns) < 2:
        raise ValueError("at least two score columns are required")
    missing = set(score_columns) - set(scores.columns)
    if missing:
        raise ValueError(f"missing score column(s): {', '.join(sorted(missing))}")
    non_numeric = [
        column
        for column in score_columns
        if not pd.api.types.is_numeric_dtype(scores[column])
    ]
    if non_numeric:
        raise TypeError("score columns must be numeric: " + ", ".join(non_numeric))
