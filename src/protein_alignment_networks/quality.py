"""Configurable quality checks for cross-method score comparisons."""

from __future__ import annotations

import math

import pandas as pd

CORRELATION_COLUMNS = ["score_a", "score_b", "spearman_rho", "n_pairs"]


def assess_score_correlations(
    correlations: pd.DataFrame,
    *,
    total_pairs: int,
    minimum_rho: float | None = None,
    minimum_overlap_fraction: float | None = None,
    minimum_pairs: int | None = None,
) -> pd.DataFrame:
    """Flag correlations using caller-supplied, optional quality thresholds.

    Thresholds default to ``None`` deliberately: the software records the
    evidence but does not invent a scientific cutoff before data exploration.
    When no thresholds are configured, valid rows receive ``not_evaluated``.
    Undefined correlations are always flagged because they cannot be assessed.
    """

    missing = set(CORRELATION_COLUMNS) - set(correlations.columns)
    if missing:
        raise ValueError(f"missing correlation column(s): {', '.join(sorted(missing))}")
    if total_pairs < 1:
        raise ValueError("total_pairs must be positive")
    if minimum_rho is not None and not -1.0 <= minimum_rho <= 1.0:
        raise ValueError("minimum_rho must be between -1 and 1")
    if minimum_overlap_fraction is not None and not 0.0 <= minimum_overlap_fraction <= 1.0:
        raise ValueError("minimum_overlap_fraction must be between 0 and 1")
    if minimum_pairs is not None and minimum_pairs < 2:
        raise ValueError("minimum_pairs must be at least 2")

    assessed = correlations[CORRELATION_COLUMNS].copy()
    assessed["overlap_fraction"] = assessed["n_pairs"] / total_pairs
    assessed["minimum_rho"] = minimum_rho
    assessed["minimum_overlap_fraction"] = minimum_overlap_fraction
    assessed["minimum_pairs"] = minimum_pairs
    configured = any(
        threshold is not None
        for threshold in (minimum_rho, minimum_overlap_fraction, minimum_pairs)
    )

    flags: list[str] = []
    statuses: list[str] = []
    for row in assessed.itertuples(index=False):
        row_flags: list[str] = []
        if pd.isna(row.spearman_rho) or not math.isfinite(float(row.spearman_rho)):
            row_flags.append("undefined_correlation")
        elif minimum_rho is not None and row.spearman_rho < minimum_rho:
            row_flags.append("poor_correlation")
        if minimum_pairs is not None and row.n_pairs < minimum_pairs:
            row_flags.append("insufficient_pairs")
        if (
            minimum_overlap_fraction is not None
            and row.overlap_fraction < minimum_overlap_fraction
        ):
            row_flags.append("poor_overlap")

        flags.append(";".join(row_flags))
        if row_flags:
            statuses.append("flag")
        elif configured:
            statuses.append("pass")
        else:
            statuses.append("not_evaluated")

    assessed["status"] = statuses
    assessed["flags"] = flags
    return assessed
