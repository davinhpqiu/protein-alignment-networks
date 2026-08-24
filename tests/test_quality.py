import math

import pandas as pd
import pytest

from protein_alignment_networks.quality import assess_score_correlations


def correlation_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "score_a": "method_a",
                "score_b": "method_b",
                "spearman_rho": 0.85,
                "n_pairs": 90,
            },
            {
                "score_a": "method_a",
                "score_b": "method_c",
                "spearman_rho": 0.40,
                "n_pairs": 60,
            },
        ]
    )


def test_unconfigured_thresholds_do_not_invent_a_quality_decision():
    result = assess_score_correlations(correlation_rows(), total_pairs=100)
    assert set(result["status"]) == {"not_evaluated"}
    assert result["flags"].eq("").all()


def test_configured_thresholds_flag_correlation_and_overlap_separately():
    result = assess_score_correlations(
        correlation_rows(),
        total_pairs=100,
        minimum_rho=0.7,
        minimum_overlap_fraction=0.8,
        minimum_pairs=50,
    )
    assert result.loc[0, "status"] == "pass"
    assert result.loc[0, "flags"] == ""
    assert result.loc[1, "status"] == "flag"
    assert result.loc[1, "flags"] == "poor_correlation;poor_overlap"


def test_threshold_boundaries_pass():
    result = assess_score_correlations(
        correlation_rows().iloc[[0]],
        total_pairs=100,
        minimum_rho=0.85,
        minimum_overlap_fraction=0.9,
        minimum_pairs=90,
    )
    assert result.loc[0, "status"] == "pass"


def test_undefined_correlation_is_always_flagged():
    rows = correlation_rows().iloc[[0]].copy()
    rows.loc[0, "spearman_rho"] = math.nan
    result = assess_score_correlations(rows, total_pairs=100)
    assert result.loc[0, "flags"] == "undefined_correlation"


@pytest.mark.parametrize(
    ("keyword", "value", "message"),
    [
        ("minimum_rho", 1.1, "between -1 and 1"),
        ("minimum_overlap_fraction", -0.1, "between 0 and 1"),
        ("minimum_pairs", 1, "at least 2"),
    ],
)
def test_invalid_thresholds_are_rejected(keyword, value, message):
    with pytest.raises(ValueError, match=message):
        assess_score_correlations(
            correlation_rows(), total_pairs=100, **{keyword: value}
        )
