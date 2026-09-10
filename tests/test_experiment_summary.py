import pandas as pd
import pytest

from scripts.summarize_resampled_graph_experiment import (
    add_analysis_eligibility,
    collection_fingerprints,
    comparison_status,
    fully_matched_target_densities,
    paired_method_differences,
    rule_label,
)


def test_collection_fingerprints_ignore_membership_order():
    membership = pd.DataFrame(
        {
            "collection_id": ["a", "a", "b", "b"],
            "node_id": ["x", "y", "y", "x"],
        }
    )
    fingerprints = collection_fingerprints(membership).set_index("collection_id")
    assert fingerprints.loc["a", "node_set_sha256"] == fingerprints.loc[
        "b", "node_set_sha256"
    ]


def test_paired_method_differences_pair_exact_graph_designs():
    rows = []
    for method, nmi, density in (("biopython", 0.8, 0.1), ("blast", 0.75, 0.08)):
        rows.append(
            {
                "collection_id": "c1",
                "clans_per_collection": 2,
                "families_per_clan": 2,
                "domains_per_family": 10,
                "node_count": 40,
                "threshold_rule": "top_k",
                "threshold": 5,
                "reference_label": "family_accession",
                "method": method,
                "nmi": nmi,
                "ami": nmi - 0.1,
                "density": density,
                "edges": 100,
                "connected_components": 1,
            }
        )
    paired = paired_method_differences(pd.DataFrame(rows))
    assert paired.loc[0, "nmi_blast_minus_biopython"] == pytest.approx(-0.05)
    assert paired.loc[0, "density_blast_minus_biopython"] == pytest.approx(-0.02)


def test_analysis_eligibility_is_applied_by_composition():
    collections = pd.DataFrame(
        {
            "collection_id": ["eligible", "audit"],
            "clans_per_collection": [2, 2],
            "families_per_clan": [2, 3],
            "domains_per_family": [10, 40],
        }
    )
    uniqueness = collections.drop(columns="collection_id").copy()
    uniqueness["analysis_eligible"] = [True, False]
    results = pd.DataFrame(
        {"collection_id": ["eligible", "eligible", "audit", "audit"]}
    )

    labelled, collection_status = add_analysis_eligibility(
        results, collections, uniqueness
    )

    assert labelled["analysis_eligible"].tolist() == [True, True, False, False]
    assert collection_status["analysis_eligible"].tolist() == [True, False]


@pytest.mark.parametrize(
    ("rule", "threshold", "expected"),
    [
        ("top_k", 5.0, "top-k 5"),
        ("percentile", 0.95, "percentile 0.95"),
        ("target_density", 0.02, "density 2%"),
    ],
)
def test_rule_label_is_human_readable(rule, threshold, expected):
    assert rule_label(rule, threshold) == expected


def test_fully_matched_target_densities_are_derived_from_shortfall():
    results = pd.DataFrame(
        {
            "threshold_rule": ["target_density"] * 4 + ["top_k"],
            "threshold": [0.01, 0.01, 0.05, 0.05, 5.0],
            "selection_shortfall": [0, 0, 0, 3, float("nan")],
        }
    )

    assert fully_matched_target_densities(results) == [0.01]


@pytest.mark.parametrize(
    ("rule", "threshold", "expected"),
    [
        ("target_density", 0.01, "primary_matched_edge_budget"),
        ("target_density", 0.05, "stress_test_edge_budget_shortfall"),
        ("percentile", 0.95, "diagnostic_unmatched_available_fraction"),
        ("top_k", 5.0, "local_neighbour_rule_realised_density"),
    ],
)
def test_comparison_status_describes_matching_scope(rule, threshold, expected):
    assert comparison_status(rule, threshold, [0.01, 0.02]) == expected
