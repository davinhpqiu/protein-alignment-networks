import pandas as pd
import pytest

from scripts.summarize_resampled_graph_experiment import (
    collection_fingerprints,
    paired_method_differences,
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
