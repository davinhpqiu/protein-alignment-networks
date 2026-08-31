import json

import pandas as pd
import pytest

from protein_alignment_networks.sampling import (
    family_sampling_frame,
    repeated_balanced_collections,
    sample_balanced_collection,
)


def panel_metadata() -> pd.DataFrame:
    rows = []
    for clan in ("C1", "C2", "C3"):
        for family_index in range(1, 4):
            family = f"{clan}F{family_index}"
            for domain_index in range(1, 6):
                rows.append(
                    {
                        "node_id": f"{family}D{domain_index}",
                        "family_accession": family,
                        "family_name": f"Family {family}",
                        "clan_accession": clan,
                        "clan_name": f"Clan {clan}",
                    }
                )
    return pd.DataFrame(rows)


def test_family_sampling_frame_records_eligibility_reasons():
    metadata = panel_metadata()
    metadata = metadata.loc[metadata["family_accession"] != "C3F3"]
    frame = family_sampling_frame(
        metadata,
        min_domains_per_family=5,
        min_families_per_clan=3,
    )

    assert frame.loc[frame["clan_accession"] == "C1", "eligible"].all()
    c3 = frame.loc[frame["clan_accession"] == "C3"]
    assert not c3["eligible"].any()
    assert set(c3["exclusion_reason"]) == {"too_few_eligible_families_in_clan"}


def test_balanced_collection_is_reproducible_and_balanced():
    metadata = panel_metadata()
    first = sample_balanced_collection(
        metadata,
        clans_per_collection=2,
        families_per_clan=2,
        domains_per_family=3,
        seed=17,
    )
    second = sample_balanced_collection(
        metadata,
        clans_per_collection=2,
        families_per_clan=2,
        domains_per_family=3,
        seed=17,
    )

    pd.testing.assert_frame_equal(first, second)
    assert len(first) == 12
    assert first["clan_accession"].nunique() == 2
    assert (first.groupby("clan_accession")["family_accession"].nunique() == 2).all()
    assert (first.groupby("family_accession")["node_id"].size() == 3).all()
    assert not first["node_id"].duplicated().any()


def test_repeated_collections_have_stable_ids_and_complete_manifests():
    metadata = panel_metadata()
    manifest, membership = repeated_balanced_collections(
        metadata,
        clans_per_collection=[2],
        families_per_clan=[2],
        domains_per_family=[2, 3],
        replicates=2,
        base_seed=42,
    )

    assert manifest["collection_id"].tolist() == [
        "c2_f2_d2_r001",
        "c2_f2_d2_r002",
        "c2_f2_d3_r001",
        "c2_f2_d3_r002",
    ]
    assert manifest["sampling_seed"].nunique() == 4
    assert manifest["node_count"].tolist() == [8, 8, 12, 12]
    assert set(membership["collection_id"]) == set(manifest["collection_id"])
    assert all(json.loads(value) for value in manifest["family_accessions"])


def test_sampling_rejects_an_impossible_design():
    with pytest.raises(ValueError, match="only 3 are eligible"):
        sample_balanced_collection(
            panel_metadata(),
            clans_per_collection=4,
            families_per_clan=2,
            domains_per_family=2,
            seed=1,
        )
