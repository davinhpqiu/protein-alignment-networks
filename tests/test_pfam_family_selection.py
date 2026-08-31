import gzip

from scripts.select_pfam_family_panel import (
    parse_pfam_clans,
    select_family_panel,
)


def test_parse_pfam_clans_excludes_unassigned_families():
    content = gzip.compress(
        b"PF00001\tCL0001\tClan_one\tFamily_1\tFirst family\n"
        b"PF00002\t\t\tFamily_2\tUnassigned family\n"
    )
    assert parse_pfam_clans(content) == [
        {
            "family_accession": "PF00001",
            "clan_accession": "CL0001",
            "clan_id": "Clan_one",
            "family_id": "Family_1",
            "family_description": "First family",
        }
    ]


def test_family_panel_selection_is_seeded_and_records_rejections():
    catalogue = []
    details = {}
    for clan_index in range(1, 4):
        for family_index in range(1, 5):
            family = f"PF{clan_index}{family_index:04d}"
            catalogue.append(
                {
                    "family_accession": family,
                    "clan_accession": f"CL{clan_index:04d}",
                    "clan_id": f"Clan_{clan_index}",
                    "family_id": f"Family_{family_index}",
                    "family_description": "Test family",
                }
            )
            details[family] = {
                "entry_type": "domain",
                "reported_seed_count": 50,
            }
    details["PF10001"]["reported_seed_count"] = 2

    first, trace = select_family_panel(
        catalogue,
        clans=2,
        families_per_clan=2,
        min_reported_seed_domains=40,
        seed=7,
        inspect_family=details.__getitem__,
    )
    second, _ = select_family_panel(
        catalogue,
        clans=2,
        families_per_clan=2,
        min_reported_seed_domains=40,
        seed=7,
        inspect_family=details.__getitem__,
    )

    assert first == second
    assert len(first) == 4
    assert len({row["clan_accession"] for row in first}) == 2
    assert all(row["reported_seed_count"] >= 40 for row in first)
    assert any(not row["eligible_family"] for row in trace) or len(trace) == 4
