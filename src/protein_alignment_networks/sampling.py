"""Reproducible sampling designs for labelled protein-domain collections."""

from __future__ import annotations

import hashlib
import json
import random
from collections.abc import Sequence
from itertools import product

import pandas as pd

REQUIRED_METADATA_COLUMNS = {
    "node_id",
    "family_accession",
    "family_name",
    "clan_accession",
    "clan_name",
}


def family_sampling_frame(
    metadata: pd.DataFrame,
    *,
    min_domains_per_family: int,
    min_families_per_clan: int,
) -> pd.DataFrame:
    """Describe which prepared Pfam families and clans can enter a design.

    The returned table retains excluded families and records why they cannot be
    sampled. This makes the sampling frame auditable rather than silently
    dropping small families or clans.
    """

    if min_domains_per_family < 1:
        raise ValueError("min_domains_per_family must be positive")
    if min_families_per_clan < 1:
        raise ValueError("min_families_per_clan must be positive")
    _validate_metadata(metadata)

    frame = (
        metadata.groupby(
            [
                "clan_accession",
                "clan_name",
                "family_accession",
                "family_name",
            ],
            as_index=False,
            sort=True,
        )
        .agg(available_domains=("node_id", "nunique"))
        .sort_values(["clan_accession", "family_accession"], ignore_index=True)
    )
    frame["family_has_enough_domains"] = (
        frame["available_domains"] >= min_domains_per_family
    )
    eligible_counts = (
        frame.loc[frame["family_has_enough_domains"]]
        .groupby("clan_accession")["family_accession"]
        .nunique()
    )
    frame["eligible_families_in_clan"] = (
        frame["clan_accession"].map(eligible_counts).fillna(0).astype(int)
    )
    frame["clan_has_enough_families"] = (
        frame["eligible_families_in_clan"] >= min_families_per_clan
    ) & (frame["clan_accession"].str.lower() != "unassigned")
    frame["eligible"] = (
        frame["family_has_enough_domains"] & frame["clan_has_enough_families"]
    )
    frame["exclusion_reason"] = ""
    frame.loc[
        ~frame["family_has_enough_domains"], "exclusion_reason"
    ] = "too_few_domains"
    frame.loc[
        frame["family_has_enough_domains"] & ~frame["clan_has_enough_families"],
        "exclusion_reason",
    ] = "too_few_eligible_families_in_clan"
    return frame


def sample_balanced_collection(
    metadata: pd.DataFrame,
    *,
    clans_per_collection: int,
    families_per_clan: int,
    domains_per_family: int,
    seed: int,
) -> pd.DataFrame:
    """Sample one balanced collection of Pfam domain instances.

    Clans, then families within clans, then nodes within families are sampled
    uniformly without replacement from sorted candidate lists. Pfam labels are
    returned for later evaluation but are not used when constructing graphs or
    inferring communities.
    """

    if clans_per_collection < 1:
        raise ValueError("clans_per_collection must be positive")
    if families_per_clan < 1:
        raise ValueError("families_per_clan must be positive")
    if domains_per_family < 1:
        raise ValueError("domains_per_family must be positive")
    _validate_metadata(metadata)

    frame = family_sampling_frame(
        metadata,
        min_domains_per_family=domains_per_family,
        min_families_per_clan=families_per_clan,
    )
    eligible = frame.loc[frame["eligible"]].copy()
    candidate_clans = sorted(eligible["clan_accession"].unique())
    if len(candidate_clans) < clans_per_collection:
        raise ValueError(
            f"requested {clans_per_collection} clans, but only "
            f"{len(candidate_clans)} are eligible"
        )

    rng = random.Random(seed)
    selected_clans = rng.sample(candidate_clans, clans_per_collection)
    rows: list[dict] = []
    for clan_rank, clan_accession in enumerate(selected_clans, start=1):
        clan_frame = eligible.loc[
            eligible["clan_accession"] == clan_accession
        ]
        candidate_families = sorted(clan_frame["family_accession"].unique())
        selected_families = rng.sample(candidate_families, families_per_clan)
        for family_rank, family_accession in enumerate(selected_families, start=1):
            family_nodes = metadata.loc[
                metadata["family_accession"] == family_accession
            ].sort_values("node_id")
            selected_indices = rng.sample(
                list(range(len(family_nodes))), domains_per_family
            )
            selected_nodes = family_nodes.iloc[selected_indices]
            for domain_rank, (_, node) in enumerate(
                selected_nodes.iterrows(), start=1
            ):
                rows.append(
                    {
                        "node_id": node["node_id"],
                        "family_accession": node["family_accession"],
                        "family_name": node["family_name"],
                        "clan_accession": node["clan_accession"],
                        "clan_name": node["clan_name"],
                        "clan_sample_rank": clan_rank,
                        "family_sample_rank": family_rank,
                        "domain_sample_rank": domain_rank,
                        "sampling_seed": seed,
                    }
                )
    return pd.DataFrame(rows)


def repeated_balanced_collections(
    metadata: pd.DataFrame,
    *,
    clans_per_collection: Sequence[int],
    families_per_clan: Sequence[int],
    domains_per_family: Sequence[int],
    replicates: int,
    base_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate a factorial series of balanced, reproducible collections.

    Returns a collection-level manifest and a long membership table. A stable
    hash derives each replicate seed from its full condition, so extending the
    design does not alter collections that were generated previously.
    """

    if replicates < 1:
        raise ValueError("replicates must be positive")
    for name, values in (
        ("clans_per_collection", clans_per_collection),
        ("families_per_clan", families_per_clan),
        ("domains_per_family", domains_per_family),
    ):
        if not values or any(value < 1 for value in values):
            raise ValueError(f"{name} must contain positive integers")

    manifest_rows: list[dict] = []
    membership_tables: list[pd.DataFrame] = []
    conditions = product(
        clans_per_collection, families_per_clan, domains_per_family
    )
    for clan_count, family_count, domain_count in conditions:
        for replicate in range(1, replicates + 1):
            collection_id = (
                f"c{clan_count}_f{family_count}_d{domain_count}_"
                f"r{replicate:03d}"
            )
            replicate_seed = _condition_seed(
                base_seed,
                clan_count,
                family_count,
                domain_count,
                replicate,
            )
            membership = sample_balanced_collection(
                metadata,
                clans_per_collection=clan_count,
                families_per_clan=family_count,
                domains_per_family=domain_count,
                seed=replicate_seed,
            )
            membership.insert(0, "collection_id", collection_id)
            membership.insert(1, "replicate", replicate)
            membership_tables.append(membership)

            clans = sorted(membership["clan_accession"].unique())
            families = sorted(membership["family_accession"].unique())
            node_count = len(membership)
            manifest_rows.append(
                {
                    "collection_id": collection_id,
                    "replicate": replicate,
                    "sampling_unit": "pfam_domain_instance",
                    "base_seed": base_seed,
                    "sampling_seed": replicate_seed,
                    "clans_per_collection": clan_count,
                    "families_per_clan": family_count,
                    "domains_per_family": domain_count,
                    "node_count": node_count,
                    "unique_pair_count": node_count * (node_count - 1) // 2,
                    "clan_accessions": json.dumps(clans),
                    "family_accessions": json.dumps(families),
                }
            )

    manifest = pd.DataFrame(manifest_rows)
    membership = pd.concat(membership_tables, ignore_index=True)
    return manifest, membership


def _condition_seed(
    base_seed: int,
    clan_count: int,
    family_count: int,
    domain_count: int,
    replicate: int,
) -> int:
    payload = (
        f"{base_seed}|{clan_count}|{family_count}|{domain_count}|{replicate}"
    ).encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def _validate_metadata(metadata: pd.DataFrame) -> None:
    missing = REQUIRED_METADATA_COLUMNS - set(metadata.columns)
    if missing:
        raise ValueError(f"metadata is missing columns: {', '.join(sorted(missing))}")
    if metadata.empty:
        raise ValueError("metadata must not be empty")
    if metadata[list(REQUIRED_METADATA_COLUMNS)].isna().any().any():
        raise ValueError("sampling metadata must not contain missing labels")
    if metadata["node_id"].duplicated().any():
        raise ValueError("node_id values must be unique")
    family_clans = metadata.groupby("family_accession")["clan_accession"].nunique()
    if (family_clans > 1).any():
        raise ValueError("each family must map to exactly one clan")
