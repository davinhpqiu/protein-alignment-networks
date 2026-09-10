"""Generate reproducible balanced Pfam-domain collections for graph experiments.

The script crosses all requested clan/family/domain design levels and creates a
stable seeded draw for every replicate.  It writes a collection manifest, long
node-membership table, eligible-family frame, and checksum manifest; these are
the authoritative definitions loaded and audited in Notebook 03.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from protein_alignment_networks.sampling import (
    family_sampling_frame,
    repeated_balanced_collections,
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path, help="prepared Pfam metadata TSV")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--clans", nargs="+", type=int, default=[2, 3])
    parser.add_argument("--families-per-clan", nargs="+", type=int, default=[2])
    parser.add_argument("--domains-per-family", nargs="+", type=int, default=[5])
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260828)
    arguments = parser.parse_args()

    metadata_path = arguments.metadata.resolve()
    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(metadata_path, sep="\t", dtype={"node_id": str})

    largest_family_count = max(arguments.families_per_clan)
    largest_domain_count = max(arguments.domains_per_family)
    frame = family_sampling_frame(
        metadata,
        min_domains_per_family=largest_domain_count,
        min_families_per_clan=largest_family_count,
    )
    collections, membership = repeated_balanced_collections(
        metadata,
        clans_per_collection=arguments.clans,
        families_per_clan=arguments.families_per_clan,
        domains_per_family=arguments.domains_per_family,
        replicates=arguments.replicates,
        base_seed=arguments.seed,
    )

    frame_path = output_directory / "family_sampling_frame.tsv"
    collections_path = output_directory / "collection_manifest.tsv"
    membership_path = output_directory / "collection_membership.tsv"
    manifest_path = output_directory / "experiment_manifest.json"
    frame.to_csv(frame_path, sep="\t", index=False)
    collections.to_csv(collections_path, sep="\t", index=False)
    membership.to_csv(membership_path, sep="\t", index=False)

    manifest = {
        "experiment": "repeated balanced Pfam-domain sampling",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "sampling_unit": "Pfam seed domain instance",
        "selection_order": "clans, then families within clans, then domains",
        "replacement": False,
        "balanced_within_collection": True,
        "source_metadata": str(metadata_path),
        "source_metadata_sha256": sha256_file(metadata_path),
        "base_seed": arguments.seed,
        "clans_per_collection": arguments.clans,
        "families_per_clan": arguments.families_per_clan,
        "domains_per_family": arguments.domains_per_family,
        "replicates_per_condition": arguments.replicates,
        "collection_count": len(collections),
        "files": {
            path.name: sha256_file(path)
            for path in (frame_path, collections_path, membership_path)
        },
        "generator": "scripts/generate_pfam_collections.py",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Generated {len(collections)} collections containing "
        f"{len(membership)} total sampled memberships"
    )
    print(f"Collection manifest: {collections_path}")
    print(f"Membership table: {membership_path}")


if __name__ == "__main__":
    main()
