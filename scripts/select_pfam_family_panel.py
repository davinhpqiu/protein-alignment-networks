"""Select a reproducible multi-clan Pfam family panel from the current release."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import random
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path

PFAM_CLANS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
    "Pfam-A.clans.tsv.gz"
)
PFAM_VERSION_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz"
)
ALLOWED_ENTRY_TYPES = {"domain", "family"}


def fetch(url: str, *, timeout: int = 300, attempts: int = 3) -> tuple[bytes, dict]:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "protein-alignment-networks/0.1 "
                "(https://github.com/davinhpqiu/protein-alignment-networks)"
            )
        },
    )
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                headers = {
                    key.lower(): value for key, value in response.headers.items()
                }
                return response.read(), headers
        except (TimeoutError, urllib.error.URLError) as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"download failed after {attempts} attempts: {url}") from last_error


def load_or_fetch(
    url: str,
    content_path: Path,
    headers_path: Path,
    *,
    refresh: bool,
) -> tuple[bytes, dict]:
    if content_path.exists() and headers_path.exists() and not refresh:
        return content_path.read_bytes(), json.loads(headers_path.read_text())
    content, headers = fetch(url)
    content_path.parent.mkdir(parents=True, exist_ok=True)
    content_path.write_bytes(content)
    record = {
        "url": url,
        "downloaded_at_utc": datetime.now(UTC).isoformat(),
        "response_headers": headers,
    }
    headers_path.write_text(json.dumps(record, indent=2) + "\n")
    return content, record


def parse_pfam_clans(content: bytes) -> list[dict[str, str]]:
    """Parse the official five-column Pfam family/clan catalogue."""

    rows = csv.reader(
        io.StringIO(gzip.decompress(content).decode("utf-8")), delimiter="\t"
    )
    catalogue = []
    for row_number, row in enumerate(rows, start=1):
        if len(row) != 5:
            raise ValueError(
                f"unexpected Pfam-A.clans row {row_number}: expected 5 columns"
            )
        family_accession, clan_accession, clan_id, family_id, description = row
        if not family_accession.startswith("PF"):
            raise ValueError(f"invalid Pfam family accession: {family_accession}")
        if not clan_accession:
            continue
        catalogue.append(
            {
                "family_accession": family_accession,
                "clan_accession": clan_accession,
                "clan_id": clan_id,
                "family_id": family_id,
                "family_description": description,
            }
        )
    return catalogue


def select_family_panel(
    catalogue: list[dict[str, str]],
    *,
    clans: int,
    families_per_clan: int,
    min_reported_seed_domains: int,
    seed: int,
    inspect_family: Callable[[str], dict],
) -> tuple[list[dict], list[dict]]:
    """Select clans and families by seeded rejection sampling.

    Candidate clans and their member families are independently shuffled. A
    clan is accepted when the first requested number of eligible families has
    been found. The inspection trace includes every family examined before the
    target panel was complete.
    """

    if clans < 1 or families_per_clan < 1 or min_reported_seed_domains < 1:
        raise ValueError("panel sizes and minimum seed count must be positive")
    by_clan: dict[str, list[dict[str, str]]] = {}
    for row in catalogue:
        by_clan.setdefault(row["clan_accession"], []).append(row)
    candidate_clans = sorted(
        clan
        for clan, rows in by_clan.items()
        if len(rows) >= families_per_clan
    )
    rng = random.Random(seed)
    rng.shuffle(candidate_clans)

    selected: list[dict] = []
    trace: list[dict] = []
    selected_clan_count = 0
    for clan_order, clan_accession in enumerate(candidate_clans, start=1):
        candidates = sorted(
            by_clan[clan_accession], key=lambda row: row["family_accession"]
        )
        rng.shuffle(candidates)
        accepted: list[dict] = []
        clan_trace_start = len(trace)
        for family_order, family in enumerate(candidates, start=1):
            details = inspect_family(family["family_accession"])
            entry_type = details["entry_type"]
            reported_seed_count = int(details["reported_seed_count"])
            eligible = (
                entry_type in ALLOWED_ENTRY_TYPES
                and reported_seed_count >= min_reported_seed_domains
            )
            trace.append(
                {
                    "candidate_clan_order": clan_order,
                    "candidate_family_order": family_order,
                    **family,
                    "entry_type": entry_type,
                    "reported_seed_count": reported_seed_count,
                    "eligible_family": eligible,
                    "selected_clan": False,
                    "selected_family": False,
                }
            )
            if eligible:
                accepted.append({**family, **details})
            if len(accepted) == families_per_clan:
                break

        if len(accepted) < families_per_clan:
            continue
        selected_clan_count += 1
        for trace_row in trace[clan_trace_start:]:
            trace_row["selected_clan"] = True
        accepted_accessions = {
            family["family_accession"] for family in accepted
        }
        for trace_row in trace[clan_trace_start:]:
            trace_row["selected_family"] = (
                trace_row["family_accession"] in accepted_accessions
            )
        for family_rank, family in enumerate(accepted, start=1):
            selected.append(
                {
                    "selected_clan_rank": selected_clan_count,
                    "selected_family_rank": family_rank,
                    **family,
                }
            )
        if selected_clan_count == clans:
            break

    if selected_clan_count < clans:
        raise ValueError(
            f"only {selected_clan_count} clans met the requested eligibility rules"
        )
    return selected, trace


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(destination, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--clans", type=int, default=8)
    parser.add_argument("--families-per-clan", type=int, default=3)
    parser.add_argument("--min-reported-seed-domains", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--refresh", action="store_true")
    arguments = parser.parse_args()

    root = arguments.project_root.resolve()
    external = root / "data/external/pfam_sampling_frame"
    processed = root / "data/processed/pfam_sampling_frame"
    external.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)

    catalogue_content, catalogue_record = load_or_fetch(
        PFAM_CLANS_URL,
        external / "Pfam-A.clans.tsv.gz",
        external / "Pfam-A.clans.headers.json",
        refresh=arguments.refresh,
    )
    version_content, _version_record = load_or_fetch(
        PFAM_VERSION_URL,
        external / "Pfam.version.gz",
        external / "Pfam.version.headers.json",
        refresh=arguments.refresh,
    )
    catalogue = parse_pfam_clans(catalogue_content)

    def inspect_family(family_accession: str) -> dict:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{family_accession}/"
        content, _record = load_or_fetch(
            url,
            external / f"{family_accession}.interpro-entry.json",
            external / f"{family_accession}.interpro-entry.headers.json",
            refresh=arguments.refresh,
        )
        metadata = json.loads(content)["metadata"]
        return {
            "entry_type": metadata["type"],
            "reported_seed_count": metadata["entry_annotations"]["alignment:seed"],
        }

    selected, trace = select_family_panel(
        catalogue,
        clans=arguments.clans,
        families_per_clan=arguments.families_per_clan,
        min_reported_seed_domains=arguments.min_reported_seed_domains,
        seed=arguments.seed,
        inspect_family=inspect_family,
    )
    selected_path = processed / "selected_families.tsv"
    trace_path = processed / "family_selection_trace.tsv"
    manifest_path = processed / "selection_manifest.json"
    write_tsv(selected_path, selected)
    write_tsv(trace_path, trace)

    version_text = gzip.decompress(version_content).decode("utf-8").strip()
    manifest = {
        "dataset": "reproducibly selected multi-clan Pfam family panel",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "pfam_version_file": version_text,
        "sampling_unit": "Pfam family nested within Pfam clan",
        "selection_rule": (
            "Exclude unassigned families and clans with fewer than the requested "
            "number of catalogue families. Shuffle clans with the fixed seed; "
            "within each clan shuffle families; accept domain/family entries "
            "whose reported Pfam seed count meets the minimum; take the first "
            "requested number of eligible families and the first requested "
            "number of eligible clans."
        ),
        "seed": arguments.seed,
        "requested_clans": arguments.clans,
        "requested_families_per_clan": arguments.families_per_clan,
        "min_reported_seed_domains": arguments.min_reported_seed_domains,
        "selected_family_count": len(selected),
        "catalogue_url": catalogue_record["url"],
        "catalogue_sha256": sha256_bytes(catalogue_content),
        "selected_families_sha256": sha256_bytes(selected_path.read_bytes()),
        "selection_trace_sha256": sha256_bytes(trace_path.read_bytes()),
        "generator": "scripts/select_pfam_family_panel.py",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Selected {len(selected)} families from {arguments.clans} clans; "
        f"inspected {len(trace)} candidate families"
    )
    print(f"Selected families: {selected_path}")
    print(f"Selection trace: {trace_path}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
