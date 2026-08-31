"""Prepare a diverse multi-family Pfam seed collection for graph exploration."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import statistics
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from Bio import AlignIO

DEFAULT_FAMILIES = (
    "PF00042",  # Globin
    "PF01152",  # Bacterial-like globin
    "PF00069",  # Protein kinase domain
    "PF07714",  # Protein tyrosine/serine-threonine kinase
    "PF00018",  # SH3 domain
    "PF08239",  # Bacterial SH3 domain
)
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
UNIPROT_FIELDS = (
    "accession",
    "id",
    "reviewed",
    "protein_name",
    "organism_name",
    "organism_id",
    "length",
    "sequence_version",
)


def sha256_bytes(content: bytes) -> str:
    """Return the SHA-256 checksum of bytes."""

    return hashlib.sha256(content).hexdigest()


def fetch(url: str, *, timeout: int = 300, attempts: int = 3) -> tuple[bytes, dict]:
    """Download one resource with bounded retries and a project user agent."""

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
    """Reuse a recorded source response unless an explicit refresh is requested."""

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


def parse_seed_alignment(content: bytes, pfam_accession: str) -> tuple[str, list[dict]]:
    """Parse valid, unique domain sequences from one Pfam seed alignment."""

    uncompressed = gzip.decompress(content).decode("utf-8")
    pfam_version = ""
    for line in uncompressed.splitlines():
        if line.startswith("#=GF AC"):
            pfam_version = line.split(maxsplit=2)[-1]
            break
    if not pfam_version.startswith(f"{pfam_accession}."):
        raise ValueError(f"unexpected Pfam accession in seed: {pfam_version}")

    alignment = AlignIO.read(io.StringIO(uncompressed), "stockholm")
    candidates: list[dict] = []
    seen_sequences: set[str] = set()
    for record in alignment:
        versioned_accession = record.annotations.get("accession", "")
        accession = versioned_accession.split(".")[0]
        if not accession:
            accession = record.id.split("/")[0]
        start = int(record.annotations["start"])
        end = int(record.annotations["end"])
        aligned_sequence = str(record.seq).upper()
        sequence = aligned_sequence.replace("-", "").replace(".", "")
        if not sequence or set(sequence) - STANDARD_AMINO_ACIDS:
            continue
        if len(sequence) != end - start + 1 or sequence in seen_sequences:
            continue
        seen_sequences.add(sequence)
        candidates.append(
            {
                "node_id": f"{accession}_{start}_{end}",
                "uniprot_accession": accession,
                "pfam_sequence_id": record.id,
                "pfam_sequence_accession": versioned_accession,
                "domain_start": start,
                "domain_end": end,
                "domain_length": len(sequence),
                "sequence": sequence,
                "aligned_sequence": aligned_sequence,
            }
        )
    if not candidates:
        raise ValueError(f"no eligible seed sequences found for {pfam_accession}")
    return pfam_version, candidates


def aligned_distance(sequence_a: str, sequence_b: str) -> float:
    """Return one minus identity over positions occupied by both sequences."""

    overlap = 0
    matches = 0
    for residue_a, residue_b in zip(sequence_a, sequence_b, strict=True):
        if residue_a in STANDARD_AMINO_ACIDS and residue_b in STANDARD_AMINO_ACIDS:
            overlap += 1
            matches += residue_a == residue_b
    return 1.0 if overlap == 0 else 1.0 - matches / overlap


def select_diverse(candidates: list[dict], count: int) -> list[dict]:
    """Select deterministic farthest-point representatives from an aligned seed."""

    if count < 1:
        raise ValueError("count must be positive")
    if len(candidates) <= count:
        return sorted(candidates, key=lambda row: row["node_id"])

    median_length = statistics.median(row["domain_length"] for row in candidates)
    first = min(
        candidates,
        key=lambda row: (abs(row["domain_length"] - median_length), row["node_id"]),
    )
    selected = [first]
    remaining = [row for row in candidates if row is not first]
    while len(selected) < count:
        next_row = max(
            remaining,
            key=lambda row: (
                min(
                    aligned_distance(
                        row["aligned_sequence"], chosen["aligned_sequence"]
                    )
                    for chosen in selected
                ),
                row["node_id"],
            ),
        )
        selected.append(next_row)
        remaining.remove(next_row)
    return selected


def uniprot_url(accessions: list[str]) -> str:
    """Build one UniProt metadata request for selected accessions."""

    query = urllib.parse.urlencode(
        {
            "accessions": ",".join(accessions),
            "format": "tsv",
            "fields": ",".join(UNIPROT_FIELDS),
        }
    )
    return f"https://rest.uniprot.org/uniprotkb/accessions?{query}"


def parse_uniprot_metadata(content: bytes) -> dict[str, dict[str, str]]:
    """Parse a UniProt TSV response by accession."""

    rows = csv.DictReader(io.StringIO(content.decode("utf-8")), delimiter="\t")
    return {row["Entry"]: row for row in rows}


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(
        sequence[index : index + width]
        for index in range(0, len(sequence), width)
    )


def relative(path: Path, root: Path) -> str:
    return str(path.resolve().relative_to(root.resolve()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    family_group = parser.add_mutually_exclusive_group()
    family_group.add_argument(
        "--families",
        nargs="+",
        default=None,
        help="Pfam accessions to include",
    )
    family_group.add_argument(
        "--family-file",
        type=Path,
        help="TSV containing a family_accession column",
    )
    parser.add_argument(
        "--dataset-id",
        default="pfam_clan_panel",
        help="safe directory and output-file stem",
    )
    parser.add_argument(
        "--per-family",
        type=int,
        default=40,
        help="maximum diverse seed domains selected from each family",
    )
    parser.add_argument("--refresh", action="store_true")
    arguments = parser.parse_args()
    if arguments.per_family < 1:
        raise ValueError("--per-family must be positive")
    if not arguments.dataset_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("--dataset-id may contain only letters, numbers, _ and -")
    if arguments.family_file:
        with arguments.family_file.open(encoding="utf-8") as source:
            arguments.families = [
                row["family_accession"] for row in csv.DictReader(source, delimiter="\t")
            ]
        if not arguments.families:
            raise ValueError("--family-file contains no families")
    elif arguments.families is None:
        arguments.families = list(DEFAULT_FAMILIES)

    root = arguments.project_root.resolve()
    external = root / "data/external" / arguments.dataset_id
    processed = root / "data/processed" / arguments.dataset_id
    external.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)

    selected_rows: list[dict] = []
    family_records: list[dict] = []
    source_checksums: dict[str, str] = {}
    interpro_versions: set[str] = set()
    for family in arguments.families:
        entry_url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{family}/"
        seed_url = entry_url + "?annotation=alignment:seed&download"
        entry_path = external / f"{family}.interpro-entry.json"
        seed_path = external / f"{family}.seed.sto.gz"
        entry_content, entry_record = load_or_fetch(
            entry_url,
            entry_path,
            external / f"{family}.interpro-entry.headers.json",
            refresh=arguments.refresh,
        )
        seed_content, _seed_record = load_or_fetch(
            seed_url,
            seed_path,
            external / f"{family}.seed.headers.json",
            refresh=arguments.refresh,
        )
        entry = json.loads(entry_content)
        metadata = entry["metadata"]
        pfam_version, candidates = parse_seed_alignment(seed_content, family)
        selected = select_diverse(candidates, arguments.per_family)
        clan = metadata.get("set_info") or {}
        for rank, row in enumerate(selected, start=1):
            selected_rows.append(
                {
                    **row,
                    "family_accession": family,
                    "family_version": pfam_version,
                    "family_name": metadata["name"]["name"],
                    "clan_accession": clan.get("accession", "unassigned"),
                    "clan_name": clan.get("name", "Unassigned"),
                    "family_selection_rank": rank,
                }
            )
        response_headers = entry_record.get("response_headers", {})
        if response_headers.get("interpro-version"):
            interpro_versions.add(response_headers["interpro-version"])
        family_records.append(
            {
                "family_accession": family,
                "family_version": pfam_version,
                "family_name": metadata["name"]["name"],
                "clan_accession": clan.get("accession"),
                "clan_name": clan.get("name"),
                "reported_seed_count": metadata["entry_annotations"][
                    "alignment:seed"
                ],
                "eligible_unique_seed_count": len(candidates),
                "selected_count": len(selected),
            }
        )
        source_checksums[relative(entry_path, root)] = sha256_bytes(entry_content)
        source_checksums[relative(seed_path, root)] = sha256_bytes(seed_content)

    accessions = sorted({row["uniprot_accession"] for row in selected_rows})
    metadata_url = uniprot_url(accessions)
    uniprot_path = external / "selected-uniprot.tsv"
    uniprot_content, uniprot_record = load_or_fetch(
        metadata_url,
        uniprot_path,
        external / "selected-uniprot.headers.json",
        refresh=arguments.refresh,
    )
    uniprot = parse_uniprot_metadata(uniprot_content)
    source_checksums[relative(uniprot_path, root)] = sha256_bytes(uniprot_content)

    fasta_parts: list[str] = []
    metadata_rows: list[dict] = []
    for selection_index, row in enumerate(selected_rows, start=1):
        accession = row["uniprot_accession"]
        record = uniprot.get(accession, {})
        fasta_parts.append(
            f">{row['node_id']} uniprot={accession} "
            f"family={row['family_accession']} clan={row['clan_accession']}\n"
            f"{wrap_fasta(row['sequence'])}\n"
        )
        metadata_rows.append(
            {
                "selection_index": selection_index,
                "node_id": row["node_id"],
                "uniprot_accession": accession,
                "uniprot_entry_name": record.get("Entry Name", ""),
                "reviewed": record.get("Reviewed", "unknown"),
                "protein_name": record.get("Protein names", ""),
                "organism": record.get("Organism", ""),
                "taxonomy_id": record.get("Organism (ID)", ""),
                "full_protein_length": record.get("Length", ""),
                "uniprot_sequence_version": record.get("Sequence version", ""),
                "family_accession": row["family_accession"],
                "family_version": row["family_version"],
                "family_name": row["family_name"],
                "clan_accession": row["clan_accession"],
                "clan_name": row["clan_name"],
                "family_selection_rank": row["family_selection_rank"],
                "pfam_sequence_id": row["pfam_sequence_id"],
                "pfam_sequence_accession": row["pfam_sequence_accession"],
                "domain_start": row["domain_start"],
                "domain_end": row["domain_end"],
                "domain_length": row["domain_length"],
                "sequence_sha256": sha256_bytes(row["sequence"].encode("ascii")),
            }
        )

    fasta_path = processed / f"{arguments.dataset_id}.fasta"
    metadata_path = processed / f"{arguments.dataset_id}_metadata.tsv"
    manifest_path = processed / f"{arguments.dataset_id}_manifest.json"
    fasta_path.write_text("".join(fasta_parts))
    with metadata_path.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(
            destination, fieldnames=list(metadata_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(metadata_rows)

    uniprot_headers = uniprot_record.get("response_headers", {})
    output_checksums = {
        relative(fasta_path, root): sha256_bytes(fasta_path.read_bytes()),
        relative(metadata_path, root): sha256_bytes(metadata_path.read_bytes()),
    }
    manifest = {
        "dataset": "Pfam clan/family exploration panel",
        "dataset_id": arguments.dataset_id,
        "created_at_utc": datetime.now(UTC).isoformat(),
        "selection_rule": (
            "For each requested Pfam family, filter the curated seed to unique "
            "standard-amino-acid domain sequences and apply deterministic "
            "farthest-point sampling using seed-alignment identity. Select at "
            "most the requested number per family; smaller seeds contribute all "
            "eligible sequences. The Pfam alignment is used only for sampling "
            "and is not supplied to the scoring methods."
        ),
        "target_per_family": arguments.per_family,
        "selected_sequence_count": len(metadata_rows),
        "unique_pair_count": len(metadata_rows) * (len(metadata_rows) - 1) // 2,
        "families": family_records,
        "interpro_releases": sorted(interpro_versions),
        "uniprot_release": uniprot_headers.get("x-uniprot-release"),
        "uniprot_release_date": uniprot_headers.get("x-uniprot-release-date"),
        "licence": "CC0 1.0 Universal for InterPro/Pfam downloadable data",
        "generator": "scripts/prepare_pfam_clan_panel.py",
        "sources": source_checksums,
        "files": output_checksums,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Prepared {len(metadata_rows)} domain sequences and "
        f"{manifest['unique_pair_count']} unique pairs"
    )
    for family in family_records:
        print(
            f"  {family['family_accession']} {family['family_name']}: "
            f"{family['selected_count']} selected"
        )
    print(f"FASTA: {fasta_path}")
    print(f"Metadata: {metadata_path}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
