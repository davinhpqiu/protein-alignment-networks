"""Download and prepare the reproducible PF00042 globin pilot dataset."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from Bio import AlignIO

PFAM_ACCESSION = "PF00042"
INTERPRO_ENTRY_URL = (
    "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/"
)
INTERPRO_SEED_URL = (
    "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/"
    "?annotation=alignment:seed&download"
)
UNIPROT_FIELDS = [
    "accession",
    "id",
    "reviewed",
    "protein_name",
    "organism_name",
    "organism_id",
    "length",
    "sequence_version",
]
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# This is a deliberately small, interpretable pilot rather than a random sample.
# All accessions are reviewed UniProtKB members of the PF00042.29 seed alignment.
SELECTION = [
    {
        "accession": "P02204",
        "subgroup": "myoglobin",
        "selection_reason": "ray-finned fish myoglobin",
    },
    {
        "accession": "P02200",
        "subgroup": "myoglobin",
        "selection_reason": "reptile myoglobin",
    },
    {
        "accession": "P02206",
        "subgroup": "myoglobin",
        "selection_reason": "cartilaginous-fish myoglobin",
    },
    {
        "accession": "P02143",
        "subgroup": "beta-like haemoglobin",
        "selection_reason": "cartilaginous-fish beta globin",
    },
    {
        "accession": "P02133",
        "subgroup": "beta-like haemoglobin",
        "selection_reason": "amphibian beta globin",
    },
    {
        "accession": "P04443",
        "subgroup": "beta-like haemoglobin",
        "selection_reason": "mammalian beta-like globin",
    },
    {
        "accession": "P02020",
        "subgroup": "alpha-like haemoglobin",
        "selection_reason": "lungfish alpha globin",
    },
    {
        "accession": "P06714",
        "subgroup": "alpha-like haemoglobin",
        "selection_reason": "mammalian theta globin from the alpha-like cluster",
    },
    {
        "accession": "P01967",
        "subgroup": "alpha-like haemoglobin",
        "selection_reason": "mammalian alpha globin",
    },
    {
        "accession": "P09187",
        "subgroup": "divergent globin",
        "selection_reason": "plant leghemoglobin",
    },
    {
        "accession": "P24232",
        "subgroup": "divergent globin",
        "selection_reason": "bacterial flavohemoglobin globin domain",
    },
    {
        "accession": "P30627",
        "subgroup": "divergent globin",
        "selection_reason": "nematode globin",
    },
]


def sha256_bytes(content: bytes) -> str:
    """Return the SHA-256 checksum of bytes."""

    return hashlib.sha256(content).hexdigest()


def fetch(url: str, *, timeout: int = 300, attempts: int = 3) -> tuple[bytes, dict]:
    """Download a resource with a descriptive user agent and bounded retries."""

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
                headers = {key.lower(): value for key, value in response.headers.items()}
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
    """Reuse a recorded source file or download it with response metadata."""

    if content_path.exists() and headers_path.exists() and not refresh:
        return content_path.read_bytes(), json.loads(headers_path.read_text())

    content, headers = fetch(url)
    content_path.parent.mkdir(parents=True, exist_ok=True)
    content_path.write_bytes(content)
    header_record = {
        "url": url,
        "downloaded_at_utc": datetime.now(UTC).isoformat(),
        "response_headers": headers,
    }
    headers_path.write_text(json.dumps(header_record, indent=2) + "\n")
    return content, header_record


def uniprot_url(accessions: list[str]) -> str:
    """Build the UniProt metadata request for the selected accessions."""

    query = urllib.parse.urlencode(
        {
            "accessions": ",".join(accessions),
            "format": "tsv",
            "fields": ",".join(UNIPROT_FIELDS),
        }
    )
    return f"https://rest.uniprot.org/uniprotkb/accessions?{query}"


def parse_seed_alignment(content: bytes) -> tuple[str, dict[str, dict]]:
    """Parse a gzipped Pfam Stockholm seed alignment by UniProt accession."""

    uncompressed = gzip.decompress(content).decode("utf-8")
    pfam_version = ""
    for line in uncompressed.splitlines():
        if line.startswith("#=GF AC"):
            pfam_version = line.split(maxsplit=2)[-1]
            break
    if not pfam_version.startswith(f"{PFAM_ACCESSION}."):
        raise ValueError(f"unexpected Pfam accession in seed alignment: {pfam_version}")

    alignment = AlignIO.read(io.StringIO(uncompressed), "stockholm")
    records: dict[str, dict] = {}
    for record in alignment:
        versioned_accession = record.annotations.get("accession", "")
        accession = versioned_accession.split(".")[0]
        sequence = str(record.seq).replace("-", "").replace(".", "").upper()
        start = int(record.annotations["start"])
        end = int(record.annotations["end"])
        if len(sequence) != end - start + 1:
            raise ValueError(f"domain coordinates do not match sequence length: {record.id}")
        records[accession] = {
            "pfam_sequence_id": record.id,
            "pfam_sequence_accession": versioned_accession,
            "domain_start": start,
            "domain_end": end,
            "sequence": sequence,
        }
    return pfam_version, records


def parse_uniprot_metadata(content: bytes) -> dict[str, dict[str, str]]:
    """Parse the selected UniProt metadata TSV by accession."""

    rows = csv.DictReader(io.StringIO(content.decode("utf-8")), delimiter="\t")
    return {row["Entry"]: row for row in rows}


def wrap_fasta(sequence: str, width: int = 80) -> str:
    """Wrap a FASTA sequence to a stable line width."""

    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def relative(path: Path, root: Path) -> str:
    """Return a portable project-relative path."""

    return str(path.relative_to(root))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="project root containing data/ (default: inferred from this script)",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="redownload source files instead of reusing the recorded copies",
    )
    arguments = parser.parse_args()
    root = arguments.project_root.resolve()
    external = root / "data/external/pfam_pf00042"
    processed = root / "data/processed/pfam_pf00042_globin_pilot"
    external.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)

    entry_content, entry_record = load_or_fetch(
        INTERPRO_ENTRY_URL,
        external / "PF00042.interpro-entry.json",
        external / "PF00042.interpro-entry.headers.json",
        refresh=arguments.refresh,
    )
    seed_content, _seed_record = load_or_fetch(
        INTERPRO_SEED_URL,
        external / "PF00042.seed.sto.gz",
        external / "PF00042.seed.headers.json",
        refresh=arguments.refresh,
    )
    accessions = [item["accession"] for item in SELECTION]
    metadata_url = uniprot_url(accessions)
    uniprot_content, uniprot_record = load_or_fetch(
        metadata_url,
        external / "PF00042.selected-uniprot.tsv",
        external / "PF00042.selected-uniprot.headers.json",
        refresh=arguments.refresh,
    )

    entry = json.loads(entry_content)
    pfam_version, seed_records = parse_seed_alignment(seed_content)
    uniprot_records = parse_uniprot_metadata(uniprot_content)
    missing_seed = sorted(set(accessions) - set(seed_records))
    missing_metadata = sorted(set(accessions) - set(uniprot_records))
    if missing_seed or missing_metadata:
        raise ValueError(
            f"missing selected records; seed={missing_seed}, metadata={missing_metadata}"
        )

    fasta_parts: list[str] = []
    metadata_rows: list[dict] = []
    seen_sequences: set[str] = set()
    for index, selection in enumerate(SELECTION, start=1):
        accession = selection["accession"]
        seed = seed_records[accession]
        metadata = uniprot_records[accession]
        sequence = seed["sequence"]
        invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
        if invalid:
            raise ValueError(f"{accession} contains non-standard residues: {invalid}")
        if sequence in seen_sequences:
            raise ValueError(f"selected domain sequence is duplicated: {accession}")
        if len(sequence) > 511:
            raise ValueError(f"selected domain exceeds DEDAL limit: {accession}")
        if metadata["Reviewed"] != "reviewed":
            raise ValueError(f"selected UniProt entry is not reviewed: {accession}")
        seen_sequences.add(sequence)

        seed_version = seed["pfam_sequence_accession"].rsplit(".", 1)[-1]
        current_version = metadata["Sequence version"]
        fasta_parts.append(
            f">{accession} pfam_id={seed['pfam_sequence_id']} "
            f"subgroup={selection['subgroup'].replace(' ', '_')}\n"
            f"{wrap_fasta(sequence)}\n"
        )
        metadata_rows.append(
            {
                "selection_index": index,
                "uniprot_accession": accession,
                "uniprot_entry_name": metadata["Entry Name"],
                "reviewed": metadata["Reviewed"],
                "protein_name": metadata["Protein names"],
                "organism": metadata["Organism"],
                "taxonomy_id": metadata["Organism (ID)"],
                "full_protein_length": metadata["Length"],
                "uniprot_sequence_version": current_version,
                "pfam_accession": PFAM_ACCESSION,
                "pfam_version": pfam_version,
                "pfam_sequence_id": seed["pfam_sequence_id"],
                "pfam_sequence_accession": seed["pfam_sequence_accession"],
                "domain_start": seed["domain_start"],
                "domain_end": seed["domain_end"],
                "domain_length": len(sequence),
                "seed_and_current_sequence_version_match": seed_version
                == current_version,
                "subgroup": selection["subgroup"],
                "selection_reason": selection["selection_reason"],
                "sequence_sha256": sha256_bytes(sequence.encode("ascii")),
            }
        )

    fasta_path = processed / "PF00042_globin_pilot.fasta"
    metadata_path = processed / "PF00042_globin_pilot_metadata.tsv"
    manifest_path = processed / "PF00042_globin_pilot_manifest.json"
    fasta_path.write_text("".join(fasta_parts))
    with metadata_path.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(
            destination, fieldnames=list(metadata_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(metadata_rows)

    entry_headers = entry_record.get("response_headers", {})
    uniprot_headers = uniprot_record.get("response_headers", {})
    manifest = {
        "dataset": "PF00042 globin pilot",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "selection_rule": (
            "Select three distinct reviewed seed members from each of the "
            "myoglobin, beta-like haemoglobin, and alpha-like haemoglobin "
            "groups, plus three reviewed divergent seed members representing "
            "plant, bacterial, and nematode globins. Require standard amino "
            "acids, unique domain sequences, exact Pfam domain boundaries, and "
            "length at most 511 residues."
        ),
        "pfam_accession": PFAM_ACCESSION,
        "pfam_version": pfam_version,
        "interpro_release": entry_headers.get("interpro-version"),
        "interpro_entry_name": entry["metadata"]["name"]["name"],
        "interpro_seed_count": entry["metadata"]["entry_annotations"][
            "alignment:seed"
        ],
        "uniprot_release": uniprot_headers.get("x-uniprot-release"),
        "uniprot_release_date": uniprot_headers.get("x-uniprot-release-date"),
        "licence": "CC0 1.0 Universal for InterPro/Pfam downloadable data",
        "selected_sequence_count": len(metadata_rows),
        "sources": {
            "interpro_entry": INTERPRO_ENTRY_URL,
            "pfam_seed_alignment": INTERPRO_SEED_URL,
            "uniprot_metadata": metadata_url,
        },
        "files": {
            relative(external / "PF00042.interpro-entry.json", root): sha256_bytes(
                entry_content
            ),
            relative(external / "PF00042.seed.sto.gz", root): sha256_bytes(
                seed_content
            ),
            relative(external / "PF00042.selected-uniprot.tsv", root): sha256_bytes(
                uniprot_content
            ),
            relative(fasta_path, root): sha256_bytes(fasta_path.read_bytes()),
            relative(metadata_path, root): sha256_bytes(metadata_path.read_bytes()),
        },
        "generator": "scripts/prepare_pfam_globin_pilot.py",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"Prepared {len(metadata_rows)} sequences from {pfam_version}")
    print(f"FASTA: {fasta_path}")
    print(f"Metadata: {metadata_path}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
