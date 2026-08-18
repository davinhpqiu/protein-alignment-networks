"""Input, output, and validation helpers."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

AMINO_ACID_CODES = frozenset("ACDEFGHIKLMNPQRSTVWY")


def validate_protein_sequence(sequence: str, *, allow_empty: bool = True) -> str:
    """Normalise and validate a sequence of the 20 standard amino acids."""

    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    normalised = "".join(sequence.split()).upper()
    if not normalised and not allow_empty:
        raise ValueError("protein sequence cannot be empty")
    invalid = sorted(set(normalised) - AMINO_ACID_CODES)
    if invalid:
        raise ValueError(f"invalid amino-acid code(s): {', '.join(invalid)}")
    return normalised


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a simple FASTA file into an insertion-ordered ID-to-sequence map.

    The first whitespace-delimited field of each header is used as the unique
    identifier. Sequences are validated against the 20 standard amino acids.
    """

    records: dict[str, list[str]] = {}
    current_id: str | None = None
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                current_id = line[1:].split(maxsplit=1)[0]
                if not current_id:
                    raise ValueError(f"empty FASTA header at line {line_number}")
                if current_id in records:
                    raise ValueError(f"duplicate FASTA identifier: {current_id}")
                records[current_id] = []
            elif current_id is None:
                raise ValueError(
                    f"sequence data before the first FASTA header at line {line_number}"
                )
            else:
                records[current_id].append(line)
    if not records:
        raise ValueError("FASTA file contains no records")
    return {
        identifier: validate_protein_sequence("".join(parts), allow_empty=False)
        for identifier, parts in records.items()
    }


def write_fasta(records: dict[str, str], path: str | Path, line_width: int = 80) -> None:
    """Write validated sequences to FASTA format."""

    if line_width < 1:
        raise ValueError("line_width must be positive")
    lines: list[str] = []
    for identifier, sequence in records.items():
        if not identifier or any(character.isspace() for character in identifier):
            raise ValueError("FASTA identifiers must be non-empty and contain no spaces")
        clean_sequence = validate_protein_sequence(sequence, allow_empty=False)
        lines.append(f">{identifier}")
        lines.extend(_chunks(clean_sequence, line_width))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _chunks(sequence: str, width: int) -> Iterable[str]:
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]

