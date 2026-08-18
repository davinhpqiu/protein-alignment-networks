import pytest

from protein_alignment_networks.io import read_fasta, validate_protein_sequence


def test_sequence_validation_normalises_case_and_whitespace():
    assert validate_protein_sequence(" acd\nEFG ") == "ACDEFG"


def test_sequence_validation_rejects_nonstandard_codes():
    with pytest.raises(ValueError, match="invalid amino-acid"):
        validate_protein_sequence("ACDX")


def test_read_fasta(tmp_path):
    path = tmp_path / "example.fasta"
    path.write_text(">alpha description\nACD\nEFG\n>beta\nMKT\n", encoding="utf-8")
    assert read_fasta(path) == {"alpha": "ACDEFG", "beta": "MKT"}

