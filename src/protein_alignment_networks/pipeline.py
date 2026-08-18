"""Composable steps for the sequence-to-network pipeline."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
from Bio.Align import PairwiseAligner

from .io import validate_protein_sequence


def pairwise_score_matrix(
    sequences: Mapping[str, str],
    *,
    aligner: PairwiseAligner,
) -> pd.DataFrame:
    """Use a configured Biopython aligner to score every sequence pair."""

    if not sequences:
        raise ValueError("at least one sequence is required")
    names = list(sequences)
    if len(names) != len(set(names)):
        raise ValueError("sequence identifiers must be unique")
    clean_sequences = {
        name: validate_protein_sequence(sequence, allow_empty=False)
        for name, sequence in sequences.items()
    }
    matrix = pd.DataFrame(0.0, index=names, columns=names)
    for i, name_a in enumerate(names):
        for j in range(i, len(names)):
            name_b = names[j]
            score = float(
                aligner.score(clean_sequences[name_a], clean_sequences[name_b])
            )
            matrix.iat[i, j] = score
            matrix.iat[j, i] = score
    matrix.index.name = "protein_id"
    matrix.columns.name = "protein_id"
    return matrix
