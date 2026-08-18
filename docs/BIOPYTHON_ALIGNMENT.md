# Using Biopython for local protein alignment

Biopython performs all sequence alignment in this project. We do not maintain a
separate implementation of Smith-Waterman.

Our project code begins only after the aligner has been configured: it applies
Biopython to every sequence pair, labels the resulting matrix, and converts that
matrix into a graph.

## 1. Configure Biopython

Run this in a notebook cell:

```python
from Bio import Align
from Bio.Align import substitution_matrices

aligner = Align.PairwiseAligner(
    mode="local",
    substitution_matrix=substitution_matrices.load("BLOSUM62"),
    open_gap_score=-11,
    extend_gap_score=-1,
)

print(aligner.algorithm)
```

Biopython reports a local Gotoh algorithm. Gotoh is the affine-gap extension of
Smith-Waterman: opening a gap and extending an existing gap have different
costs.

The parameters above are sensible starting values, not universal biological
truths. Record them and test sensitivity to alternatives in real experiments.

## 2. Align one pair directly

```python
sequence_a = "PAWHEAE"
sequence_b = "HEAGAWGHEE"

print("Score:", aligner.score(sequence_a, sequence_b))

alignments = aligner.align(sequence_a, sequence_b)
if len(alignments):
    print(alignments[0])
```

Use `aligner.score()` when only the numerical score is required. Use
`aligner.align()` when the actual matches and gaps need to be inspected.

## 3. Compute the complete pairwise matrix

Pass the configured Biopython object directly into the project pipeline:

```python
from protein_alignment_networks import pairwise_score_matrix

sequences = {
    "protein_a": "PAWHEAE",
    "protein_b": "HEAGAWGHEE",
    "protein_c": "MKTAYIAKQRQISFVKSHFSRQ",
}

scores = pairwise_score_matrix(sequences, aligner=aligner)
scores
```

There is no custom alignment layer between the pipeline and Biopython.

## 4. Build a graph carefully

```python
from protein_alignment_networks import score_matrix_to_graph

graph = score_matrix_to_graph(scores, threshold=20)
```

Here, `20` is only a programming example. A raw score depends on length and
composition. Before selecting graph edges scientifically, add alignment
coverage, normalised scores, and either empirical shuffled-sequence significance
or BLAST E-values. Test multiple thresholds rather than choosing one because its
graph looks convenient.

