import pandas as pd
import pytest
from Bio import Align
from Bio.Align import substitution_matrices

from protein_alignment_networks.graphs import score_matrix_to_graph
from protein_alignment_networks.pipeline import pairwise_score_matrix


@pytest.fixture
def protein_aligner():
    return Align.PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,
        extend_gap_score=-1,
    )


def test_pairwise_matrix_is_labelled_and_symmetric(protein_aligner):
    sequences = {"a": "MKT", "b": "MNT", "c": "GGG"}
    matrix = pairwise_score_matrix(sequences, aligner=protein_aligner)
    assert list(matrix.index) == ["a", "b", "c"]
    pd.testing.assert_frame_equal(matrix, matrix.T.rename_axis("protein_id"))
    assert matrix.loc["a", "b"] == matrix.loc["b", "a"]
    assert matrix.loc["a", "b"] == protein_aligner.score("MKT", "MNT")


def test_threshold_graph_keeps_isolates_and_score_attributes():
    matrix = pd.DataFrame(
        [[6.0, 4.0, 0.0], [4.0, 6.0, 1.0], [0.0, 1.0, 6.0]],
        index=["a", "b", "c"],
        columns=["a", "b", "c"],
    )
    graph = score_matrix_to_graph(matrix, threshold=4.0)
    assert set(graph.nodes) == {"a", "b", "c"}
    assert set(graph.edges) == {("a", "b")}
    assert graph["a"]["b"]["score"] == 4.0


def test_graph_rejects_asymmetric_matrix():
    matrix = pd.DataFrame([[1.0, 2.0], [0.0, 1.0]], index=["a", "b"], columns=["a", "b"])
    with pytest.raises(ValueError, match="symmetric"):
        score_matrix_to_graph(matrix, threshold=1.0)
