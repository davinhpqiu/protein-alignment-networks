"""Tools for building networks from protein sequence-alignment scores."""

from .blast import blastp_all_vs_all, blastp_version, summarize_blast_pairs
from .comparison import (
    add_descending_ranks,
    canonical_pair_index,
    merge_pair_tables,
    pair_rank_disagreements,
    score_matrix_to_pairs,
    spearman_correlations,
)
from .dedal import dedal_all_vs_all
from .graphs import (
    build_similarity_graph,
    graph_summary,
    save_graph_bundle,
    score_matrix_to_graph,
)
from .io import read_fasta, validate_protein_sequence
from .pipeline import pairwise_score_matrix
from .quality import assess_score_correlations

__all__ = [
    "add_descending_ranks",
    "assess_score_correlations",
    "blastp_all_vs_all",
    "blastp_version",
    "build_similarity_graph",
    "canonical_pair_index",
    "dedal_all_vs_all",
    "graph_summary",
    "merge_pair_tables",
    "pair_rank_disagreements",
    "pairwise_score_matrix",
    "read_fasta",
    "save_graph_bundle",
    "score_matrix_to_graph",
    "score_matrix_to_pairs",
    "spearman_correlations",
    "summarize_blast_pairs",
    "validate_protein_sequence",
]
