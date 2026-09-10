# Running BLAST locally

The project uses NCBI BLAST+ directly. It does not reimplement BLAST and does
not send sequences to the NCBI web service.

## Installation

On macOS with Homebrew:

```bash
brew install blast
blastp -version
```

BLAST+ is a system program rather than a Python package, so it is intentionally
not listed in `requirements.txt`.

Homebrew's current BLAST+ 2.17.0 bottle on this Mac prints an external MBEDTLS
header/runtime version warning. Database creation, real searches, parsing, and
the automated smoke test all complete successfully despite that packaging
warning.

## Python usage

```python
from protein_alignment_networks import (
    blastp_all_vs_all,
    blastp_version,
    summarize_blast_directionality,
    summarize_blast_pairs,
)

sequences = {
    "protein_a": "PAWHEAE",
    "protein_b": "PAWHDQE",
    "protein_c": "MKTAYIAKQRQISFVKSHFSRQ",
}

print(blastp_version())
hits = blastp_all_vs_all(
    sequences,
    output_path="outputs/tables/blast_toy_hits.tsv",
    threads=4,
)
pairs = summarize_blast_pairs(hits)
directionality = summarize_blast_directionality(hits)
```

`hits` preserves every high-scoring segment pair reported by BLAST. `pairs`
removes self-hits, treats the protein relationship as undirected, and keeps the
highest-bit-score hit for each unique protein pair. `directionality` records how
many unordered pairs appeared in both query directions or only one, plus the
defined relative bit-score spread for bidirectional pairs.

The bit score used as the graph weight is
$S'=(\lambda S-\ln K)/\ln 2$: raw alignment score $S$ is rescaled by the
scoring system's Karlin--Altschul statistical parameters $\lambda$ and $K$.
This definition follows Karlin and Altschul (1990) and Altschul et al. (1990),
recorded as `karlin1990methods` and `altschul1990basic` in
`references/references.bib`.

The runner explicitly records standard BLASTP scoring settings: BLOSUM62,
gap-opening cost 11, gap-extension cost 1, and composition-based statistics
mode 2. The default reporting threshold is E-value 10. A protein pair absent
from the table was not reported by BLAST at that threshold; absence is not a
raw score of zero. The local runner sets `-max_target_seqs` to the full input
sequence count, preventing the program's usual target limit from silently
truncating a larger all-versus-all search. The run manifest records that cap
and the directionality audit derived from the saved HSP table.

## Important comparison detail

BLAST's documented gap convention charges the existence cost and an extension
cost for every residue in a gap. Biopython's `PairwiseAligner` charges its open
score for the first residue and its extension score for each later residue.
Therefore, Biopython `open_gap_score=-12` and `extend_gap_score=-1` reproduce
the numerical gap cost described by BLAST as 11 and 1. The project's recorded
Biopython baseline instead uses `-11` and `-1`, so the two methods share the
standard BLOSUM62 parameter names but are not presented as numerically
identical scoring systems.
