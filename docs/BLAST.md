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
```

`hits` preserves every high-scoring segment pair reported by BLAST. `pairs`
removes self-hits, treats the protein relationship as undirected, and keeps the
highest-bit-score hit for each unique protein pair.

The runner explicitly records standard BLASTP scoring settings: BLOSUM62,
gap-opening cost 11, gap-extension cost 1, and composition-based statistics
mode 2. The default reporting threshold is E-value 10. A protein pair absent
from the table was not reported by BLAST at that threshold; absence is not a
raw score of zero.

## Important comparison detail

BLAST's documented gap convention charges the existence cost and an extension
cost for every residue in a gap. Biopython's `PairwiseAligner` charges its open
score for the first residue and its extension score for each later residue.
Therefore, Biopython `open_gap_score=-12` and `extend_gap_score=-1` reproduce
the numerical gap cost described by BLAST as 11 and 1.
