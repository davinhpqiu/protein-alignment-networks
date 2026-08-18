# Running pretrained DEDAL locally

DEDAL runs locally through its published TensorFlow Hub model. No sequence is
sent to a web service. The model is downloaded once and cached under
`models/dedal/`, which is excluded from version control.

## Isolated setup

DEDAL's TensorFlow dependencies are deliberately separated from the normal
project environment:

```bash
/opt/homebrew/bin/python3.11 -m venv .venv-dedal
.venv-dedal/bin/python -m pip install -r requirements-dedal.txt
mkdir -p external models/dedal
git clone --depth 1 --filter=blob:none --sparse \
  https://github.com/google-research/google-research.git external/google-research
git -C external/google-research sparse-checkout set dedal
```

The first inference downloads `https://tfhub.dev/google/dedal/3`; subsequent
runs use the persistent cache.

The tested local configuration is Python 3.11.16 with TensorFlow macOS 2.15.0,
TensorFlow Hub 0.15.0, and TensorFlow Probability 0.23.0 on Apple Silicon.

## Python usage

Run this from the normal notebook kernel:

```python
from protein_alignment_networks import dedal_all_vs_all

results = dedal_all_vs_all(
    {
        "protein_a": "MKTAYIAKQRQISFVKSHFSRQ",
        "protein_b": "MKTAYIAKQRQISFVKSHFSRN",
    },
    output_path="outputs/tables/dedal_toy.tsv",
)
```

The main Python process launches one local worker in `.venv-dedal`, loads the
model once, aligns all requested pairs, and reads the resulting table back into
pandas. It preserves the learned Smith-Waterman score, homology logit, logistic
probability, reconstructed alignment, identity, similarity, gaps, and local
alignment coordinates.

## Limits

- Inputs may contain at most 511 residues because the published tensor has 512
  positions including the EOS token. Longer inputs are rejected, not cropped.
- DEDAL is substantially slower and more memory-intensive than BLAST or
  Biopython. On the current Mac, loading the 371 MB SavedModel took about three
  minutes and the first pair took about 32 seconds. Begin with a very small
  number of Pfam domains and measure subsequent-pair throughput before scaling.
- The logistic transformation of a homology logit is convenient, but it is not
  automatically a calibrated probability for every new dataset. Preserve the
  original logit for analysis.
