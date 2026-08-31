# Notebooks

Notebooks are numbered in intended reading order:

- `01_goal_1_pipeline.ipynb`: small end-to-end pipeline using Biopython local
  alignment with BLOSUM62 and affine gap penalties.
- `02_pf00042_complete_pipeline.ipynb`: complete PF00042 pipeline from protein
  set through selectable methods and configurable checks to separate saved
  method graphs.
- `03_large_scale_graph_exploration.ipynb`: first a 205-domain controlled Pfam
  panel, then the 960-domain, 24-family repeated-sampling experiment, then the
  complete 21,600-graph design. It covers collection composition, graph rules,
  sparsity, visualization, unsupervised community recovery, family/clan NMI
  distributions, the score-percentile density artefact, matched edge budgets,
  the collection-uniqueness audit, and the relative size of each experimental
  factor. Sections 13 onwards read tables written by
  `scripts/summarize_resampled_graph_experiment.py`, so run that script first.
  DEDAL is optional and disabled by default.

The project uses Biopython directly rather than maintaining its own alignment
algorithm. See `docs/BIOPYTHON_ALIGNMENT.md` for the minimal configuration.

Keep notebooks reproducible: state their purpose, use relative paths, fix random
seeds, avoid hidden state, and restart the kernel before running all cells.
Move reusable code into the package under `src/` and test it there.
