# Notebooks

Notebooks are numbered in intended reading order:

- `01_goal_1_pipeline.ipynb`: small end-to-end pipeline using Biopython local
  alignment with BLOSUM62 and affine gap penalties.
- `02_pf00042_complete_pipeline.ipynb`: complete PF00042 pipeline from protein
  set through selectable methods and configurable checks to separate saved
  method graphs.
- `03_pf00042_data_exploration.ipynb`: score distributions, missing BLAST hits,
  rank correlations, and pair-level disagreements used to justify later
  quality and graph thresholds.

The project uses Biopython directly rather than maintaining its own alignment
algorithm. See `docs/BIOPYTHON_ALIGNMENT.md` for the minimal configuration.

Keep notebooks reproducible: state their purpose, use relative paths, fix random
seeds, avoid hidden state, and restart the kernel before running all cells.
Move reusable code into the package under `src/` and test it there.
