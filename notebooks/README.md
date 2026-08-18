# Notebooks

Notebooks are numbered in intended reading order:

- `01_goal_1_pipeline.ipynb`: small end-to-end pipeline using Biopython local
  alignment with BLOSUM62 and affine gap penalties.

The project uses Biopython directly rather than maintaining its own alignment
algorithm. See `docs/BIOPYTHON_ALIGNMENT.md` for the minimal configuration.

Keep notebooks reproducible: state their purpose, use relative paths, fix random
seeds, avoid hidden state, and restart the kernel before running all cells.
Move reusable code into the package under `src/` and test it there.
