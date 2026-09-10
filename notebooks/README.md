# Notebooks

The notebooks form one narrative, but they have different roles:

1. `01_goal_1_pipeline.ipynb` is a tiny teaching example. It explains local
   alignment, an all-pairs score matrix, and the score-to-graph step without
   making biological claims.
2. `02_pf00042_complete_pipeline.ipynb` is the complete engineering example.
   It loads a small curated Pfam globin-domain pilot, combines Biopython,
   BLAST, and optional DEDAL scores, writes separate method graphs, and shows
   exactly where the stored pair, graph, quality, runtime, and provenance
   results live.
3. `03_large_scale_graph_exploration.ipynb` is the main scientific analysis
   and can be read independently. Part A develops the analysis on a controlled
   panel; Part B loads the repeated-sampling experiment, audits the collections,
   compares graph rules and methods, and derives the report-facing figures and
   conclusions.

Empirical results must not be copied into Markdown. Code cells load the
authoritative TSV/JSON artifacts under `data/processed/` and `outputs/tables/`,
perform the analysis, and display the result. Markdown explains the question,
method, interpretation, limitation, and path to the source artifact. Each
notebook ends with a provenance map; Notebook 03 also derives a compact
`headline_results` table from the loaded files.

Reusable implementation belongs in `src/protein_alignment_networks/` and is
tested under `tests/`. Scripts under `scripts/` orchestrate complete recorded
runs and create the artifacts consumed by the notebooks. In particular:

- `run_pair_scores.py` creates the canonical all-pairs score table;
- `run_resampled_graph_experiment.py` creates raw graph/community evaluations;
- `run_score_permutation_null.py` creates the deterministic negative-control
  evaluations;
- `summarize_resampled_graph_experiment.py` audits samples and creates summary
  tables and figures.

Run notebooks from the repository root or the `notebooks/` directory using the
**Protein Alignment Networks** kernel. Restart the kernel and run all cells in
order; there should be no hidden state. Notebook 03 expects the completed
large-panel tables described in the root `README.md`. DEDAL is optional and is
not regenerated during the scaled experiment.

Scientific citations appear inline and their canonical records are stored in
`references/references.bib`. The project uses maintained standard alignment
implementations rather than maintaining its own Smith--Waterman or BLAST code.
