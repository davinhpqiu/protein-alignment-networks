# Protein Sequence Alignment Networks

A reproducible research project for comparing networks built from protein
sequence-alignment scores. The immediate focus is Goal 1 from the project
proposal:

```text
protein sequences -> pairwise alignment-score matrix -> graph
```

Pairwise local alignment is delegated directly to Biopython, using BLOSUM62 and
affine gap penalties. NCBI BLAST+ and the pretrained DEDAL model are integrated
as additional local methods; DEDAL uses a separate optional environment.

## Quick start

Python 3.11 or 3.12 is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e .
brew install blast
python -m ipykernel install --user --name protein-alignment-networks \
  --display-name "Protein Alignment Networks"
jupyter lab
```

Open `notebooks/01_goal_1_pipeline.ipynb` and select the **Protein Alignment
Networks** kernel. Run all cells from top to bottom.

### Running the notebook in VS Code

The project environment has already been created locally as `.venv`.

1. Open the complete `Protein Seq Alignment` folder in VS Code, not only the
   notebook file.
2. Open `notebooks/01_goal_1_pipeline.ipynb`.
3. Click **Select Kernel** in the top-right corner.
4. Choose **Python Environments...**, as shown in the kernel menu.
5. Select the interpreter whose path ends with
   `Protein Seq Alignment/.venv/bin/python`. It may be labelled simply
   `.venv (Python 3.12)`.
6. Click **Run All** near the top of the notebook, or run one cell at a time with
   the triangular button beside each cell.

If `.venv` is not listed, choose **Enter interpreter path...** and enter:

```text
<project folder>/.venv/bin/python
```

A successful run prints a local-alignment score of `12`, followed later by
`Nodes: 4` and `Edges: 3`. If VS Code asks to install the Python or Jupyter
extension, install the official Microsoft extension it names and repeat the
kernel selection.

### Running the notebook in a browser

From a Terminal opened in the project folder:

```bash
source .venv/bin/activate
jupyter lab
```

Open the notebook from Jupyter's file browser, select the project kernel, and
use **Run > Run All Cells**. Stop Jupyter later by returning to the Terminal and
pressing `Control-C` twice.

Run the checks with:

```bash
python -m pytest
```

See `docs/BLAST.md` and `docs/DEDAL.md` for the optional local methods and their
recorded outputs.

## Documentation

- `docs/PROJECT.md`: the research question, scope, accepted decisions, current
  evidence, progress, and next actions. This is the only project-status record.
- `docs/BIOPYTHON_ALIGNMENT.md`, `docs/BLAST.md`, and `docs/DEDAL.md`:
  method-specific setup and parameter guides.
- `data/README.md`, `notebooks/README.md`, and `outputs/README.md`: concise
  guides for those folders.
- `references/`: BibTeX records and the literature tracker.

The original proposal PDF and `Notes.docx` are retained unchanged as source
documents. Local research records are intentionally not published to GitHub.

## Prepare the pilot dataset

The first real-data experiment uses 12 reviewed protein-domain sequences from
the curated Pfam `PF00042` globin seed alignment. It deliberately includes
three myoglobins, three alpha-like haemoglobins, three beta-like haemoglobins,
and three more divergent globins so that the three methods have both easier and
harder relationships to compare.

With the main environment active, run:

```bash
python scripts/prepare_pfam_globin_pilot.py
```

The first run downloads official Pfam/InterPro and UniProt records; later runs
reuse those exact recorded sources. The generated FASTA, metadata table, and
provenance manifest are stored under
`data/processed/pfam_pf00042_globin_pilot/` and are intentionally ignored by
Git. See `data/README.md` for the exact selection and recorded versions.

## Run the complete PF00042 pipeline

The second notebook runs the protein set through Biopython, BLAST, and DEDAL,
records configurable correlation checks, and writes a separate graph for each
method under explicit graph settings.

Generate the quick Biopython and BLAST results with:

```bash
python scripts/run_pf00042_correlation.py --methods biopython blast --threads 4
```

Add the slower local DEDAL results when its optional environment is ready:

```bash
python scripts/run_pf00042_correlation.py --methods dedal
```

Completed outputs are reused, so the two commands can be run separately. Then
open `notebooks/02_pf00042_complete_pipeline.ipynb` and run all cells. The
method, graph threshold rule, graph threshold, and correlation-check thresholds
are visible configuration values near the top. The included `top_n=13` graphs
only verify the pipeline; scientific settings remain unchosen until data
exploration.

Notebook 02 remains the small end-to-end pipeline check. Its `top_n=13` rule
means retaining exactly the 13 strongest available pairs and should not be used
as a scientific default.

## Run the large-scale graph exploration

Open `notebooks/03_large_scale_graph_exploration.ipynb` and run all cells. Part
A retains the 205-domain, six-family controlled panel used to establish the
scaled graph workflow. Part B loads the main repeated-sampling experiment: 960
Pfam seed-domain instances from 24 reproducibly selected families in eight
clans, 900 balanced sampled collections, and NMI/AMI recovery results after
unsupervised Louvain community detection. Part C, sections 13 to 17, loads the
complete design of 21,600 graphs, and covers the score-percentile density
artefact, matched edge budgets, the collection-uniqueness audit, and the
relative size of each experimental factor. DEDAL is optional and remains
disabled by default.

The collection can also be prepared independently with:

```bash
python scripts/prepare_pfam_clan_panel.py
```

Reproduce the larger sampling experiment with:

```bash
python scripts/select_pfam_family_panel.py \
  --clans 8 --families-per-clan 3 \
  --min-reported-seed-domains 40 --seed 20260828
python scripts/prepare_pfam_clan_panel.py \
  --family-file data/processed/pfam_sampling_frame/selected_families.tsv \
  --dataset-id pfam_large_panel --per-family 40
python scripts/generate_pfam_collections.py \
  data/processed/pfam_large_panel/pfam_large_panel_metadata.tsv \
  data/processed/pfam_large_panel/resampling \
  --clans 2 4 8 --families-per-clan 2 3 \
  --domains-per-family 10 20 40 --replicates 50 --seed 20260828
python scripts/run_pair_scores.py \
  --fasta data/processed/pfam_large_panel/pfam_large_panel.fasta \
  --dataset-manifest \
    data/processed/pfam_large_panel/pfam_large_panel_manifest.json \
  --output-directory outputs/tables/pfam_large_panel_scores \
  --methods biopython blast --threads 8
```

The first recovery slice used 50 fixed-design collections and can be
regenerated with:

```bash
python scripts/run_resampled_graph_experiment.py \
  data/processed/pfam_large_panel/resampling/collection_manifest.tsv \
  data/processed/pfam_large_panel/resampling/collection_membership.tsv \
  data/processed/pfam_large_panel/pfam_large_panel_metadata.tsv \
  outputs/tables/pfam_large_panel_scores/paired_scores.tsv \
  outputs/tables/pfam_large_panel_graph_recovery_initial \
  --methods biopython --top-k 1 2 5 10 \
  --clans 4 --families-per-clan 2 --domains-per-family 20
```

### Reproduce the complete design

The complete predeclared experiment is two runs over the same 900 collections,
so that every method comparison is paired on identical protein sets. Both are
checkpointed every 25 collections and accept `--resume`, so an interrupted run
continues rather than restarting. Together they take roughly an hour.

```bash
python scripts/run_resampled_graph_experiment.py \
  data/processed/pfam_large_panel/resampling/collection_manifest.tsv \
  data/processed/pfam_large_panel/resampling/collection_membership.tsv \
  data/processed/pfam_large_panel/pfam_large_panel_metadata.tsv \
  outputs/tables/pfam_large_panel_scores/paired_scores.tsv \
  outputs/tables/pfam_large_panel_graph_recovery_full \
  --methods biopython blast \
  --top-k 1 2 5 10 --percentiles 0.90 0.95 0.98 0.99 --resume

python scripts/run_resampled_graph_experiment.py \
  data/processed/pfam_large_panel/resampling/collection_manifest.tsv \
  data/processed/pfam_large_panel/resampling/collection_membership.tsv \
  data/processed/pfam_large_panel/pfam_large_panel_metadata.tsv \
  outputs/tables/pfam_large_panel_scores/paired_scores.tsv \
  outputs/tables/pfam_large_panel_graph_recovery_matched_density \
  --methods biopython blast --top-k \
  --target-densities 0.005 0.01 0.02 0.05 --resume
```

The bare `--top-k` in the second command clears the default `1 2 5 10`, so
that run produces target-density graphs only. The `target_density` rule gives
both methods the same edge budget on the same nodes. Equal score percentiles do not: BLAST's
no-hit pairs shrink its available-score pool, so the same nominal percentile
produces a much sparser BLAST graph. See `docs/PROJECT.md` for the effect this
has on the results.

Then derive the summary tables and figures:

```bash
MPLCONFIGDIR=/tmp/mplcache python scripts/summarize_resampled_graph_experiment.py \
  outputs/tables/pfam_large_panel_graph_recovery_full/graph_recovery_results.tsv \
  outputs/tables/pfam_large_panel_graph_recovery_matched_density/graph_recovery_results.tsv \
  data/processed/pfam_large_panel/resampling/collection_manifest.tsv \
  data/processed/pfam_large_panel/resampling/collection_membership.tsv \
  outputs/tables/pfam_large_panel_graph_summary \
  outputs/figures/pfam_large_panel_full
```

Setting `MPLCONFIGDIR` to a writable folder matters: Matplotlib can hang while
building its font cache in a protected macOS location, which silently produces
tables but no figures. Parts B and C of Notebook 03 read the tables this script
writes, so run it before the notebook.

The family selector saves every inspected candidate and rejection. The
collection generator saves every selected clan, family, domain, and random
seed. Pfam labels are used to balance the sampled input, but graph construction
and community detection do not receive them; labels are revealed only for NMI
and AMI evaluation.

In the scaled notebook, `top_k=5` means that each protein nominates its five
strongest available neighbours and the undirected graph retains the union of
those nominations. This scales with the collection and avoids the arbitrary
edge count used by the earlier smoke test.

## Repository map

```text
data/                         Local data, separated by processing stage
docs/                         Project overview and method-specific guides
notebooks/                    Numbered, reproducible analyses
outputs/                      Generated figures, tables, and graph files
scripts/                      Reproducible data-preparation and method workers
src/protein_alignment_networks/
                              Reusable, tested Python code
tests/                        Automated correctness checks
```

`data/raw/` and `data/external/` are treated as immutable inputs. Derived data
belongs in `data/interim/` or `data/processed/`; figures and result files belong
under `outputs/`. Large/generated files are ignored by Git, while the folder
guides remain tracked.

## Working conventions

- Put exploratory work in a numbered notebook, but move reusable logic into
  `src/protein_alignment_networks/` and add a test.
- Do not manually edit raw or external data.
- Record data sources and download dates in `data/README.md`.
- Record method versions, parameters, and dataset identifiers with important
  outputs.
- Use a fixed random seed in every stochastic experiment.

## Scope

Goal 1 covers alignment scores and graph construction. Goal 2 will add real and
synthetic protein collections and methods for comparing the resulting graphs.
