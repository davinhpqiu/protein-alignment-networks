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

## Repository map

```text
data/                         Local data, separated by processing stage
docs/                         Public method setup and usage guides
notebooks/                    Numbered, reproducible analyses
outputs/                      Generated figures, tables, and graph files
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
