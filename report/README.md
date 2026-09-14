# Report

`report.tex` adapts supplied Oxford thesis styling to a single-column research
report. Thesis front matter and chapter structure are omitted.

Empirical values in prose come from `generated_results.tex`, created by
`scripts/build_report_assets.py` from saved TSV and JSON artifacts. Notebook 03,
Section B8 runs same script. It writes canonical figures to
`outputs/figures/pfam_large_panel_report/`, copies identical files to
`report/figures/`, and records source tables and SHA-256 checksums in
`report/figure_provenance.tsv`. It also copies static image assets and
`references.bib`.

Python dependencies are listed in root `requirements.txt`. LaTeX requires
`latexmk` plus packages imported at top of `report.tex`; a full TeX Live
installation includes them.

Generate assets from project root:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/protein-report-mpl \
  python scripts/build_report_assets.py
```

Compile from `report/`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=protein_alignment_networks_report \
  -output-directory=../output/pdf report.tex
```

Output: `output/pdf/protein_alignment_networks_report.pdf`.
