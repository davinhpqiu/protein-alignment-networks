# References

This folder separates formal citation metadata from research notes:

- `references.bib` is the canonical BibTeX database used by the report and any
  citation-aware Markdown workflow.
- `LITERATURE_TRACKER.md` records local reading status, relevance, limitations,
  and project-specific notes.
- `GRAPH_CONSTRUCTION_LITERATURE.md` is a local working review of how protein
  sequence-similarity graphs are constructed.

`references.bib` and this guide are published with the repository so notebook
citations resolve for other readers. The two working-note files remain local.

## Workflow for a new paper

1. Add or import a verified BibTeX record into `references.bib`.
2. Give it a stable key in the form `firstauthorYYYYkeyword`.
3. Add one row to the reading table.
4. When reading, add structured notes using the supplied template.
5. Cite it in Markdown as `[@citation_key]` if using Pandoc or Quarto.
6. Change the tracker status to **Cited** only when it appears in a deliverable.

Check author names, title, venue, year, pages, DOI, and URL against the
publisher or another authoritative record. Do not store copyrighted paper PDFs
in Git unless their licence clearly permits redistribution.

The original proposal remains at the repository root. Part of it is explicitly
confidential, so it is ignored by Git and must not be redistributed.
