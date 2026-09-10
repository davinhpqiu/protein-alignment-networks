# Data guide

Data files are intentionally not committed. For every dataset, record its name,
source URL or accession, version, download date, licence, filters, checksums if
available, and the script/notebook that produced each derivative.

- `raw/`: original data collected directly for this project; never edit in place.
- `external/`: third-party datasets retained exactly as received.
- `interim/`: temporary cleaned or transformed data.
- `processed/`: analysis-ready data produced by a documented pipeline.

Suggested dataset record:

```text
Dataset:
Source/accession:
Version and download date:
Licence:
Original filename:
SHA-256:
Filters or inclusion criteria:
Produced derivatives:
```

## PF00042 globin pilot

- **Dataset:** PF00042 globin pilot (12 protein-domain sequences)
- **Source/accession:** Pfam `PF00042` (Globin), obtained from the official
  InterPro API seed alignment; reviewed-entry metadata obtained from UniProtKB
- **Versions and download date:** Pfam `PF00042.29`, InterPro `109.0`, UniProt
  `2026_02`; downloaded 2026-08-20
- **Licence:** InterPro/Pfam downloadable data are CC0 1.0 Universal
- **Selection:** Three myoglobins, three beta-like haemoglobins, three
  alpha-like haemoglobins, and three divergent globins (plant, bacterial, and
  nematode) from the curated Pfam seed. All are reviewed UniProt records with
  unique, standard-amino-acid domain sequences no longer than 511 residues.
- **Input representation:** Ungapped Pfam domain sequences, not complete
  proteins. The seed alignment is retained only as the immutable source and is
  not supplied to any comparison method.
- **Generator:** `scripts/prepare_pfam_globin_pilot.py`
- **Raw sources:** `data/external/pfam_pf00042/`
- **Analysis-ready files:**
  `data/processed/pfam_pf00042_globin_pilot/PF00042_globin_pilot.fasta`,
  `PF00042_globin_pilot_metadata.tsv`, and
  `PF00042_globin_pilot_manifest.json`
- **Checksums:** Recorded in the generated manifest. The prepared FASTA SHA-256
  is `e01429d29d2b5aebe71b188bc4a72117bd6cd4daf8a741a5dbe3acc11514cb9e`.

Reproduce the dataset from the project root with:

```bash
python scripts/prepare_pfam_globin_pilot.py
```

Add `--refresh` only when intentionally checking for a newer upstream release;
the resulting version and checksums may then change.

## Pfam clan/family exploration panel

- **Dataset:** Pfam clan/family exploration panel (205 protein-domain
  sequences; 20,910 unique pairs)
- **Families:** `PF00042` and `PF01152` from globin clan `CL0090`; `PF00069`
  and `PF07714` from kinase clan `CL0016`; `PF00018` and `PF08239` from SH3
  clan `CL0010`
- **Versions and download date:** Pfam versions recorded per family in the
  manifest, InterPro `109.0`, UniProt `2026_02`; downloaded 2026-08-26
- **Selection:** Up to 40 valid unique domain sequences per curated Pfam seed,
  selected deterministically by farthest-point sampling on seed-alignment
  identity. Seeds with fewer than 40 eligible sequences contribute all of
  them.
- **Use of Pfam alignment:** Selection only. The scoring methods receive the
  ungapped domain sequences and never receive the seed alignment.
- **Generator:** `scripts/prepare_pfam_clan_panel.py`
- **Raw sources:** `data/external/pfam_clan_panel/`
- **Analysis-ready files:** `data/processed/pfam_clan_panel/`
- **Checksums and exact counts:**
  `data/processed/pfam_clan_panel/pfam_clan_panel_manifest.json`

Reproduce the collection from the project root with:

```bash
python scripts/prepare_pfam_clan_panel.py
```

## Reproducibly sampled large Pfam panel

- **Dataset:** Large Pfam family/clan panel (960 seed-domain instances; 460,320
  unique master pairs)
- **Release:** Pfam `38.2`, based on UniProtKB `2025_03`; selected and downloaded
  2026-08-29
- **Family sampling frame:** Official current-release
  `Pfam-A.clans.tsv.gz`. Unassigned families were excluded. With fixed seed
  `20260828`, clans and their member families were shuffled; the first eight
  clans containing three domain/family entries with at least 40 reported seed
  domains were retained.
- **Selected panel:** Eight clans, three families per clan, and 40 valid,
  diversity-selected seed-domain instances per family. Exact accessions are in
  `data/processed/pfam_sampling_frame/selected_families.tsv`; all 415 inspected
  candidates and rejections are in `family_selection_trace.tsv`.
- **Repeated collections:** A factorial balanced design with 2, 4, or 8 clans;
  2 or 3 families per clan; 10, 20, or 40 domain instances per family; and 50
  replicates per condition. This gives 900 predeclared collections. Sampling is
  without replacement within each collection.
- **Generators:** `scripts/select_pfam_family_panel.py`,
  `scripts/prepare_pfam_clan_panel.py`, and
  `scripts/generate_pfam_collections.py`
- **Raw sources:** `data/external/pfam_sampling_frame/` and
  `data/external/pfam_large_panel/`
- **Analysis-ready files:** `data/processed/pfam_sampling_frame/` and
  `data/processed/pfam_large_panel/`
- **Provenance:** Source and output checksums, release information, parameters,
  fixed seeds, selected membership, and selection traces are stored in the
  generated manifests and TSV tables.

The first diagnostic slice fixes each collection at four clans, two families
per clan, and 20 domains per family. The complete factorial design has since
been evaluated. Notebook 03 loads the manifests and raw result tables to derive
the realised collection, graph, eligibility, and shortfall counts.
