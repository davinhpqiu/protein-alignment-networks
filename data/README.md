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
