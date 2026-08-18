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

