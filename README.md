````markdown
# rnacentral-parser

Lightweight utilities to parse RNAcentral FASTA releases and construct sequence-level and TaxID-aware mappings without relying on the official PostgreSQL dump.

## Rationale

The recommended PostgreSQL setup from the RNAcentral FTP archive is currently fragile:  
the provided dumps reference schemas that are not defined in the archive itself and cannot be reliably reconstructed from external sources.

The most robust approach is therefore to **read the raw FTP files directly** and construct the required mappings explicitly (“hard-coded” from the authoritative flat files). This repository follows that strategy.

## Data source

RNAcentral FTP archive:  
https://rnacentral.org/help/ftp

The code assumes you are working from the **raw `current_release` directory**, without any intermediate database layer.

## Environment setup

Create and activate a clean python environment.
````

Install required Python dependencies:

```bash
pip install pandas pyarrow biopython duckdb
```

## Required FTP directories

From the RNAcentral FTP `current_release`, download the following directories **as-is**:

```
go_annotations/
id_mapping/
md5/
rfam/
sequences/
```

Example local layout:

```
rnacentral_data/
└── current_release/
    ├── go_annotations/
    ├── id_mapping/
    ├── md5/
    ├── rfam/
    └── sequences/
```

In particular:

* `sequences/` contains the compressed FASTA files (URS → sequence)
* `id_mapping/` provides URS ↔ external IDs ↔ TaxID relationships
* `rfam/` and `go_annotations/` support functional annotation joins

## Conceptual model

* **URS IDs are sequence-centric**: one URS corresponds to exactly one RNA sequence.
* **Taxonomic assignment is external**: the same URS can map to multiple TaxIDs.
* Species information is **not inferred from FASTA headers** and must be joined explicitly.

This parser enforces that separation to avoid silent biological assumptions.

## Scope

This repository is intended for:

* FASTA parsing (URS, RNA class, sequence)
* Construction of explicit URS–TaxID mappings
* Integration with interaction or annotation datasets
* Scalable downstream processing (DuckDB / Arrow)

It is **not** a reimplementation of the RNAcentral database.

## Notes

If you need organism-specific analyses, always model data as:

```
(URS sequence) × (TaxID context)
```

Never assume URS implies a unique organism.
