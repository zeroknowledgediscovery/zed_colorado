# General ZeBRA + genomics comparison pipeline

This directory is the disease-agnostic counterpart of `../zebra_comp/`. The original `zebra_comp/` remains the frozen Colorado ILD/IPF manuscript analysis. This package makes the disease-specific pieces declarative so the same analysis sequence can be rerun for ILD/IPF, ADRD, HFrEF, or another disease/gene pair without editing analysis code.

**Parity status.** This package reproduces the core scientific analysis logic, but it is not yet a byte-for-byte replacement for the historical notebook-32 implementation (including its original model-search/SHAP details). Until the ILD reference configuration is regression-checked against the frozen manuscript outputs, `../zebra_comp/` remains authoritative for published/manuscript numbers.

## What is configured per disease

A JSON config specifies:

1. **Genomic matrix**: file, patient ID, and which raw genomic columns are included.
2. **Phenotype**: label column and how positive/negative/missing values are encoded.
3. **ZeBRA score**: prediction file, patient ID, and score column.
4. **Gene signal**: the binary candidate-gene state used for local/LR/rescue analyses (for example MUC5B T-carrier, APOE4 carrier, or TTN pathogenic-variant carrier).

The broad genomic feature panel and the focal gene signal are deliberately separate. Global analyses can use hundreds/thousands of genomic columns while local likelihood-ratio analyses focus on one interpretable binary genomic channel.

## Run sequence

`run_all.sh` executes the same sequence for every config:

1. `01_GLOBAL_COMPARISON` — repeated held-out ZeBRA vs genomics-only vs combined comparison plus stringent operating points.
2. `02_INCREMENTAL_LOGISTIC` — paired regularized ZeBRA-only vs ZeBRA+genomics incremental analysis.
3. `03_GENE_SCORE_ASSOCIATION` — pooled and disease-stratified ZeBRA-to-gene association/enrichment.
4. `04_LOCAL_INFORMATION` — local nested-model LR/deviance statistics and local effect sizes along the ZeBRA axis.
5. `05_DIRECT_LR` — cross-fitted score-level `Lambda_Z` and `Lambda_{G|Z}`, finite-stage `b_delta`, clinical-tail summaries, and decision-flip localization.
6. `06_BOUNDARY_RESCUE` — leakage-free gene rescue between lower/upper ZeBRA thresholds compared with a ZeBRA-only threshold matched on training-set FPR.

Each run writes a frozen config/parameter/version manifest and the encoded feature list under `RESULTS/<analysis_name>/`.

## Basic use

```bash
cd zebra_comp_general
pip install -r requirements.txt
./run_all.sh configs/ild_muc5b.json fast
./run_all.sh configs/ild_muc5b.json full
```

The included ILD/MUC5B config points back to the current `../zebra_comp/` processed files and is intended as the regression/reference configuration for the generalized code.

For a new disease, copy one of the templates under `configs/`, set the three input files, phenotype column, score column, genomic selector, and focal gene signal, then run the same command.

## Genomic column selection

`genomic_features` accepts any of the following selectors:

### All genomic columns

```json
{"type": "all", "exclude_columns": ["target"]}
```

### Explicit columns

```json
{"type": "columns", "columns": ["rs123_A", "rs456_G"]}
```

### A contiguous section of a wide genomic table

```json
{"type": "slice", "start_column": "FIRST_GENOMIC_COLUMN", "end_column": "LAST_GENOMIC_COLUMN"}
```

### Regex/prefix selection

```json
{"type": "regex", "pattern": "^(rs|JHU_)"}
```

### A discovered disease-gene column file

```json
{"type": "file", "path": "../ADRD_genomic_header_matches.csv", "column_field": "header_name"}
```

Multiple selectors can be unioned with `"include": [...]`.

Selected categorical genotype columns are one-hot encoded automatically; numeric genotype/dosage columns remain numeric. Invariant encoded columns are removed. Missing numeric values are median-imputed within each training split.

## Gene signal encodings

The focal gene state must be binary for the local/direct-LR/rescue chain. Supported forms are:

- `binary_column`: already encoded 0/1 carrier/status column.
- `dosage_column`: positive if dosage is at or above a threshold.
- `genotype_column`: raw categorical genotype with configured positive values.
- `onehot_any`: current MUC5B-style one-hot state encoding; positive if any configured risk-state column is active.

The ILD config demonstrates `onehot_any`. ADRD and HFrEF templates assume a precomputed `APOE4_carrier` or `TTN_pathogenic_carrier`; change those definitions to the representation actually present in the new cohort.

## Disease gene-column discovery

The existing Colorado lineage already used IPF, ADRD, and HFrEF gene panels. This generalized package moves those panels into `gene_panels.json` and provides reusable discovery scripts.

### ADRD

```bash
scripts/get_adrd_gene_columns.sh genomicdataheader.csv
```

writes `ADRD_genomic_header_matches.csv`.

### HFrEF

```bash
scripts/get_hfref_gene_columns.sh genomicdataheader.csv
```

writes `HFREF_genomic_header_matches.csv`.

Or call the generic tool:

```bash
python scripts/find_phenotype_genomic_headers.py genomicdataheader.csv --panel ADRD
python scripts/find_phenotype_genomic_headers.py genomicdataheader.csv --panel HFREF
python scripts/find_phenotype_genomic_headers.py genomicdataheader.csv --panel IPF
```

The utility resolves GRCh37 gene intervals and overlapping rsIDs through the Ensembl GRCh37 REST API and caches the result under `.ensembl_cache/`. Output has exactly `gene_name,header_name`, so it can be fed directly into a config using the `file` selector.

The ADRD panel includes APOE, APP, PSEN1/2, TREM2, SORL1, ABCA7, BIN1, CR1, CLU, PICALM, CD33, INPP5D, PLCG2 and additional dementia genes. The HFrEF panel is centered on cardiomyopathy/heart-failure genes including TTN, MYBPC3, BAG3, FLNC, LMNA, RBM20, DSP, PLN, MYH7, TNNT2/TNNI3/TNNC1, SCN5A, ACTN2 and FHOD3.

## ADRD example workflow

```bash
# 1. Find all ADRD-panel columns represented in the raw genomic header.
scripts/get_adrd_gene_columns.sh /path/to/genomicdataheader.csv ADRD_genomic_header_matches.csv

# 2. Copy and edit the template.
cp configs/adrd_apoe.template.json configs/adrd_apoe.json
# Fill input paths, ADRD label column, and actual APOE4 carrier representation.

# 3. Validate quickly, then run fully.
./run_all.sh configs/adrd_apoe.json fast
./run_all.sh configs/adrd_apoe.json full
```

## HFrEF example workflow

```bash
scripts/get_hfref_gene_columns.sh /path/to/genomicdataheader.csv HFREF_genomic_header_matches.csv
cp configs/hfref_ttn.template.json configs/hfref_ttn.json
# Fill input paths, HFrEF label, and actual TTN carrier representation.
./run_all.sh configs/hfref_ttn.json fast
./run_all.sh configs/hfref_ttn.json full
```

## Important interpretation constraint

The direct LR analysis is exact for the observed score-level factorization `Lambda_{Z,G}=Lambda_Z Lambda_{G|Z}`. It does not assume that ZeBRA is sufficient for full clinical history, and a finite observed `b_delta` must not be described as a universal genome-wide or sequential bound. The generalized package preserves the same distinction as the ILD manuscript analysis.
