# General ZeBRA + genomics comparison pipeline

This directory is the disease-agnostic counterpart of `../zebra_comp/`. The original `zebra_comp/` remains the frozen Colorado ILD/IPF manuscript analysis and provenance snapshot. This package makes the disease-specific pieces declarative so the same analysis sequence can be rerun for ILD/IPF, ADRD, HFrEF, or another disease/gene pair without editing analysis code.

**Validation status.** The generalized pipeline is validated against the frozen Colorado ILD/MUC5B manuscript analysis. A full reference run at code baseline `8c5db26623c3bc84ac2b0d5bae9bd109752e5dc8` produced `PASS=38, WARN=0, FAIL=0, INFO=2`. All manuscript-critical local-information, direct-LR, finite-stage `b_delta`, estimator-robustness, incremental-value, and rescue comparisons passed. The two INFO rows are the intentionally implementation-dependent global nonlinear genomics-only and combined-model AUCs. See `ILD_REFERENCE_VALIDATION.md` for the frozen validation record.

Use `../zebra_comp/` when reproducing the exact current ILD manuscript provenance. Use this generalized package for new disease/gene analyses.

## What is configured per disease

A JSON config specifies:

1. **Genomic matrix**: file, patient ID, and which raw genomic columns are included.
2. **Phenotype**: label column and how positive/negative/missing values are encoded, including whether absence from a sparse phenotype table means `negative` or `drop`.
3. **ZeBRA score**: prediction file, patient ID, and score column.
4. **Gene signal**: the binary candidate-gene state used for local/LR/rescue analyses (for example MUC5B T-carrier, APOE4 carrier, or TTN truncating-variant carrier).

The broad genomic feature panel and the focal gene signal are deliberately separate. Global analyses can use hundreds/thousands of genomic columns while local likelihood-ratio analyses focus on one interpretable binary genomic channel.

## Run sequence

`run_all.sh` executes the same sequence for every config:

1. `01_GLOBAL_COMPARISON` — repeated held-out ZeBRA vs genomics-only vs combined comparison plus stringent operating points.
2. `02_INCREMENTAL_LOGISTIC` — paired ZeBRA-only vs ZeBRA+genomics incremental analysis using the validated nested elastic-net design.
3. `03_GENE_SCORE_ASSOCIATION` — pooled and disease-stratified ZeBRA-to-gene association/enrichment.
4. `04_LOCAL_INFORMATION` — local nested-model LR/deviance statistics and local effect sizes along the ZeBRA axis.
5. `05_DIRECT_LR` — cross-fitted score-level `Lambda_Z` and `Lambda_{G|Z}`, finite-stage `b_delta`, clinical-tail summaries, decision-flip localization, and the prespecified 3x3 spline/regularization sensitivity grid used in the ILD analysis.
6. `06_BOUNDARY_RESCUE` — leakage-free gene rescue between lower/upper ZeBRA thresholds compared with a ZeBRA-only threshold matched on training-set FPR.

Each run writes a frozen config/parameter/version manifest and the encoded feature list under `RESULTS/<analysis_name>/`.

## ILD reference validation

Run the validated ILD/MUC5B reference comparison with:

```bash
cd zebra_comp_general
./ild_reference.sh --run full
```

The expected cohort invariants are:

- full scored genomic cohort: 12,825;
- full FILD/FILA cases: 254;
- MUC5B-called cohort: 12,766;
- MUC5B-called FILD/FILA cases: 252.

The ILD config contains explicit cohort guards so a sparse-phenotype merge error stops before the expensive analyses run.

## Basic use for a new disease

```bash
cd zebra_comp_general
pip install -r requirements.txt
cp configs/adrd_apoe.template.json configs/adrd_apoe.json
# Edit paths, phenotype semantics, genomic selector, and focal gene signal.
./run_all.sh configs/adrd_apoe.json fast
./run_all.sh configs/adrd_apoe.json full
```

## Genomic column selection

`genomic_features` accepts any of the following selectors.

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

Selected categorical genotype columns are one-hot encoded automatically; numeric genotype/dosage columns remain numeric. Invariant encoded columns are removed. Missing numeric values are imputed inside the training split.

## Gene signal encodings

The focal gene state must be binary for the local/direct-LR/rescue chain. Supported forms are:

- `binary_column`: already encoded 0/1 carrier/status column;
- `dosage_column`: positive if dosage is at or above a threshold;
- `genotype_column`: raw categorical genotype with configured positive values;
- `onehot_any`: current MUC5B-style one-hot state encoding.

The ILD config demonstrates `onehot_any`. ADRD and HFrEF templates assume a precomputed `APOE4_carrier` or `TTNtv_carrier`; change those definitions to the representation actually present in the new cohort.

**Important:** the broad disease-gene header matcher only identifies genomic columns that lie in or map to candidate genes. It does not define the focal biological carrier variable. For ADRD, `APOE4_carrier` should be derived from the actual APOE genotype/haplotype representation (classically rs429358 and rs7412, or an already curated APOE field), not from the presence of an arbitrary APOE-region variant. For HFrEF, `TTNtv_carrier` should be based on an appropriately annotated truncating/pathogenic TTN variant definition; any TTN-region variant is not equivalent to TTNtv carrier status.

## Disease gene-column discovery

The existing Colorado lineage used IPF, ADRD, and HFrEF gene panels. They now live in `gene_panels.json`.

### ADRD

```bash
scripts/get_adrd_gene_columns.sh genomicdataheader.csv ADRD_genomic_header_matches.csv
python scripts/find_apoe4_columns.py genomicdataheader.csv
```

The first command returns the broad ADRD candidate-gene feature set. The second locates columns containing the two APOE-defining SNP IDs (`rs429358`, `rs7412`) so the actual allele/genotype encoding can be inspected before constructing `APOE4_carrier`.

### HFrEF

```bash
scripts/get_hfref_gene_columns.sh genomicdataheader.csv HFREF_genomic_header_matches.csv
```

Or call the generic panel matcher:

```bash
python scripts/find_phenotype_genomic_headers.py genomicdataheader.csv --panel ADRD
python scripts/find_phenotype_genomic_headers.py genomicdataheader.csv --panel HFREF
python scripts/find_phenotype_genomic_headers.py genomicdataheader.csv --panel IPF
```

The utility resolves GRCh37 gene intervals and overlapping rsIDs through the Ensembl GRCh37 REST API and caches the result under `.ensembl_cache/`. Output has exactly `gene_name,header_name`, so it can be fed directly into a config using the `file` selector.

## ADRD/APOE next-analysis workflow

See `ADRD_APOE_RUNBOOK.md`. The intended sequence is:

1. discover the broad ADRD genomic feature columns;
2. locate and inspect `rs429358` and `rs7412` (or use an existing curated APOE genotype field);
3. construct a defensible binary `APOE4_carrier` variable;
4. explicitly specify ADRD phenotype-table semantics (`absent_row=negative` only if absence truly means control; otherwise `drop`);
5. run the unchanged validated pipeline first in `fast` mode and then in `full` mode.

## HFrEF workflow

```bash
scripts/get_hfref_gene_columns.sh /path/to/genomicdataheader.csv HFREF_genomic_header_matches.csv
cp configs/hfref_ttn.template.json configs/hfref_ttn.json
# Fill input paths, HFrEF label, cohort semantics, and actual TTNtv carrier representation.
./run_all.sh configs/hfref_ttn.json fast
./run_all.sh configs/hfref_ttn.json full
```

## Important interpretation constraint

The direct LR analysis is exact for the observed score-level factorization `Lambda_{Z,G}=Lambda_Z Lambda_{G|Z}`. It does not assume that ZeBRA is sufficient for full clinical history, and a finite observed `b_delta` must not be described as a universal genome-wide or sequential bound. The generalized package preserves the same distinction as the ILD manuscript analysis.
