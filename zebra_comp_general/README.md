# General ZeBRA + genomics comparison pipeline

This directory is the disease-agnostic counterpart of `../zebra_comp/`. The original `zebra_comp/` remains the frozen Colorado ILD/IPF manuscript provenance snapshot. Use `zebra_comp_general/` for new disease/gene analyses.

## Validation status

The generalized pipeline is validated against the frozen Colorado ILD/MUC5B analysis. A full reference run at code baseline `8c5db26623c3bc84ac2b0d5bae9bd109752e5dc8` produced:

```text
PASS=38, WARN=0, FAIL=0, INFO=2
```

All manuscript-critical local-information, direct-LR, finite-stage `b_delta`, estimator-robustness, incremental-value, and rescue comparisons passed. The two INFO rows are the intentionally implementation-dependent global nonlinear genomics-only and combined-model AUCs. See `ILD_REFERENCE_VALIDATION.md`.

## Analysis convention: one phenotype per run

Each disease config should contain one primary binary phenotype.

- ILD: `FILD_FILA` + ZeBRA-ILD + ILD genomic panel + MUC5B T-carrier
- ADRD: `ADRD` + ZeBRA-ADRD + ADRD genomic panel + APOE4 carrier
- HFrEF: `HFrEF` + ZeBRA-HFrEF + HF/cardiomyopathy genomic panel + TTNtv carrier

The current `configs/ild_muc5b.json` uses only `FILD or FILA ADJUDICATED`. The historical nonfibrotic ILD/ILA endpoints are not run by the generalized ILD config.

## Replacing ZeBRA risk values

If only the ZeBRA score changes for the same disease/cohort, do not modify analysis code. Change only the `inputs.zebra` block in the config:

```json
"zebra": {
  "path": "/path/to/new_predictions.parquet",
  "id_column": "patient_id",
  "score_column": "predicted_risk",
  "duplicate_policy": "identical"
}
```

If the score is not bounded in `[0,1]`, also change/remove `score_range`.

Everything downstream is recalculated automatically: AUC, incremental genomic value, local FPR bands, `Lambda_Z`, `Lambda_{G|Z}`, `b_hat_delta`, tail summaries, decision-flip localization, and rescue.

For the validated ILD reference cohort the guards are:

- full scored genomic cohort: 12,825
- full FILD/FILA cases: 254
- MUC5B-called cohort: 12,766
- MUC5B-called FILD/FILA cases: 252

## Run sequence

```text
01_GLOBAL_COMPARISON
02_INCREMENTAL_LOGISTIC
03_GENE_SCORE_ASSOCIATION
04_LOCAL_INFORMATION
05_DIRECT_LR
06_BOUNDARY_RESCUE
```

Run with:

```bash
./run_all.sh CONFIG.json fast
./run_all.sh CONFIG.json full
```

For ILD regression validation:

```bash
./ild_reference.sh --run full
```

## Genomic feature selection

The config can select genomic columns by all columns, explicit names, regex/prefix, a contiguous column slice, or a CSV containing `header_name`.

Example disease-panel file:

```json
"genomic_features": {
  "type": "file",
  "path": "../ADRD_genomic_header_matches.csv",
  "column_field": "header_name"
}
```

The disease gene panels are defined in `gene_panels.json`.

Generate matched columns from a genomic header file with:

```bash
scripts/get_adrd_gene_columns.sh genomicdataheader.csv ADRD_genomic_header_matches.csv
scripts/get_hfref_gene_columns.sh genomicdataheader.csv HFREF_genomic_header_matches.csv
```

For ADRD, locate the APOE-defining SNP columns with:

```bash
python scripts/find_apoe4_columns.py genomicdataheader.csv --output APOE4_defining_columns.csv
```

The broad gene-column match is only the feature panel. It does not define the focal carrier variable.

- `APOE4_carrier` must come from verified APOE genotype encoding, normally rs429358 + rs7412 or an existing curated APOE field.
- `TTNtv_carrier` must come from an annotated TTN truncating/pathogenic-variant definition; any TTN-region variant is not sufficient.

# TODO: exact steps for the next runs

Do not change the validated analysis modules. Only produce the required input files/columns and edit the configs.

## ADRD / APOE4

1. Generate the ADRD genomic-column list from the genomic header:

```bash
cd zebra_comp_general
scripts/get_adrd_gene_columns.sh \
    /path/to/genomicdataheader.csv \
    ADRD_genomic_header_matches.csv
```

Output required:

```text
ADRD_genomic_header_matches.csv
```

2. Find the APOE-defining columns:

```bash
python scripts/find_apoe4_columns.py \
    /path/to/genomicdataheader.csv \
    --output APOE4_defining_columns.csv

cat APOE4_defining_columns.csv
```

Output required:

```text
APOE4_defining_columns.csv
```

It should identify the representations of `rs429358` and `rs7412`.

3. From the actual genotype encoding, create one patient-level binary column:

```text
APOE4_carrier
```

Coding:

```text
1 = >=1 epsilon-4 allele
0 = confidently called non-carrier
NA = uncalled/ambiguous
```

Add this column to the genomic matrix used by the run.

4. Obtain/provide the ADRD ZeBRA score file. It must contain at least:

```text
patient_id
predicted_risk
```

(or use the actual ID/score column names in the config).

5. Copy the config:

```bash
cp configs/adrd_apoe.template.json configs/adrd_apoe.json
```

6. In `configs/adrd_apoe.json`, set only these cohort-specific fields:

```text
inputs.genomics.path
inputs.genomics.id_column
inputs.phenotypes.path
inputs.phenotypes.id_column
inputs.zebra.path
inputs.zebra.id_column
inputs.zebra.score_column
phenotypes[0].column
phenotypes[0].missing
phenotypes[0].absent_row
genomic_features.path
gene_signal.column
```

Expected settings for the generated files are:

```json
"genomic_features": {
  "type": "file",
  "path": "../ADRD_genomic_header_matches.csv",
  "column_field": "header_name"
},
"gene_signal": {
  "enabled": true,
  "name": "APOE4_carrier",
  "type": "binary_column",
  "column": "APOE4_carrier"
}
```

7. Run:

```bash
./run_all.sh configs/adrd_apoe.json fast
cat RESULTS/adrd_apoe/RUN_MANIFEST.json
./run_all.sh configs/adrd_apoe.json full
```

That is the complete ADRD execution path.

## HFrEF / TTNtv

1. Generate the HFrEF genomic-column list:

```bash
cd zebra_comp_general
scripts/get_hfref_gene_columns.sh \
    /path/to/genomicdataheader.csv \
    HFREF_genomic_header_matches.csv
```

Output required:

```text
HFREF_genomic_header_matches.csv
```

2. Create one patient-level binary focal-gene column from the annotated genomic data:

```text
TTNtv_carrier
```

Coding:

```text
1 = >=1 qualifying TTN truncating/pathogenic variant
0 = confidently callable non-carrier
NA = uncalled/ambiguous
```

Add this column to the genomic matrix. Do not define `TTNtv_carrier` from arbitrary TTN-region variation.

3. Obtain/provide the HFrEF ZeBRA score file. It must contain at least:

```text
patient_id
predicted_risk
```

(or use the actual ID/score column names in the config).

4. Copy the config:

```bash
cp configs/hfref_ttn.template.json configs/hfref_ttn.json
```

5. In `configs/hfref_ttn.json`, set only these cohort-specific fields:

```text
inputs.genomics.path
inputs.genomics.id_column
inputs.phenotypes.path
inputs.phenotypes.id_column
inputs.zebra.path
inputs.zebra.id_column
inputs.zebra.score_column
phenotypes[0].column
phenotypes[0].missing
phenotypes[0].absent_row
genomic_features.path
gene_signal.column
```

Expected settings for the generated files are:

```json
"genomic_features": {
  "type": "file",
  "path": "../HFREF_genomic_header_matches.csv",
  "column_field": "header_name"
},
"gene_signal": {
  "enabled": true,
  "name": "TTNtv_carrier",
  "type": "binary_column",
  "column": "TTNtv_carrier"
}
```

6. Run:

```bash
./run_all.sh configs/hfref_ttn.json fast
cat RESULTS/hfref_ttn/RUN_MANIFEST.json
./run_all.sh configs/hfref_ttn.json full
```

That is the complete HFrEF execution path.

## Important interpretation constraint

The direct LR analysis is exact for the observed score-level factorization `Lambda_{Z,G}=Lambda_Z Lambda_{G|Z}`. It does not assume ZeBRA is sufficient for the full clinical history, and a finite observed `b_delta` is not a universal genome-wide or sequential bound.
