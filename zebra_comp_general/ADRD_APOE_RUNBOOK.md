# ADRD/APOE runbook

This is the next disease/gene application of the ILD-validated `zebra_comp_general` pipeline. The statistical analysis code should not be modified for ADRD; only the disease/gene configuration and the cohort-specific input preparation should change.

## Required inputs

1. **Genomic matrix** with one row per patient and a stable patient ID.
2. **ADRD phenotype table** with an ADRD case/control label or a sparse case/adjudication table whose absence semantics are explicitly known.
3. **ADRD ZeBRA predictions** with one risk score per patient for the intended prediction horizon.
4. **Genomic header list** corresponding to the genomic matrix.
5. **APOE4 carrier variable**, either already curated or derivable from verified rs429358/rs7412 genotype encoding.

Do not begin the full analysis until the phenotype-table semantics and APOE encoding are verified.

## 1. Discover the broad ADRD genomic feature panel

From `zebra_comp_general/`:

```bash
scripts/get_adrd_gene_columns.sh /path/to/genomicdataheader.csv ADRD_genomic_header_matches.csv
```

This produces a `gene_name,header_name` file spanning the configured ADRD candidate-gene panel. It is used for the broad genomic feature set in the global and incremental analyses.

Review its size and contents before proceeding:

```bash
head ADRD_genomic_header_matches.csv
wc -l ADRD_genomic_header_matches.csv
```

## 2. Locate the APOE-defining SNP columns

```bash
python scripts/find_apoe4_columns.py /path/to/genomicdataheader.csv --output APOE4_defining_columns.csv
cat APOE4_defining_columns.csv
```

The expected SNP IDs are:

- `rs429358`
- `rs7412`

The script only locates matching columns. It deliberately does not infer APOE alleles or APOE4 status.

## 3. Construct and validate APOE4_carrier

Preferred order:

1. use an existing curated APOE genotype/carrier field if one exists;
2. otherwise derive APOE genotype from rs429358 and rs7412 after verifying the export's reference/alternate allele coding and genotype representation;
3. define `APOE4_carrier=1` for subjects carrying at least one epsilon-4 allele and `0` for confidently called non-carriers;
4. leave uncalled/ambiguous APOE genotypes missing rather than forcing them to zero.

Before running the analysis, report:

- total patients;
- APOE-called patients;
- APOE4 carriers;
- APOE4 carrier prevalence;
- ADRD cases among APOE-called patients;
- carrier prevalence in cases and controls.

These become ADRD cohort invariants analogous to the ILD MUC5B checks.

## 4. Decide ADRD phenotype-table semantics

The template defaults to:

```json
"missing": "drop",
"absent_row": "drop"
```

This is intentionally conservative.

Change `absent_row` to `negative` only if the ADRD phenotype file is known to be a sparse case/adjudication table in which genomic-cohort patients absent from the file are valid controls. Do not infer that from sparsity alone.

If the phenotype table explicitly contains every cohort member and uses a blank cell for negative, set `missing` according to that documented encoding.

## 5. Create the runnable ADRD config

```bash
cp configs/adrd_apoe.template.json configs/adrd_apoe.json
```

Edit:

- `inputs.genomics.path`
- `inputs.genomics.id_column`
- `inputs.phenotypes.path`
- `inputs.phenotypes.id_column`
- `inputs.zebra.path`
- `inputs.zebra.id_column`
- `inputs.zebra.score_column`
- phenotype `column`
- `missing` and `absent_row`
- `gene_signal.column` or encoding if `APOE4_carrier` is not already a binary column

The default broad genomic selector expects:

```json
{
  "type": "file",
  "path": "../ADRD_genomic_header_matches.csv",
  "column_field": "header_name"
}
```

## 6. Smoke test

```bash
./run_all.sh configs/adrd_apoe.json fast
```

Inspect:

```bash
cat RESULTS/adrd_apoe/RUN_MANIFEST.json
cat RESULTS/adrd_apoe/01_GLOBAL_COMPARISON/auc_distribution_summary.csv
cat RESULTS/adrd_apoe/02_INCREMENTAL_LOGISTIC/paired_incremental_summary.csv
cat RESULTS/adrd_apoe/05_DIRECT_LR/stochastic_bounds.csv
```

The manifest should show a plausible cohort size, case count, and `gene_called_n`. If those are surprising, resolve the cohort construction before a full run.

## 7. Full ADRD/APOE run

```bash
./run_all.sh configs/adrd_apoe.json full
```

The primary questions are the same as in ILD:

1. Does the broad ADRD genomic panel improve global discrimination beyond ZeBRA?
2. Does APOE4 retain conditional predictive information given ZeBRA?
3. Where along the ZeBRA score axis is that information concentrated?
4. What is the empirical score-level residual-genomic scale `b_hat_delta`?
5. Does APOE4 change decisions primarily near operating boundaries?
6. Does an APOE4 rescue rule improve sensitivity at matched FPR?

## 8. Interpretation

The ADRD result should be compared to ILD at the level of geometry, not by requiring identical numerical values. In particular, compare:

- global incremental delta AUC;
- local `LRT_G_given_Z` profile;
- APOE4 odds ratio across score/FPR bands;
- `b_hat_0.05`, `b_hat_0.01`, and robustness range;
- fraction of decisions changed near the boundary;
- rescue sensitivity gain at matched FPR.

The direct LR analysis remains a score-level empirical analogue. A finite `b_hat_delta` for APOE4 does not establish a universal genome-wide or sequential bound.
