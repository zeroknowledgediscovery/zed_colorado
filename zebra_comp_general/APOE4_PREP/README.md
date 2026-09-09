# APOE4 preparation for ADRD

This folder isolates the APOE4 preparation step from the main `zebra_comp_general` analysis.

The goal is to create exactly one patient-level column:

```text
APOE4_carrier
```

with:

```text
1   confidently carries >=1 APOE epsilon-4 allele
0   confidently called non-carrier
NA  uncalled / ambiguous / noncanonical
```

Do not use arbitrary APOE-region variants for this. APOE epsilon status is defined from rs429358 and rs7412 (or an already curated APOE genotype field).

## Step 1: find the two defining columns

Run against the full genomic header list:

```bash
python 01_find_apoe_columns.py /path/to/genomicdataheader.csv
```

This writes:

```text
APOE_candidates.csv
```

The script uses exact rsID matching and GRCh37 coordinate fallback. It will not confuse `rs7412` with `rs74127625`, `rs74123749`, etc.

Expected GRCh37 targets:

```text
rs429358  chr19:45411941  T/C
rs7412    chr19:45412079  C/T
```

## Step 2: inspect how those columns are encoded

```bash
python 02_inspect_apoe_columns.py \
    /path/to/genomicdata.csv \
    APOE_candidates.csv \
    --id-column FID
```

This writes:

```text
APOE_column_value_counts.csv
APOE_first_rows.csv
```

Inspect these outputs. Determine which column represents rs429358, which represents rs7412, and which allele the 0/1/2 dosage counts.

For example, if a column is named `rs429358_C` and values are 0/1/2, it may be C dosage. Verify this from the export documentation/header convention before proceeding.

## Step 3: create APOE4_carrier

After the dosage allele is known, run explicitly:

```bash
python 03_make_apoe4_carrier.py \
    /path/to/genomicdata.csv \
    --id-column FID \
    --rs429358-column 'EXACT_COLUMN_NAME' \
    --rs429358-dosage-allele C \
    --rs7412-column 'EXACT_COLUMN_NAME' \
    --rs7412-dosage-allele T \
    --output APOE4_carrier.csv
```

`--rs429358-dosage-allele` must be `C` or `T`.
`--rs7412-dosage-allele` must be `C` or `T`.

The script converts both loci to C-allele dosage internally and uses a conservative two-SNP APOE mapping.

By default, the double-heterozygous state `(rs429358 C dosage=1, rs7412 C dosage=1)` is left missing because phase is not established from unphased dosage alone. If the dataset documentation/curated genotype convention establishes that this state should be treated as epsilon-2/epsilon-4, rerun with:

```bash
--double-het-e2e4
```

The output contains:

```text
FID
rs429358_C_dosage
rs7412_C_dosage
APOE_genotype
APOE4_carrier
```

## Step 4: validate counts

```bash
python 04_check_apoe4.py APOE4_carrier.csv --id-column FID
```

Optionally include the ADRD phenotype file:

```bash
python 04_check_apoe4.py APOE4_carrier.csv \
    --id-column FID \
    --phenotype-file /path/to/ADRD_labels.csv \
    --phenotype-id-column FID \
    --phenotype-column ADRD \
    --positive-values 1,Y,y
```

Only after these counts look plausible should `APOE4_carrier` be merged into the genomic matrix/config used by `zebra_comp_general`.
