# Curated Colorado manuscript analysis

This directory is the compact reproducibility workspace for the ZeBRA + genomics analyses used in `tex/zebra_genomics_boundary_theorem_empirical.tex`.

The core notebooks are copied from `ZEBRA_GENDRIVERS_CLASSIFIER/v1` without changing their analysis logic. Supporting v1 inputs are copied here so the downstream notebooks can be run from this directory. The raw Colorado genotype matrix is intentionally kept separate under `colorado_data/` and is not included yet.

## Core notebook lineage

Run from the `curated/` directory.

1. `01_GATHER_GENDRIVER_DATA.ipynb`
   - Builds the targeted genomic matrix used by the v1 analysis.
   - Matches the expanded `MORE_LOCI.csv` list to available genotype columns, unions those columns with the pre-existing biological-driver panel, and writes `ILD_TOP_DRIVERS_DATA.csv`.
   - The archived notebook contains the original cluster-specific `GENDIR` path. For a fresh reconstruction, set `GENDIR` in that notebook to `./colorado_data` after placing the Colorado genomic files there.

2. `02_GEN_DRIVER_CLASSIFIERS.ipynb`
   - Fits the genomic-only classifiers for the three adjudicated ILD/FILA targets.
   - Produces the genomic-only portion of the fixed-split AUC comparison.

3. `03_COMBINED_CLASSIFIERS.ipynb`
   - Fits/evaluates ZeBRA and combined ZeBRA + genomic models using the archived notebook-03 target construction.
   - Uses `PREDICTIONS_104W_PRED_WINDOW.parquet`.
   - Produces the fixed-split ZeBRA and combined results reported in the manuscript.

4. `032_COMBINED_CLASSIFIERS_REPEATED_SPLITS.ipynb`
   - Repeats the notebook-03 target logic over randomized outer splits and evaluates ZeBRA, genomic-only, and combined models.
   - Uses `PREDICTIONS_52W_PRED_WINDOW.parquet`. That prediction file was stored under v2 in the source repository, so it is copied here solely to satisfy this archived v1 notebook dependency; the notebook itself remains the v1 version.

The larger legacy `32_COMBINED_CLASSIFIERS.ipynb` and `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb` are intentionally not copied into this curated core because the current empirical TeX and the v1 provenance note identify notebooks 01/02/03/032 as the relevant manuscript lineage.

## Included supporting inputs

- `ILD_TOP_DRIVERS_DATA.csv` — archived v1 derived genomic matrix. This allows notebooks 02/03/032 to run even before the raw Colorado genotype matrix is staged locally.
- `MORE_LOCI.csv` — expanded candidate-locus list used by notebook 01.
- `REPHENOTYPES FOR IC.csv` — adjudicated target labels.
- `PREDICTIONS_104W_PRED_WINDOW.parquet` — input used by notebook 03.
- `PREDICTIONS_52W_PRED_WINDOW.parquet` — input required by notebook 032.
- `SNP_PANEL_USED.csv` — reconstructed 171-locus panel documented in the manuscript supplement.
- `PHENOTYPES_AND_GENOMIC_PANEL.md` — archived provenance note for target construction and genomic panel.

## Curated gene-to-header mapping

`./scripts/find_phenotype_genomic_headers.py` maps curated disease-associated gene panels to the actual variant names in the one-column Colorado `genomicdataheader.csv`. It currently provides panels for:

- `IPF` — IPF/fibrotic ILD
- `ADRD` — Alzheimer's disease and related dementias
- `HFREF` — HFrEF/dilated-cardiomyopathy genetics

Example:

```bash
python scripts/find_phenotype_genomic_headers.py \
    colorado_data/genomicdataheader.csv \
    --panel IPF \
    --output-dir colorado_data
```

The output is a two-column CSV:

```text
gene_name,header_name
```

Run all three panels with `--panel ALL`.

## Colorado raw data

Place the full Colorado genomic inputs under `colorado_data/`. See `colorado_data/README.md` for expected files.

The raw patient-level genomic data should only be placed in Git if its governance and sharing permissions allow that. The curated workflow itself does not require committing the raw matrix: notebooks can be run with locally staged data.

## Python environment

A minimal dependency list is provided in `requirements.txt`.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Manuscript linkage

The current empirical manuscript is:

`../tex/zebra_genomics_boundary_theorem_empirical.tex`

The archived fixed-split table in that manuscript reports approximately:

| Target | ZeBRA | Genomic | Combined |
| --- | ---: | ---: | ---: |
| FILD/FILA | 0.804 | 0.643 | 0.800 |
| Nonfibrotic ILD | 0.838 | 0.629 | 0.855 |
| Nonfibrotic ILA | 0.809 | 0.573 | 0.806 |

Use these values as regression checks when rerunning the archived fixed-split workflow.
