# Curated analyses for the ZeBRA–genomics manuscript

This folder contains the analysis lineage needed for the current manuscript framing. The canonical combined-classifier analysis is `32_COMBINED_CLASSIFIERS.ipynb`; the earlier fixed-split `03_COMBINED_CLASSIFIERS.ipynb` is intentionally not retained.

## Manuscript result structure

The retained analyses support three linked conclusions:

1. Broadly adding genomic features to ZeBRA does not materially improve global discrimination.
2. Genomic information can add local information in selected regions of the ZeBRA score distribution.
3. MUC5B can improve classification when used selectively at a ZeBRA decision boundary, rather than being added indiscriminately to the global predictor.

The separate `test_fixed_muc5b_switch_regions.py` experiment remains deliberately excluded. The local-information and ROC-hull analyses are retained as supporting analyses; the matched-FPR high/low rescue analysis is the primary decision-boundary result.

## Core notebooks

### `01_GATHER_GENDRIVER_DATA.ipynb`
Upstream provenance for construction of the genomic feature matrix. Downstream analyses can start from the included processed `ILD_TOP_DRIVERS_DATA.csv`.

### `02_GEN_DRIVER_CLASSIFIERS.ipynb`
Genomics-only comparator models.

### `32_COMBINED_CLASSIFIERS.ipynb`
Canonical/latest combined-classifier notebook for this v0 lineage. It uses the 104-week ZeBRA predictions, repeated randomized outer splits, the notebook-03 phenotype/target logic, and repeated ZeBRA/genomic/combined comparisons with zedstat/calibration/SHAP analyses.

For FILD/FILA, the retained repeated-split summary is approximately:

- ZeBRA mean AUC: 0.8141
- genomics mean AUC: 0.5623
- combined mean AUC: 0.7888

Thus broad genomic augmentation does not improve the principal global ZeBRA result.

### `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb`
Conditional incremental-value analysis on the notebook-32 cohort. Adding the genomic panel after ZeBRA produces only a very small mean AUC change (about +0.0007, median 0) and does not provide robust global improvement.

## Local genomic information and decision-boundary analyses

### `test_local_zebra_genomics_predictive_curves.py`
Maps local ZeBRA information and local MUC5B information across overlapping windows of the empirical ZeBRA score distribution. This provides the descriptive/local-information evidence motivating the view that genomics can be useful in restricted regions even when it is weak globally.

Outputs are written automatically to:

`RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/`

### `test_muc5b_zebra_rescue_matched_fpr_lr.py`
Primary decision-boundary analysis. For each repeated 40/60 train/test split:

1. A stringent high ZeBRA threshold is derived from training negatives at a target FPR of 0.5% or 1%.
2. A less stringent low ZeBRA threshold is derived at a target FPR of 2%, 5%, or 10%.
3. Patients above the high threshold are ZeBRA-positive without genomics.
4. Patients below the low threshold remain negative.
5. Only patients in the intermediate band are eligible for MUC5B-based rescue/reclassification.
6. Performance is compared with a ZeBRA-only classifier matched to the empirical FPR of the rescue rule.

Representative 100-split results include:

| high-threshold FPR | low-threshold FPR | rescue sensitivity | matched ZeBRA sensitivity | sensitivity gain | rescue FPR | matched ZeBRA FPR |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5% | 5% | 0.3486 | 0.3126 | +0.0360 | 0.01405 | 0.01430 |
| 0.5% | 10% | 0.4012 | 0.3584 | +0.0428 | 0.02326 | 0.02374 |
| 1% | 5% | 0.3529 | 0.3320 | +0.0209 | 0.01837 | 0.01870 |

Outputs are written automatically to:

`RESULTS/MUC5B_ZEBRA_RESCUE_RULE/`

### `test_zebra_hybrid_roc_convex_hull_zedstat.py`
Supporting ROC/zedstat analysis. It compares ZeBRA, the cross-fitted MUC5B hybrid, and the upper ROC convex hull formed from the union of both empirical ROC curves. The hull should be interpreted as the attainable operating envelope obtained by selecting between the two classifiers across operating regions, not as a third scalar classifier.

This script uses the fixed hybrid windows encoded in the script for its ROC-envelope analysis; it is retained as supporting geometric/operating-point evidence rather than as the primary boundary rule.

Outputs are written automatically to:

`RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL/`

## Retained result directories

- `RESULTS/GENDRIVERS/` — genomic-only comparator.
- `RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/` — compact repeated-split/global summaries from notebook 32.
- `RESULTS/033_INCREMENTAL_LOGISTIC_32_COHORT/` — compact conditional incremental-value summaries from notebook 33.
- `RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/` — local-information curves and candidate local regions.
- `RESULTS/MUC5B_ZEBRA_RESCUE_RULE/` — high/low-threshold matched-FPR rescue analysis.
- `RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL/` — ZeBRA/hybrid ROC and zedstat convex-hull analysis.

## Environment

From `curated_1/`:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-analysis.txt
```

The requirements file installs the scientific Python stack used here and the ZKDL `zedstat` package from its GitHub repository.

## Recreate the three local/boundary result folders

All three scripts use relative paths and can be run directly from this directory. The required processed inputs are already present:

- `ILD_TOP_DRIVERS_DATA.csv`
- `REPHENOTYPES FOR IC.csv`
- `PREDICTIONS_104W_PRED_WINDOW.parquet`

For a clean regeneration, either run:

```bash
bash run_boundary_analyses.sh
```

or explicitly:

```bash
rm -rf RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES
rm -rf RESULTS/MUC5B_ZEBRA_RESCUE_RULE
rm -rf RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL

python test_local_zebra_genomics_predictive_curves.py
python test_muc5b_zebra_rescue_matched_fpr_lr.py
python test_zebra_hybrid_roc_convex_hull_zedstat.py
```

Each script creates its own result directory automatically.

## Recreate the notebook-derived global results

Starting from the included processed input files, execute the notebooks in this order:

```bash
jupyter nbconvert --to notebook --execute --inplace 02_GEN_DRIVER_CLASSIFIERS.ipynb
jupyter nbconvert --to notebook --execute --inplace 32_COMBINED_CLASSIFIERS.ipynb
jupyter nbconvert --to notebook --execute --inplace 33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb
```

The important notebook-derived result directories are:

- `RESULTS/GENDRIVERS/`
- `RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/`
- `RESULTS/033_INCREMENTAL_LOGISTIC_32_COHORT/`

`01_GATHER_GENDRIVER_DATA.ipynb` is upstream provenance only. Re-running it from scratch requires the restricted/raw Colorado genotype source that is not part of this curated package.

## Recommended full reproduction order

1. `02_GEN_DRIVER_CLASSIFIERS.ipynb`
2. `32_COMBINED_CLASSIFIERS.ipynb`
3. `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb`
4. `test_local_zebra_genomics_predictive_curves.py`
5. `test_muc5b_zebra_rescue_matched_fpr_lr.py`
6. `test_zebra_hybrid_roc_convex_hull_zedstat.py`

## Deliberately excluded

- `03_COMBINED_CLASSIFIERS.ipynb`: superseded by notebook 32 for this curated analysis.
- `test_fixed_muc5b_switch_regions.py` and `RESULTS/FIXED_MUC5B_SWITCH_REGIONS/`: not part of the intended manuscript analysis.
- the original `test_muc5b_zebra_rescue_rule.py`: superseded by the matched-FPR/LR implementation.
- MUC5B error-strata scripts: audit/exploratory analyses not required for the manuscript claim.
- older 032 variants, ExtraTrees/isotonic branches, checkpoints, backups, fitted model pickles, and unrelated exploratory analyses.
