# Curated analyses for the ZeBRA–genomics manuscript

This folder contains only the analysis lineage needed for the current manuscript framing. The canonical combined-classifier analysis is notebook `32_COMBINED_CLASSIFIERS.ipynb`; notebook `03_COMBINED_CLASSIFIERS.ipynb` is an earlier fixed-split analysis and is intentionally not retained here.

## Manuscript result structure

The retained analyses support two linked conclusions:

1. Broadly adding genomic features to ZeBRA does not materially improve global discrimination.
2. Genomic information can nevertheless improve classification locally at a ZeBRA decision boundary when it is used selectively between a stringent high threshold and a lower rescue threshold.

The second result is the boundary result retained here. The separate fixed-percentile MUC5B switch experiment is deliberately excluded.

## Core notebooks

### `01_GATHER_GENDRIVER_DATA.ipynb`
Upstream provenance for construction of the genomic feature matrix. Downstream analyses can start from the included processed `ILD_TOP_DRIVERS_DATA.csv`.

### `02_GEN_DRIVER_CLASSIFIERS.ipynb`
Genomics-only comparator models.

### `32_COMBINED_CLASSIFIERS.ipynb`
**Canonical/latest combined-classifier notebook for this v0 lineage.** It supersedes notebook 03 for the curated analysis. It uses the 104-week ZeBRA predictions, repeated randomized outer splits, and the notebook-03 phenotype/target logic, with repeated ZeBRA/genomic/combined comparisons plus the zedstat/calibration/SHAP analyses.

`32_COMBINED_CLASSIFIERS.ipynb` is the same Git blob as the later-named `032_COMBINED_CLASSIFIERS_REPEATED_SPLITS_ZEDSTAT_FIXED.ipynb`; only `32` is retained to avoid duplication.

For FILD/FILA across the retained repeated splits, the summary is approximately:
- ZeBRA mean AUC: 0.8141
- genomics mean AUC: 0.5623
- combined mean AUC: 0.7888

Thus indiscriminate genomic augmentation does not improve the principal global ZeBRA result.

### `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb`
Conditional incremental-value analysis on the notebook-32 cohort. Adding the genomic panel after ZeBRA produces only a very small mean AUC change (about +0.0007, median 0) and does not provide robust global improvement.

## Decision-boundary analysis

### `test_muc5b_zebra_rescue_matched_fpr_lr.py`
**Primary retained boundary analysis.** This is the analysis corresponding to the manuscript concept that genomics can improve prediction selectively between a high and low ZeBRA threshold.

For each repeated 40/60 train/test split:

1. A stringent **high ZeBRA threshold** is derived from training negatives at a target FPR of 0.5% or 1%.
2. A less stringent **low ZeBRA threshold** is derived at a target FPR of 2%, 5%, or 10%.
3. Patients above the high threshold are ZeBRA-positive without using genomics.
4. Patients below the low threshold remain negative.
5. Only patients in the intermediate band are eligible for MUC5B-based rescue/reclassification.
6. Performance is compared with a ZeBRA-only classifier matched to the empirical FPR of the rescue rule.

This directly tests whether genomic information improves the classifier at the decision boundary rather than whether it raises global AUC.

Representative 100-split results from `RESULTS/MUC5B_ZEBRA_RESCUE_RULE/MATCHED_FPR_RESCUE_VS_ZEBRA.csv` include:

| high-threshold FPR | low-threshold FPR | rescue sensitivity | matched ZeBRA sensitivity | sensitivity gain | rescue FPR | matched ZeBRA FPR |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5% | 5% | 0.3486 | 0.3126 | +0.0360 | 0.01405 | 0.01430 |
| 0.5% | 10% | 0.4012 | 0.3584 | +0.0428 | 0.02326 | 0.02374 |
| 1% | 5% | 0.3529 | 0.3320 | +0.0209 | 0.01837 | 0.01870 |

The interpretation is not that genomics globally improves ZeBRA. Rather, selective MUC5B use inside a predefined uncertainty/rescue band can recover additional cases at nearly matched false-positive burden.

## Retained result directories

- `RESULTS/GENDRIVERS/` — genomic-only comparator.
- `RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/` — compact repeated-split/global summaries from notebook 32.
- `RESULTS/033_INCREMENTAL_LOGISTIC_32_COHORT/` — compact conditional incremental-value summaries from notebook 33.
- `RESULTS/MUC5B_ZEBRA_RESCUE_RULE/` — the high/low-threshold decision-boundary experiment, including matched-FPR sensitivity, LR+, LR−, diagnostic-odds-ratio, and rescue summaries.

## Deliberately excluded

- `03_COMBINED_CLASSIFIERS.ipynb` and its fixed-split-only result provenance: superseded in this curated package by notebook 32.
- `test_fixed_muc5b_switch_regions.py` and `RESULTS/FIXED_MUC5B_SWITCH_REGIONS/`: not part of the intended manuscript analysis.
- the original `test_muc5b_zebra_rescue_rule.py`: superseded by the matched-FPR/LR implementation retained above.
- `test_local_zebra_genomics_predictive_curves.py` and its exploratory window-search results: not required for the final high/low boundary claim.
- `testinverse_muc5b_error_strata.py` and its error-strata outputs: audit/exploratory analysis, not required for the manuscript claim.
- older 032 variants, ExtraTrees/isotonic branches, checkpoints, backup files, fitted model pickles, and unrelated exploratory analyses.

## Reproduction order

Starting from the included processed inputs:

1. `02_GEN_DRIVER_CLASSIFIERS.ipynb`
2. `32_COMBINED_CLASSIFIERS.ipynb`
3. `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb`
4. `test_muc5b_zebra_rescue_matched_fpr_lr.py`

`01_GATHER_GENDRIVER_DATA.ipynb` is retained only for upstream data provenance.
