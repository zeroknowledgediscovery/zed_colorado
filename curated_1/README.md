# Curated v0 analyses for the ZeBRA–genomics manuscript

This folder isolates the v0 files that directly support the current manuscript results and framing. It intentionally excludes checkpoints, superseded model variants, and exploratory analyses that are not needed for the present claims.

## Bottom line

The current evidence supports a two-part conclusion:

1. **Genomics does not materially improve global discrimination when added broadly to ZeBRA.** This is supported by the fixed-split analyses (`02`, `03`), the repeated-split analysis (`32`), and the conditional incremental analysis (`33`).
2. **Genomics can add local information near selected parts of the ZeBRA decision surface.** The strongest current evidence is the 100-split fixed MUC5B switch analysis. MUC5B is used only in pre-specified ZeBRA-percentile regions, while calibrated ZeBRA is used elsewhere.

The high/low-threshold MUC5B rescue analysis is retained because it is the direct analysis lineage behind the idea that genotype can selectively reclassify patients between a stringent high threshold and a lower threshold. It should be treated as supporting/exploratory; the fixed-switch experiment is the cleaner result for the manuscript.

## Core notebooks retained

### `01_GATHER_GENDRIVER_DATA.ipynb`
Constructs the processed genomic feature matrix. The downstream processed file `ILD_TOP_DRIVERS_DATA.csv` is included, so the downstream manuscript analyses do not require rebuilding the raw genotype matrix. Re-running `01` from scratch still requires the restricted/raw Colorado genotype source used originally.

### `02_GEN_DRIVER_CLASSIFIERS.ipynb`
Genomics-only LightGBM models for the adjudicated phenotypes. Supplies the genomic-only comparator and feature-importance analyses.

### `03_COMBINED_CLASSIFIERS.ipynb`
Original fixed 40% train / 60% held-out comparison of ZeBRA, genomics, and combined ZeBRA+genomics. This is the provenance of the fixed-split manuscript comparison.

### `32_COMBINED_CLASSIFIERS.ipynb`
**Required.** Final v0 repeated-split analysis using `PREDICTIONS_104W_PRED_WINDOW.parquet`, with repeated ZeBRA/genomic/combined comparisons and zedstat performance/calibration analyses. This file has the same Git blob as `032_COMBINED_CLASSIFIERS_REPEATED_SPLITS_ZEDSTAT_FIXED.ipynb`; only `32` is retained here to avoid duplication.

For FILD/FILA over 30 repeated splits, the retained summary gives approximately:
- ZeBRA mean AUC: 0.8141
- genomics mean AUC: 0.5623
- combined mean AUC: 0.7888

Thus broad genomic augmentation does not improve the global ZeBRA ranking for the principal FILD/FILA outcome.

### `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb`
Conditional incremental-value analysis asking the more appropriate question: after ZeBRA is known, does genomics add predictive information? Across 30 splits, the all-genomics-plus-ZeBRA model has only a very small mean AUC delta (~+0.0007; median 0), while average precision is lower on average and Brier score/log loss are worse. This supports the manuscript framing of little robust *global* incremental value from broad genomic augmentation.

## Local genomic information / decision-boundary analyses

### `test_local_zebra_genomics_predictive_curves.py`
Exploratory localization analysis. It maps where local MUC5B discrimination/information can exceed local ZeBRA information across the ZeBRA score distribution. The outputs in `RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/` provide the provenance for selecting fixed regions for the subsequent validation experiment. These results are hypothesis-generating, not the final inferential result.

### `test_muc5b_zebra_rescue_rule.py`
Original high/low-threshold rescue formulation. A stringent ZeBRA threshold defines baseline positives; a lower ZeBRA threshold defines an intermediate rescue band; MUC5B T-carriers in that band can be promoted. It uses 100 repeated 40/60 splits and considers main FPR targets of 0.5% and 1% with lower-band targets of 2%, 5%, and 10%.

### `test_muc5b_zebra_rescue_matched_fpr_lr.py`
Refined version of the high/low-threshold rescue experiment with a fair matched-empirical-FPR ZeBRA comparator and LR+/LR-/diagnostic-odds-ratio calculations. This is the preferred script for interpreting the high/low rescue concept. Its outputs are retained under `RESULTS/MUC5B_ZEBRA_RESCUE_RULE/`.

### `test_fixed_muc5b_switch_regions.py`
**Primary boundary result for the current manuscript.** Uses the exact notebook-32 FILD/FILA cohort, 104-week ZeBRA predictions, 40% training / 60% held-out testing, and 100 repeated splits. The switch windows are fixed at:

- P00–P50
- P65–P75
- P85–P95

Percentiles are defined from the training distribution only. Outside these regions, calibrated ZeBRA is used. Inside them, exact ZeBRA ranking is deliberately replaced by training-estimated MUC5B genotype-specific risk.

The current manuscript result at approximately 5% FPR is:

| metric | ZeBRA | fixed MUC5B hybrid |
|---|---:|---:|
| global AUC | 0.8160 | 0.8135 |
| sensitivity | 0.4576 | 0.4969 |
| realized FPR | 0.0586 | 0.0591 |
| LR+ | 7.89 | 8.46 |
| LR- | 0.576 | 0.535 |

Thus sensitivity increases by about 0.0393 at essentially unchanged FPR, even though global AUC is not improved. This is the clearest empirical support for the manuscript's central framing: genomic information can be useful **locally at the decision boundary without improving global discrimination**.

Mean local MUC5B odds ratios in the three fixed switch regions are approximately 2.88, 4.40, and 4.07, respectively.

## Error-strata audit

### `testinverse_muc5b_error_strata.py`
Retained as an audit/supporting analysis only. It asks whether MUC5B is enriched among ZeBRA false negatives/false positives at selected thresholds. The false-negative enrichment is not strong enough to serve as the headline justification for the boundary claim (for example, at the ~5% operating point the FN-vs-TP MUC5B OR is only about 1.3 and is not statistically compelling). The manuscript should therefore base the boundary claim on the pre-specified repeated switch experiment rather than wording that implies a strong generic MUC5B enrichment among all ZeBRA false negatives.

## Result directories retained

- `RESULTS/GENDRIVERS/validation_performance.csv`: fixed-split genomic-only performance.
- `RESULTS/COMBINED/validation_performance.csv`: fixed-split ZeBRA/combined performance.
- `RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/`: selected compact summaries from the final repeated-split/zedstat analysis.
- `RESULTS/033_INCREMENTAL_LOGISTIC_32_COHORT/`: selected summaries from the conditional incremental analysis.
- `RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/`: exploratory localization/provenance for the fixed regions.
- `RESULTS/MUC5B_ZEBRA_RESCUE_RULE/`: high/low-threshold rescue analyses.
- `RESULTS/FIXED_MUC5B_SWITCH_REGIONS/`: complete primary 100-split boundary analysis and figures.
- `RESULTS/MUC5B_ZEBRA_ERROR_STRATA/`: error-strata audit.

## Files deliberately not retained

- `032_COMBINED_CLASSIFIERS_REPEATED_SPLITS_ZEDSTAT_FIXED.ipynb`: exact duplicate of `32_COMBINED_CLASSIFIERS.ipynb`.
- `032_COMBINED_CLASSIFIERS_REPEATED_SPLITS.ipynb` and `032_COMBINED_CLASSIFIERS_REPEATED_SPLITS_ZEDSTAT.ipynb`: earlier stages superseded by `32` for the present 104-week repeated-split analysis.
- `033_COMBINED_CLASSIFIERS_REPEATED_SPLITS_EXTRATREES.ipynb`: alternative-model robustness experiment, not needed for current claims.
- `321_COMBINED_CLASSIFIERS_ISOTONIC_CALIBRATED.ipynb`: alternative calibration/modeling branch, not needed for current claims.
- `04_INCREMENTAL_GENETIC_LOGIT.ipynb`: superseded for current purposes by the repeated-split `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb`.
- `test_zebra_hybrid_roc_convex_hull_zedstat.py`: interesting ROC-envelope analysis, but not used in the current manuscript result.
- `test_zebra_predicts_ipf_prs.py`, `testinverse.py`, `testinverse_muc5b_extended.py`, and `testinverse_muc5b_stratified_by_fild.py`: exploratory/tangential to the present manuscript argument.
- all `.ipynb_checkpoints`, backup files, fitted model pickles, and large checkpoint/intermediate outputs.

## Reproduction order for the manuscript

For the main manuscript analysis, start from the included processed inputs and run:

1. `02_GEN_DRIVER_CLASSIFIERS.ipynb`
2. `03_COMBINED_CLASSIFIERS.ipynb`
3. `32_COMBINED_CLASSIFIERS.ipynb`
4. `33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb`
5. `test_local_zebra_genomics_predictive_curves.py` (provenance/exploration)
6. `test_muc5b_zebra_rescue_matched_fpr_lr.py` (high/low rescue support)
7. `test_fixed_muc5b_switch_regions.py` (primary boundary result)
8. `testinverse_muc5b_error_strata.py` only as an audit/sensitivity analysis

`01_GATHER_GENDRIVER_DATA.ipynb` is upstream provenance for reconstructing the processed genomic matrix, but downstream manuscript analyses can begin with `ILD_TOP_DRIVERS_DATA.csv`.

## Manuscript provenance correction

For the boundary-localization and fixed-switch sections, the code provenance is notebook **32**, not merely the earlier `03/032` lineage. The v0 boundary scripts explicitly reconstruct the notebook-32 cohort and use the 104-week ZeBRA prediction file. The manuscript wording should be audited so that these analyses point to notebook 32 / the 104-week cohort consistently.
