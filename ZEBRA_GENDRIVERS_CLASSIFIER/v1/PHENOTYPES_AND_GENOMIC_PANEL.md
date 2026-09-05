# v1 clinical targets and genomic panel

This note documents the three clinical prediction targets and the genomic feature panel used in the archived `ZEBRA_GENDRIVERS_CLASSIFIER/v1` ZeBRA + genomics analyses.

## Clinical prediction targets

The phenotype table is `REPHENOTYPES FOR IC.csv`. The three target columns used in notebook 03/032 are:

1. `FILD or FILA ADJUDICATED` — adjudicated fibrotic ILD or fibrotic interstitial lung abnormality (FILA).
2. `NONFIBROTIC ILD ADJUDICATED` — adjudicated nonfibrotic interstitial lung disease.
3. `NON FIBROTIC ILA ADJUDICATED` — adjudicated nonfibrotic interstitial lung abnormality.

The legacy notebook-03/032 target geometry maps `Y/y` to 1 and `N/n` to 0, converts the target to numeric, assigns missing target entries to 0, and then retains observations with non-missing ZeBRA `predicted_risk`.

Under this construction the ZeBRA-available cohort contains 12,825 observations for each target:

| Target | N | Positives | ZeBRA AUC on the archived notebook-03 split |
| --- | ---: | ---: | ---: |
| FILD or FILA ADJUDICATED | 12,825 | 254 | 0.8041 |
| NONFIBROTIC ILD ADJUDICATED | 12,825 | 27 | 0.8379 |
| NON FIBROTIC ILA ADJUDICATED | 12,825 | 74 | 0.8092 |

These are the stored adjudicated labels used by the analysis. This repository note does not infer additional adjudication criteria that are not encoded in the source files.

## Genomic panel construction

The v1 data-construction notebook `01_GATHER_GENDRIVER_DATA.ipynb` combines two sources of genomic variables:

- a pre-existing biological-driver panel already available in the genomic data object; and
- an expanded candidate-locus list from `MORE_LOCI.csv`.

`MORE_LOCI.csv` contains 3,113 candidate records. Matching this expanded list to the available genomic header identifies 148 genotype columns. The union of these matched columns with the pre-existing biological-driver set, followed by de-duplication, yields **171 distinct variant loci** in the archived v1 combined-model feature matrix.

Before phenotype/ZeBRA cohort restriction, this genomic matrix contains 19,651 rows.

Each locus is represented as a categorical state 0, 1, or 2 and was one-hot encoded for the discriminative LightGBM analyses. Thus a locus such as `rs35705950.1_G` appears in the model matrix as three indicator columns (`_0`, `_1`, `_2`). For generative/Quasinet analyses, the compact 0/1/2 locus representation is preferable to the one-hot representation.

The complete reconstructed locus list is stored in [`SNP_PANEL_USED.csv`](SNP_PANEL_USED.csv). It contains all 171 distinct locus/allele prefixes found in the archived v1 combined-model matrix, including the MUC5B promoter variant `rs35705950.1_G`.

## Archived fixed-split AUC comparison

The archived `RESULTS/COMBINED/validation_performance.csv` reports:

| Target | ZeBRA | Genomic | Combined |
| --- | ---: | ---: | ---: |
| FILD/FILA | 0.8041 | 0.6426 | 0.8000 |
| Nonfibrotic ILD | 0.8379 | 0.6290 | 0.8547 |
| Nonfibrotic ILA | 0.8092 | 0.5728 | 0.8060 |

The genomic-only rows in that archived file were evaluated on their stored genomic validation geometry, whereas ZeBRA and combined rows correspond to the target-specific ZeBRA-available validation geometry. The table should therefore be treated as descriptive rather than as a formal paired comparison.

## Boundary-localized MUC5B result

The repeated-split fixed-region experiment for the FILD/FILA target used MUC5B information in ZeBRA percentile windows 0–50, 65–75, and 85–95, while retaining calibrated ZeBRA outside those regions.

Across 100 repeated outer splits, at the nominal 5% FPR operating point:

| Metric | ZeBRA | Fixed regional MUC5B hybrid | Difference |
| --- | ---: | ---: | ---: |
| Global AUC | 0.8160 | 0.8135 | -0.0025 |
| Sensitivity | 0.4576 | 0.4969 | +0.0393 |
| Realized FPR | 0.0586 | 0.0591 | +0.0005 |
| LR+ | 7.89 | 8.46 | +0.57 |
| LR- | 0.576 | 0.535 | -0.042 |
| PPV | 0.1367 | 0.1453 | +0.0086 |

The corresponding median training-split MUC5B odds ratios were approximately 2.88 in P00–P50, 4.40 in P65–P75, and 4.07 in P85–P95.

This is an exploratory/post-hoc regional rule: the percentile windows were selected after inspection of local-information curves. Independent or fully prespecified validation is required for confirmatory inference.
