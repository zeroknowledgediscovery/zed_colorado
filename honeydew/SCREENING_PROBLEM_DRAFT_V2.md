# Screening problem draft v2: when does genomics add predictive value beyond longitudinal clinical history?

## Scientific setting

You are evaluating whether germline genomic information provides predictive information about future interstitial lung disease beyond what can already be extracted from longitudinal clinical history.

The central issue is not simply whether a clinical model or a genomic model can classify disease. Incremental genomic value is conditional on the quality of the clinical representation used as the baseline. A weak clinical model can make genomics appear strongly complementary even when a richer representation of the same longitudinal record captures most of the predictive information relevant to discrimination.

The supplied benchmarking data are de-identified and use privacy-preserving synthetic/derived representations. They contain no direct patient identifiers or re-identification keys and are intended for analysis without patient privacy concerns. The de-identification and synthetic construction do not make the data toy-like: the benchmark is designed to preserve deep biological and clinical mechanistic realism, including longitudinal disease trajectories, diagnosis/medication/procedure dependencies, genotype-phenotype relationships, temporal structure, phenotype heterogeneity, and higher-order interactions characteristic of complex biomedical systems.

The supplied data contain raw longitudinal electronic-health-record information, germline genotypes for a fixed SNP panel, a predefined diagnosis-code set for a broad ILD phenotype, an independently adjudicated FILD/FILA phenotype, and an authoritative patient-level train/test split with a patient-specific historical prediction time.

There are two outcome definitions.

1. **Diagnosis-defined ILD.** A patient is positive when a qualifying diagnosis in `ild_codes.csv` first occurs in the diagnosis stream under the supplied cohort design. `patient_split.csv` supplies the patient-specific prediction/index time `t0`. The cohort was constructed so that the qualifying diagnosis for cases is 364 days after `t0`. Clinical predictors may use only information in the 1,000-day observation window ending at `t0`, inclusive. Information after `t0` must not be used as a predictor.

2. **Adjudicated FILD/FILA.** `rephenotypes.csv` contains an independent adjudication. The first column is patient ID. The column `FILD or FILA ADJUDICATED` is the outcome: `Y`/`y` is positive, `N`/`n` is negative, and blank, `NA`, or other values are not valid adjudicated labels. For this phenotype, preserve the same patient-level train/test assignment and the same patient-specific `t0` supplied in `patient_split.csv`; restrict the existing split to patients with a valid adjudicated label and the required data rather than constructing a new split.

## Supplied files

- `clinical_data.tgz`: gzip-compressed tar archive containing one longitudinal patient JSON file. The records contain demographics and dated diagnosis, medication, and procedure histories.
- `genotypes.csv`: fixed germline genomic feature panel. Genotype features are represented as per-variant genotype-state columns with suffixes `_0`, `_1`, and `_2`. Complete genotype-state triplets may be transformed to another defensible representation using training data only. The panel includes rs35705950 at the MUC5B locus.
- `rephenotypes.csv`: adjudicated FILD/FILA phenotype table described above.
- `ild_codes.csv`: exact diagnosis codes defining the broad diagnosis-based ILD outcome.
- `patient_split.csv`: authoritative patient-level evaluation manifest with the fields `patient_id`, `split`, and `t0`. Use these assignments exactly; do not resplit the cohort.

The genomic panel was fixed before this task. This is a prediction and incremental-information problem, not a genome-wide discovery exercise. Conclusions must follow from the supplied data rather than prior knowledge about individual variants.

## Practitioner decision

A research group has access to longitudinal EHR history for prediction and is deciding whether germline genotype information adds enough predictive value to justify incorporating it into the prediction system.

The key decision is therefore not whether genomics can predict disease in isolation. It is whether any apparent genomic gain remains after the clinical history has been modeled adequately, and whether the genomic channel can still make reproducible contributions inside a combined model even when global held-out discrimination changes little.

The analysis must distinguish at least three concepts:

1. standalone genomic discrimination;
2. incremental held-out predictive performance beyond clinical history;
3. evidence that genomic variables are actually used or remain conditionally informative inside a combined model.

These quantities need not agree.

## Required analysis

For each phenotype, build and compare clinical-only, genomics-only, and combined predictors on the supplied held-out patients. The clinical part of the workflow must explicitly address **clinical baseline adequacy** rather than treating one arbitrary feature representation as definitive.

At minimum:

- construct a simple, conventional clinical baseline from the longitudinal record;
- construct a stronger clinical representation that makes substantive use of the available longitudinal information, including diagnoses, medications, procedures, demographics, and temporal structure where useful;
- select and tune both clinical approaches using training patients only;
- construct a genomics-only predictor using the fixed SNP panel;
- add genomics separately to the simple and stronger clinical predictors;
- evaluate all comparisons on the same eligible held-out patients for the phenotype being analyzed;
- quantify held-out ROC AUC and the paired incremental change in AUC when genomics is added to each clinical baseline, including uncertainty;
- assess whether conclusions about incremental genomic value change as the clinical model becomes stronger;
- use a model-appropriate feature-attribution or conditional-contribution analysis to determine whether genomic variables, including rs35705950 when supported by the data, remain reproducible contributors inside the stronger combined model;
- explicitly distinguish persistent genomic attribution/conditional contribution from actual improvement in global predictive discrimination.

Do not choose a stronger clinical model, genomic representation, fusion strategy, feature set, or hyperparameter because it performs better on the held-out test outcomes. All such decisions must be made using training data only. The held-out test set is for final evaluation of the frozen analysis choices.

Use `patient_split.csv` as the single authoritative outer split for all modalities and both outcomes. When restricting to the adjudicated phenotype, preserve each patient's supplied train/test membership. A patient must never appear in both train and test data. Clinical information after `t0` is unavailable to the prospective predictor.

The task is not to reproduce any particular proprietary or named clinical model. The choice of statistical and machine-learning methods is yours, but the clinical representation must be rich enough that conclusions about genomics are not merely artifacts of an obviously underfit diagnosis-only baseline.

## Screening deliverable

Provide a reproducible analysis and a concise machine-readable summary containing, separately for each phenotype:

- held-out AUC for the simple clinical baseline;
- held-out AUC for the stronger longitudinal clinical model;
- held-out genomics-only AUC;
- held-out combined AUC when genomics is added to the simple clinical baseline;
- held-out combined AUC when genomics is added to the stronger clinical model;
- paired incremental AUC estimates and uncertainty intervals for genomics relative to both clinical baselines;
- evidence on whether important genomic variables remain reproducible contributors in the stronger combined model;
- a final decision about whether genomics provides practically meaningful incremental predictive value once longitudinal clinical history has been modeled adequately;
- a short scientific explanation of any difference between predictive performance and genomic attribution/conditional contribution, and of any dependence on phenotype definition.

ROC AUC is the primary common discrimination metric. Additional diagnostics are welcome when they materially support the practitioner decision. Do not infer incremental predictive value from standalone genomic performance, marginal association, statistical significance, or feature attribution alone.

For this screening run, do not assume that stronger clinical modeling must either eliminate or preserve the apparent genomic gain. Let the supplied data determine the result.