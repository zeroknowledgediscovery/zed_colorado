# ZeBRA--genomics manuscript

This directory contains the canonical theory + empirical manuscript for the ZeBRA--genomics likelihood-ratio paper.

## Canonical manuscript

Build the current manuscript from:

```text
zebra_genomics_boundary_theorem_empirical.tex
```

with:

```bash
cd tex
./build_full_manuscript.sh
```

or directly:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
    zebra_genomics_boundary_theorem_empirical.tex
```

The canonical master is deliberately modular. Its narrative is assembled from:

```text
manuscript_body_intro.tex
manuscript_body_theory.tex
manuscript_body_interpretation.tex
manuscript_body_empirical_a.tex
manuscript_body_empirical_b.tex
manuscript_body_end.tex
```

The older `zebra_genomics_boundary_theorem.tex` and `empirical_section_update.tex` are retained as historical sources and are **not** part of the canonical build.

## Figure provenance

All manuscript result figures are generated natively with PGFPlots/TikZ. The canonical figure wrapper is:

```text
empirical_figures.tex
```

which loads:

```text
empirical_figures_data.tex
empirical_figures_primary.tex
empirical_figures_lr.tex
```

`empirical_figures_data.tex` defines a single results root:

```text
../zebra_comp/RESULTS/
```

and every plotted value is read from the canonical analysis outputs under that tree at compile time. There is deliberately:

- no `tex/data/` directory,
- no copied manuscript-side CSV layer,
- no pre-rendered PNG/PDF result figure in the manuscript.

The main figure inputs are:

1. Global repeated-split discrimination
   - `032_REPEATED_SPLITS_03_TARGET_LOGIC/auc_distribution_summary.csv`
   - `032_REPEATED_SPLITS_03_TARGET_LOGIC/paired_auc_delta_summary.csv`
2. Paired incremental analysis
   - `033_INCREMENTAL_LOGISTIC_32_COHORT/paired_incremental_summary.csv`
3. Disease-stratified ZeBRA--MUC5B analysis
   - `MUC5B_STRATIFIED_BY_FILD/MUC5B_ZEBRA_STRATIFIED_SUMMARY.csv`
4. Local conditional-information analysis
   - `LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/FPR_DEFINED_LOCAL_INFORMATION_BANDS.csv`
   - `LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/LOCAL_ZEBRA_MUC5B_INFORMATION_WINDOWS.csv`
5. Matched-operational-FPR MUC5B rescue
   - `MUC5B_ZEBRA_RESCUE_RULE/REPEATED_SPLIT_MUC5B_RESCUE_SUMMARY.csv`
6. Stringent operating-point behavior
   - `032_REPEATED_SPLITS_03_TARGET_LOGIC/zedstat_operating_points_summary.csv`
7. Direct likelihood-ratio robustness analysis
   - `EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS/LR_ESTIMATOR_SENSITIVITY_GRID.csv`
   - `EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS/CLINICAL_TAIL_ROBUSTNESS.csv`
   - `EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS/ROBUSTNESS_RANGE_SUMMARY.csv`

Implementation-specific source mappings and analysis-script names are kept in `zebra_genomics_boundary_theorem_empirical_supplement.tex`, not in the scientific figure captions or main narrative.

## Current theory

The primary theoretical condition is stochastic rather than a literal universal bound. A1* requires the residual genomic log-likelihood contribution to remain uniformly tight as clinical history accumulates:

```text
sup_n P0(|R_G,n| > b_delta) <= delta.
```

The main theorem gives high-probability boundary localization and high-probability clinical-tail dominance. The earlier deterministic two-sided bound is retained as a special-case corollary. A separate proposition gives the genomic-only power bound under a marginal genomic likelihood-ratio bound.

The direct Colorado LR analysis estimates the score-level factorization `Lambda_Z,G = Lambda_Z Lambda_G|Z`; it provides empirical support for the theory's geometry but does not prove a universal full-genome or full-history bound.

## Publication layout

The PGFPlots sources use compact categorical labels, in-axis panel letters, increased spacing between multi-panel axes, and logarithmic FPR axes where necessary to prevent label collisions. Wide floats are permitted to move naturally under IEEEtran; forced float barriers are intentionally avoided to prevent large blank regions.

## Supplement

The computational/provenance supplement is:

```text
zebra_genomics_boundary_theorem_empirical_supplement.tex
```

It records exact result-file mappings, direct-LR estimator specifications, calibration diagnostics, and supporting ROC-hull information that are intentionally omitted from publication-facing captions and prose.
