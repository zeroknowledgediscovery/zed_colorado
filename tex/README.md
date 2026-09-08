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

The canonical master is modular. Its narrative is assembled from:

```text
manuscript_body_intro.tex
manuscript_body_theory.tex
manuscript_body_interpretation.tex
manuscript_body_empirical_a.tex
manuscript_body_empirical_b.tex
manuscript_body_end.tex
```

The older `zebra_genomics_boundary_theorem.tex` and `empirical_section_update.tex` are historical sources and are **not** part of the canonical build.

## Figure provenance

All manuscript result figures are generated natively with PGFPlots/TikZ. The canonical wrapper is:

```text
empirical_figures.tex
```

which loads exactly:

```text
empirical_figures_data.tex
empirical_figures_primary.tex
empirical_figures_lr.tex
```

There are no secondary figure overrides. `empirical_figures_primary.tex` contains Figures 1, 2, 3, and 5; `empirical_figures_lr.tex` contains Figure 4. `empirical_figures_data.tex` defines a single results root:

```text
../zebra_comp/RESULTS/
```

Every plotted value is read from canonical analysis output under that tree at compile time. There is deliberately no `tex/data/` directory and no pre-rendered raster result figure in the manuscript.

### Current main-figure mapping

1. **Figure 1: Global discrimination and incremental genomic value**
   - `032_REPEATED_SPLITS_03_TARGET_LOGIC/auc_distribution_summary.csv`
   - `032_REPEATED_SPLITS_03_TARGET_LOGIC/paired_auc_delta_summary.csv`
   - `033_INCREMENTAL_LOGISTIC_32_COHORT/paired_incremental_summary.csv`
2. **Figure 2: ZeBRA and MUC5B as distinct predictive channels**
   - `MUC5B_STRATIFIED_BY_FILD/MUC5B_ZEBRA_STRATIFIED_SUMMARY.csv`
3. **Figure 3: Localized complementarity across the ZeBRA risk axis**
   - `LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/FPR_DEFINED_LOCAL_INFORMATION_BANDS.csv`
   - `LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/LOCAL_ZEBRA_MUC5B_INFORMATION_WINDOWS.csv`
4. **Figure 4: Direct empirical likelihood-ratio support**
   - `EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS/LR_ESTIMATOR_SENSITIVITY_GRID.csv`
   - `EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS/CLINICAL_TAIL_ROBUSTNESS.csv`
   - `EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS/ROBUSTNESS_RANGE_SUMMARY.csv`
5. **Figure 5: Decision-level consequences**
   - `MUC5B_ZEBRA_RESCUE_RULE/REPEATED_SPLIT_MUC5B_RESCUE_SUMMARY.csv`
   - `MUC5B_ZEBRA_RESCUE_RULE/REPEATED_SPLIT_ADJUDICATED_MATCHED_FPR_SUMMARY.csv`
   - `032_REPEATED_SPLITS_03_TARGET_LOGIC/zedstat_operating_points_summary.csv`

Implementation-specific script names and detailed provenance are kept in `zebra_genomics_boundary_theorem_empirical_supplement.tex`, not in publication-facing captions or narrative.

## Current theory

The theory separates two residual-genomic regularity conditions. At a fixed clinical-history stage `n`, the finite-stage stochastic scale `b_{delta,n}` is sufficient for the boundary-localization result. A stronger optional condition, A1-U, requires a single `b_delta` that works uniformly over clinical-history length and is used only for the sequential tail-dominance statement. A2 remains an explicit clinical-tail assumption rather than an empirical claim of mathematical unboundedness.

The direct Colorado analysis instantiates the factorization at the observed score level, `Lambda_Z,G = Lambda_Z Lambda_G|Z`. It supports finite-stage boundary geometry for the measured MUC5B channel but does not prove A1-U, A2, a full-history sufficiency relation, or a universal genome-wide bound.

## Publication layout

All multi-panel figures use lower-case panel labels `(a)`, `(b)`, etc. positioned outside the upper-left of each panel. Legends are placed above axes where needed to avoid covering data. Categorical labels are compact and horizontal, and logarithmic FPR axes are used where appropriate. Figure 3 deliberately omits nominal significance-threshold lines because the local-window analyses are descriptive rather than confirmatory.

## Supplement

The computational/provenance supplement is:

```text
zebra_genomics_boundary_theorem_empirical_supplement.tex
```

It records exact result-file mappings, direct-LR estimator specifications, calibration diagnostics, exploratory local-test details, the adjudicated-negative matched-comparator analysis, and supporting ROC-hull information that are intentionally omitted from publication-facing captions and prose.
