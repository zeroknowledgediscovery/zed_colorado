# ZeBRA--genomics manuscript

This directory contains the manuscript sources for the ZeBRA--genomics likelihood-ratio paper.

## Canonical manuscript

The current complete theory + empirical manuscript is built from:

```text
zebra_genomics_boundary_theorem_empirical.tex
```

Build with:

```bash
cd tex
./build_full_manuscript.sh
```

or directly:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
    zebra_genomics_boundary_theorem_empirical.tex
```

The canonical manuscript uses `zebra_genomics_boundary_theorem.tex` as the theory core, `empirical_section_update.tex` for the current Colorado evidence section, and `empirical_figures.tex` for all result figures.

## Figure provenance

All manuscript result figures are generated natively with PGFPlots/TikZ. The figure code reads the canonical analysis outputs directly from:

```text
../zebra_comp/RESULTS/
```

There is deliberately no manuscript-side `data/` directory, no copied CSV layer, and no dependency on pre-rendered PNG result figures.

Current figure inputs are:

1. Global repeated-split discrimination
   - `../zebra_comp/RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/auc_distribution_summary.csv`
   - `../zebra_comp/RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/paired_auc_delta_summary.csv`
2. Paired incremental analysis
   - `../zebra_comp/RESULTS/033_INCREMENTAL_LOGISTIC_32_COHORT/paired_incremental_summary.csv`
3. Disease-stratified ZeBRA--MUC5B analysis
   - `../zebra_comp/RESULTS/MUC5B_STRATIFIED_BY_FILD/MUC5B_ZEBRA_STRATIFIED_SUMMARY.csv`
4. Local conditional-information analysis
   - `../zebra_comp/RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/FPR_DEFINED_LOCAL_INFORMATION_BANDS.csv`
   - `../zebra_comp/RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/LOCAL_ZEBRA_MUC5B_INFORMATION_WINDOWS.csv`
5. Matched-operational-FPR MUC5B rescue
   - `../zebra_comp/RESULTS/MUC5B_ZEBRA_RESCUE_RULE/REPEATED_SPLIT_MUC5B_RESCUE_SUMMARY.csv`
6. Stringent operating-point behavior
   - `../zebra_comp/RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/zedstat_operating_points_summary.csv`

The local-information figure therefore reproduces from canonical CSV outputs the FPR-band conditional LRTs, moving-window conditional-information curves, MUC5B carrier/non-carrier prevalence curves, and local adjusted odds-ratio curves.

## Theory/empirical distinction

The theorem is conditional on two explicit regularity assumptions: bounded residual genomic likelihood-ratio evidence (A1) and an extended clinical likelihood-ratio tail (A2). The canonical wrapper makes the measure-theoretic convention explicit: likelihood ratios are Radon--Nikodym derivatives, the A1 bounds are interpreted almost surely, and A2 is an essential-unboundedness statement. The Colorado analyses are correspondence tests, not proofs of A1 or A2.

In the local empirical analysis, "conditional LRT" means a nested-model likelihood-ratio/deviance statistic for adding one predictor after the other. It is evidence of conditional predictive contribution and is not an estimate of the individual-level factor `Lambda_{G|C}` in the theorem.

The historical ZeBRA field `predicted_risk` is also not treated as an absolute posterior probability in the Colorado analysis because it is not calibrated as such in this cohort.

## Layout

The canonical IEEE build uses permissive two-column float settings and intentionally avoids forced float barriers. This keeps the wide PGFPlots figures and tables close to their discussion without the large blank regions produced by earlier forced float placement. The four-panel local-information figure is slightly reduced in height at build time to avoid an overfull two-column float, while retaining all canonical RESULTS data.

The analysis pipeline that regenerates the corresponding `zebra_comp/RESULTS` folders is documented in `../zebra_comp/README.md` and run with `../zebra_comp/run_all.sh`.
