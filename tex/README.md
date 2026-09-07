This directory contains the ZeBRA--genomics manuscript sources.

Canonical full manuscript:

```bash
cd tex
latexmk -pdf -interaction=nonstopmode -halt-on-error \
    zebra_genomics_boundary_theorem_empirical.tex
```

`zebra_genomics_boundary_theorem_empirical.tex` is the complete theory+empirical manuscript. It imports the theoretical core from `zebra_genomics_boundary_theorem.tex`, injects the Colorado empirical section, and adds the supplementary provenance/reproducibility material.

Empirical figure source is modularized in:

```text
empirical_figures.tex
```

All manuscript result figures are generated natively with PGFPlots/TikZ. The figure code reads the canonical analysis outputs directly from:

```text
../zebra_comp/RESULTS/
```

There is deliberately no manuscript-side `data/` directory and no copied CSV layer. The manuscript figures also do not depend on pre-rendered PNG result plots. This keeps the plotted values tied directly to the outputs of the analysis pipeline.

Current empirical figure inputs include:

1. Global repeated-split discrimination:
   `../zebra_comp/RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/auc_distribution_summary.csv`
   `../zebra_comp/RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/paired_auc_delta_summary.csv`
2. Paired incremental analysis:
   `../zebra_comp/RESULTS/033_INCREMENTAL_LOGISTIC_32_COHORT/paired_incremental_summary.csv`
3. Disease-stratified ZeBRA--MUC5B analysis:
   `../zebra_comp/RESULTS/MUC5B_STRATIFIED_BY_FILD/MUC5B_ZEBRA_STRATIFIED_SUMMARY.csv`
4. Local genomic information:
   `../zebra_comp/RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/FPR_DEFINED_LOCAL_INFORMATION_BANDS.csv`
5. Matched-operational-FPR MUC5B rescue:
   `../zebra_comp/RESULTS/MUC5B_ZEBRA_RESCUE_RULE/REPEATED_SPLIT_MUC5B_RESCUE_SUMMARY.csv`
6. Low-FPR ROC geometry and stringent-tail operating behavior:
   `../zebra_comp/RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL/ROC_ZEBRA.csv`
   `../zebra_comp/RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL/ROC_HYBRID.csv`
   `../zebra_comp/RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL/ROC_UNION_ZEDSTAT_CONVEX_HULL.csv`
   `../zebra_comp/RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/zedstat_operating_points_summary.csv`

The empirical section contains a theory--evidence map tying the figures and tables to the likelihood-ratio theorem. Wide tables and multi-panel figures use IEEE `table*`/`figure*` floats.

The analysis pipeline that regenerates the corresponding `zebra_comp/RESULTS` folders is documented in `../zebra_comp/README.md` and run with `../zebra_comp/run_all.sh`.
