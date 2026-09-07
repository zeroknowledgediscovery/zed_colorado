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

The manuscript deliberately uses the original `zebra_comp` result plots when they contain richer split-level or operating-point information than a reduced replot, while retaining PGFPlots for compact manuscript-native summaries. This makes the source traceable to the exact analyses that produced each result.

Current main empirical figures:

1. Notebook-32 repeated-split AUC boxplots:
   `../zebra_comp/RESULTS/032_REPEATED_SPLITS_03_TARGET_LOGIC/auc_distributions_three_targets.png`
2. Notebook-33 paired incremental AUC:
   PGFPlot in `empirical_figures.tex`, derived from
   `../zebra_comp/RESULTS/033_INCREMENTAL_LOGISTIC_32_COHORT/paired_incremental_summary.csv`
3. Apparent versus disease-stratified ZeBRA--MUC5B association:
   `../zebra_comp/RESULTS/ZEBRA_MUC5B_APPARENT_ASSOCIATION/ZEBRA_MUC5B_NONLINEAR_ASSOCIATION.png`
   and
   `../zebra_comp/RESULTS/MUC5B_STRATIFIED_BY_FILD/MUC5B_enrichment_by_FILD_stratum.png`
4. Local genomic information:
   `../zebra_comp/RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/LOCAL_INFORMATION_CROSSOVER.png`
   and
   `../zebra_comp/RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES/FPR_DEFINED_INFORMATION_BANDS.png`
5. Matched-FPR decision-boundary rescue:
   `../zebra_comp/RESULTS/MUC5B_ZEBRA_RESCUE_RULE/MATCHED_FPR_sensitivity_rescue_vs_zebra.png`
   `../zebra_comp/RESULTS/MUC5B_ZEBRA_RESCUE_RULE/MUC5B_RESCUE_delta_sensitivity_vs_delta_FPR.png`
   `../zebra_comp/RESULTS/MUC5B_ZEBRA_RESCUE_RULE/MATCHED_FPR_LRplus_rescue_vs_zebra.png`
6. Low-FPR ROC geometry:
   `../zebra_comp/RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL/ROC_ZEBRA_HYBRID_ZEDSTAT_CONVEX_HULL_LOW_FPR.png`

The empirical section also contains a theory--evidence map tying the figures/tables to the likelihood-ratio theorem. Wide tables and multi-panel figures use IEEE `table*`/`figure*` floats to avoid column and margin overlap.

The analysis pipeline that regenerates the corresponding `zebra_comp/RESULTS` folders is documented in `../zebra_comp/README.md` and run with `../zebra_comp/run_all.sh`.
