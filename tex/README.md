This directory contains the ZeBRA-genomics manuscript sources.

Primary empirical manuscript:

```bash
cd tex
latexmk -pdf -interaction=nonstopmode -halt-on-error zebra_genomics_boundary_theorem_empirical.tex
```

The empirical wrapper imports `zebra_genomics_boundary_theorem.tex` and adds the current Colorado results. The principal tables and plots are generated directly in LaTeX/PGFPlots so that the manuscript remains self-contained and publication-scalable.

The current empirical figures summarize:

- repeated-split global AUC for ZeBRA, genomics, and their global combination;
- paired incremental AUC from adding the genomic panel to ZeBRA;
- the apparent versus disease-stratified ZeBRA--MUC5B relationship;
- local conditional MUC5B information and matched-FPR decision-boundary rescue.

Wide tables and multi-panel figures use IEEE two-column `table*`/`figure*` floats to avoid column and margin overlap.
