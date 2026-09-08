# ILD/MUC5B reference validation

Validation date: 2026-09-08

Validated generalized-pipeline code baseline:

`8c5db26623c3bc84ac2b0d5bae9bd109752e5dc8`

Reference: frozen Colorado ILD/MUC5B analysis under `../zebra_comp/RESULTS`.

Command used:

```bash
cd zebra_comp_general
./ild_reference.sh --run full
```

## Result

`PASS=38, WARN=0, FAIL=0, INFO=2`

The generalized pipeline therefore reproduces all parity-critical manuscript analyses. The two informational rows are the global nonlinear genomics-only and combined-model AUCs; these are intentionally not required to be numerically identical because the generalized nonlinear model path does not reproduce notebook 32's historical model-search implementation.

## Cohort invariants

- Full scored genomic cohort: 12,825
- Full FILD/FILA cases: 254
- MUC5B-called cohort: 12,766
- MUC5B-called FILD/FILA cases: 252

## Headline parity checks

| Quantity | Frozen reference | Generalized | Difference |
|---|---:|---:|---:|
| Mean ZeBRA AUC | 0.814152 | 0.814573 | +0.000421 |
| Mean incremental delta AUC | +0.0004349 | -0.0001172 | -0.0005521 |
| Mean incremental delta AP | -0.0311057 | -0.0314756 | -0.0003699 |
| Mean incremental delta Brier | 0.000312378 | 0.000315556 | +0.00000318 |
| Mean incremental delta log loss | 0.00467993 | 0.00490429 | +0.00022436 |
| b_hat(0.05) | 0.869869 | 0.869869 | 0 |
| b_hat(0.01) | 0.909919 | 0.909919 | 0 |
| b_hat(0.005) | 0.924098 | 0.924098 | 0 |
| E0[LR_Z] | 1.000000 | 1.000000 | 0 |
| E0[LR_ZG] | 1.01487 | 1.01487 | 0 |

All four control-FPR local-information bands reproduced exactly for the reported LRT, ZeBRA AUC, and MUC5B AUC quantities. The complete 3x3 direct-LR estimator-sensitivity grid reproduced exactly for the checked b_delta quantities.

The six matched-FPR rescue configurations all passed the prespecified parity tolerances; sensitivity differences were below 0.003 and FPR differences were on the order of 1e-4.

## Interpretation

This validation establishes the generalized implementation as the reusable reference for the manuscript-critical analysis chain:

1. ZeBRA discrimination;
2. paired incremental genomic value;
3. local conditional genomic information;
4. direct score-level likelihood-ratio estimation and finite-stage b_delta;
5. estimator robustness;
6. matched-FPR boundary rescue.

The frozen `../zebra_comp/` directory should still be retained as the provenance snapshot for the current ILD manuscript. The generalized pipeline is the version to use for new disease/gene combinations.
