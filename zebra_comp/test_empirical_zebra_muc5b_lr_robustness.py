#!/usr/bin/env python3
"""Estimator-sensitivity analysis for the empirical ZeBRA--MUC5B LR study.

Run from zebra_comp/ after test_empirical_zebra_muc5b_likelihood_ratios.py.
Uses exactly the same processed input files and cohort construction by importing
that script, but refits the cross-fitted spline/logistic LR estimator over a
small prespecified grid:

    spline knots = 4, 5, 6
    logistic C   = 0.3, 1.0, 3.0

The purpose is not model selection. It tests whether the manuscript-relevant
conclusions are stable to reasonable smoothing/regularization choices:
  * high-probability residual-genomic scale b_delta,
  * clinical-tail dominance at ZeBRA >=95, 97.5, 99 percentiles,
  * localization of MUC5B-induced fixed-LR decision flips.

Outputs are written to:
  RESULTS/EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS/

No manuscript result should be selected from this grid by choosing the most
favorable specification; the default 5-knot/C=1 analysis remains primary.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import test_empirical_zebra_muc5b_likelihood_ratios as lr


OUTDIR = Path("./RESULTS/EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS")
OUTDIR.mkdir(parents=True, exist_ok=True)

KNOT_GRID = [4, 5, 6]
C_GRID = [0.3, 1.0, 3.0]
DELTAS = [0.05, 0.01, 0.005]
TAIL_CUTOFFS = [95.0, 97.5, 99.0]
FPRS = [0.005, 0.01, 0.02, 0.05, 0.10]


def stochastic_point_summary(s):
    controls = s.loc[s.target == 0]
    x = controls.max_abs_log_LR_G_statewise.to_numpy(float)
    rows = []
    for delta in DELTAS:
        b = float(np.quantile(x, 1.0 - delta))
        rows.append({
            "delta": delta,
            "b_delta_statewise_abs_logLR": b,
            "exp_b_delta": float(np.exp(b)),
            "empirical_tail_fraction": float(np.mean(x > b)),
        })
    return pd.DataFrame(rows)


def tail_summary(s):
    base = lr.clinical_tail(s)
    return base.loc[base.ZeBRA_percentile_cutoff.isin(TAIL_CUTOFFS)].copy()


def decision_summary(s, stoch):
    old_fprs = lr.FPRS
    try:
        lr.FPRS = FPRS
        d, _ = lr.decision_localization(s, stoch)
    finally:
        lr.FPRS = old_fprs
    return d


def main():
    print("Loading cohort once using the primary LR script...")
    d = lr.load_data()
    print(
        f"N={len(d):,}; cases={int(d.target.sum()):,}; "
        f"MUC5B carriers={int(d.MUC5B_T_carrier.sum()):,}"
    )

    original_knots = lr.N_KNOTS
    original_c = lr.LOGISTIC_C

    grid_rows = []
    tail_rows = []
    flip_rows = []

    try:
        for knots in KNOT_GRID:
            for c_value in C_GRID:
                print("=" * 78)
                print(f"Running knots={knots}, C={c_value:g}")

                lr.N_KNOTS = knots
                lr.LOGISTIC_C = c_value

                py, q1, q0, _ = lr.crossfit(d)
                s, diag = lr.construct_lr(d, py, q1, q0)
                stoch = stochastic_point_summary(s)
                tail = tail_summary(s)
                flips = decision_summary(s, stoch)

                row = {
                    "n_knots": knots,
                    "logistic_C": c_value,
                    "N": len(s),
                    "cases": int(s.target.sum()),
                    "LR_Z_P0_mean": float(diag["LR_Z_P0_mean_after_normalization"]),
                    "LR_ZG_P0_mean": float(diag["LR_ZG_P0_mean"]),
                    "max_statewise_abs_logLR_G_controls": float(
                        s.loc[s.target == 0, "max_abs_log_LR_G_statewise"].max()
                    ),
                }

                for _, r in stoch.iterrows():
                    tag = str(r.delta).replace(".", "p")
                    row[f"b_delta_{tag}"] = float(r.b_delta_statewise_abs_logLR)
                    row[f"exp_b_delta_{tag}"] = float(r.exp_b_delta)

                grid_rows.append(row)

                tail.insert(0, "logistic_C", c_value)
                tail.insert(0, "n_knots", knots)
                tail_rows.append(tail)

                flips.insert(0, "logistic_C", c_value)
                flips.insert(0, "n_knots", knots)
                flip_rows.append(flips)

    finally:
        lr.N_KNOTS = original_knots
        lr.LOGISTIC_C = original_c

    grid = pd.DataFrame(grid_rows)
    tails = pd.concat(tail_rows, ignore_index=True)
    flips = pd.concat(flip_rows, ignore_index=True)

    grid.to_csv(OUTDIR / "LR_ESTIMATOR_SENSITIVITY_GRID.csv", index=False)
    tails.to_csv(OUTDIR / "CLINICAL_TAIL_ROBUSTNESS.csv", index=False)
    flips.to_csv(OUTDIR / "DECISION_FLIP_ROBUSTNESS.csv", index=False)

    # Compact across-specification ranges for the manuscript decision.
    range_rows = []

    for delta in DELTAS:
        tag = str(delta).replace(".", "p")
        col = f"b_delta_{tag}"
        vals = grid[col].to_numpy(float)
        range_rows.append({
            "metric": f"b_delta_{delta:g}",
            "min_across_9_specs": float(vals.min()),
            "median_across_9_specs": float(np.median(vals)),
            "max_across_9_specs": float(vals.max()),
        })

    for cutoff in TAIL_CUTOFFS:
        x = tails.loc[tails.ZeBRA_percentile_cutoff == cutoff]
        for metric in [
            "median_log_LR_Z",
            "q95_max_abs_log_LR_G_statewise",
            "fraction_log_LR_Z_exceeds_statewise_genomic_scale",
        ]:
            vals = x[metric].to_numpy(float)
            range_rows.append({
                "metric": f"tail_{cutoff:g}_{metric}",
                "min_across_9_specs": float(vals.min()),
                "median_across_9_specs": float(np.median(vals)),
                "max_across_9_specs": float(vals.max()),
            })

    for fpr in FPRS:
        x = flips.loc[np.isclose(flips.requested_clinical_FPR, fpr)]
        for metric in [
            "n_flips",
            "median_boundary_distance_flips",
            "q95_boundary_distance_flips",
        ]:
            vals = x[metric].to_numpy(float)
            range_rows.append({
                "metric": f"fpr_{100*fpr:g}pct_{metric}",
                "min_across_9_specs": float(np.nanmin(vals)),
                "median_across_9_specs": float(np.nanmedian(vals)),
                "max_across_9_specs": float(np.nanmax(vals)),
            })

    ranges = pd.DataFrame(range_rows)
    ranges.to_csv(OUTDIR / "ROBUSTNESS_RANGE_SUMMARY.csv", index=False)

    # Human-readable terminal/file summary.
    default = grid.loc[(grid.n_knots == 5) & np.isclose(grid.logistic_C, 1.0)].iloc[0]
    lines = [
        "EMPIRICAL LR ESTIMATOR ROBUSTNESS",
        "Grid: n_knots={4,5,6}, logistic C={0.3,1,3}; primary specification is 5/1.",
        "This is a sensitivity analysis, not a model-selection exercise.",
        "",
        "Primary specification:",
        f"  P0 mean LR_ZG = {default.LR_ZG_P0_mean:.6f}",
        f"  max control statewise |log LR_G| = {default.max_statewise_abs_logLR_G_controls:.4f}",
        "",
        "Across-specification ranges:",
    ]
    for _, r in ranges.iterrows():
        lines.append(
            f"  {r.metric}: min={r.min_across_9_specs:.4f}, "
            f"median={r.median_across_9_specs:.4f}, max={r.max_across_9_specs:.4f}"
        )
    text = "\n".join(lines) + "\n"
    (OUTDIR / "SUMMARY.txt").write_text(text)
    print(text)
    print("Wrote", OUTDIR.resolve())


if __name__ == "__main__":
    main()
