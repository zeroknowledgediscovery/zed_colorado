#!/usr/bin/env python3
"""
Generate the missing fair comparator for the MUC5B rescue analysis:

    FPR(rescue | explicitly adjudicated negatives)
      versus
    FPR(matched-ZeBRA | explicitly adjudicated negatives)

This script is intentionally a POST-PROCESSOR for
test_muc5b_zebra_rescue_matched_fpr_lr.py.

It reuses the exact split seeds and matched-ZeBRA thresholds already written to:
    RESULTS/MUC5B_ZEBRA_RESCUE_RULE/REPEATED_SPLIT_MUC5B_RESCUE_RULE.csv

For each repeated outer split it reconstructs the untouched test set, identifies
the explicitly adjudicated-negative subjects, and applies:
  1. the original stringent ZeBRA rule,
  2. the MUC5B rescue rule, and
  3. the training-FPR-matched ZeBRA-only comparator.

No threshold is fitted on the adjudicated-negative test subset.

Run from zebra_comp/:
    python test_muc5b_rescue_adjudicated_matched_fpr.py

Outputs:
    RESULTS/MUC5B_ZEBRA_RESCUE_RULE/
        REPEATED_SPLIT_ADJUDICATED_MATCHED_FPR.csv
        REPEATED_SPLIT_ADJUDICATED_MATCHED_FPR_SUMMARY.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


# ============================================================
# Configuration -- must match the original rescue script
# ============================================================

BASE_DATA_FILE = "ILD_TOP_DRIVERS_DATA.csv"
TARGET_FILE = "REPHENOTYPES FOR IC.csv"
PRED_FILE = "PREDICTIONS_104W_PRED_WINDOW.parquet"

TARGET_NAME = "FILD or FILA ADJUDICATED"

GG_COLUMN = "rs35705950.1_G_2"
GT_COLUMN = "rs35705950.1_G_1"
TT_COLUMN = "rs35705950.1_G_0"

TRAIN_SIZE = 0.40

OUTDIR = Path("./RESULTS/MUC5B_ZEBRA_RESCUE_RULE")
SOURCE_RESULTS = OUTDIR / "REPEATED_SPLIT_MUC5B_RESCUE_RULE.csv"

DETAIL_OUT = OUTDIR / "REPEATED_SPLIT_ADJUDICATED_MATCHED_FPR.csv"
SUMMARY_OUT = OUTDIR / "REPEATED_SPLIT_ADJUDICATED_MATCHED_FPR_SUMMARY.csv"

RECONSTRUCTION_TOL = 1e-12


# ============================================================
# Helpers
# ============================================================

def safe_fpr(positive, explicit_negative):
    positive = np.asarray(positive, dtype=bool)
    explicit_negative = np.asarray(explicit_negative, dtype=bool)

    n = int(explicit_negative.sum())
    if n == 0:
        return np.nan, 0, 0

    fp = int(np.sum(positive & explicit_negative))
    return fp / n, fp, n


def summarize(df):
    rows = []

    group_cols = [
        "main_target_fpr",
        "lower_target_fpr",
    ]

    value_cols = [
        "explicit_negative_N",
        "baseline_explicit_negative_FPR",
        "rescued_explicit_negative_FPR",
        "matched_zebra_explicit_negative_FPR",
        "rescue_minus_baseline_explicit_negative_FPR",
        "matched_minus_baseline_explicit_negative_FPR",
        "rescue_minus_matched_explicit_negative_FPR",
        "baseline_explicit_negative_FP",
        "rescued_explicit_negative_FP",
        "matched_zebra_explicit_negative_FP",
        "rescue_minus_matched_explicit_negative_FP",
    ]

    for keys, sub in df.groupby(group_cols, sort=True):
        row = dict(zip(group_cols, keys))
        row["n_repeats"] = int(len(sub))

        for col in value_cols:
            vals = pd.to_numeric(sub[col], errors="coerce").dropna().to_numpy(float)
            if vals.size == 0:
                continue

            row[f"{col}_mean"] = float(np.mean(vals))
            row[f"{col}_median"] = float(np.median(vals))
            row[f"{col}_q025"] = float(np.quantile(vals, 0.025))
            row[f"{col}_q975"] = float(np.quantile(vals, 0.975))

            if col.startswith("rescue_minus_"):
                row[f"{col}_fraction_gt_zero"] = float(np.mean(vals > 0))
                row[f"{col}_fraction_lt_zero"] = float(np.mean(vals < 0))
                row[f"{col}_fraction_eq_zero"] = float(np.mean(vals == 0))

        rows.append(row)

    return pd.DataFrame(rows)


# ============================================================
# Load original repeated-split results
# ============================================================

if not SOURCE_RESULTS.exists():
    raise FileNotFoundError(
        f"Missing {SOURCE_RESULTS}\n"
        "Run test_muc5b_zebra_rescue_matched_fpr_lr.py first."
    )

results = pd.read_csv(SOURCE_RESULTS)

required_result_cols = {
    "repeat",
    "split_seed",
    "main_target_fpr",
    "lower_target_fpr",
    "main_threshold",
    "lower_threshold",
    "matched_zebra_threshold",
    "baseline_explicit_negative_FPR",
    "rescued_explicit_negative_FPR",
}

missing = sorted(required_result_cols - set(results.columns))
if missing:
    raise ValueError(
        "The existing repeated-split result file is missing required columns:\n"
        + "\n".join(missing)
    )


# ============================================================
# Reconstruct EXACT notebook-32 / rescue-script cohort
# ============================================================

print("Loading cohort...")

base_data = pd.read_csv(BASE_DATA_FILE).drop(
    columns=["target"],
    errors="ignore",
)

targets = (
    pd.read_csv(TARGET_FILE)
    .replace({
        "N": 0,
        "Y": 1,
        "n": 0,
        "y": 1,
        "na": np.nan,
    })
    .rename(columns={"arb_person_id": "patient_id"})
)

preds = pd.read_parquet(PRED_FILE)

data = (
    base_data
    .merge(
        targets[
            [
                "patient_id",
                TARGET_NAME,
            ]
        ].rename(columns={TARGET_NAME: "target_original"}),
        on="patient_id",
        how="left",
    )
    .merge(
        preds[
            [
                "patient_id",
                "predicted_risk",
            ]
        ],
        on="patient_id",
        how="left",
    )
)

data["target_was_observed"] = data["target_original"].notnull()

data["target"] = (
    pd.to_numeric(data["target_original"], errors="coerce")
    .fillna(0)
    .astype(int)
)

# Exact ZeBRA availability filter from the original analysis.
data = data.loc[data["predicted_risk"].notnull()].copy()

required_genotype_cols = [
    GG_COLUMN,
    GT_COLUMN,
    TT_COLUMN,
]

missing = [c for c in required_genotype_cols if c not in data.columns]
if missing:
    raise ValueError(
        "Missing expected MUC5B columns:\n" + "\n".join(missing)
    )

for col in required_genotype_cols:
    data[col] = pd.to_numeric(data[col], errors="coerce")

dummy_sum = data[required_genotype_cols].sum(axis=1)

if (dummy_sum > 1).any():
    raise ValueError("Some MUC5B rows are multi-hot.")

# Exact called-genotype restriction from the original rescue script.
data = data.loc[dummy_sum == 1].copy()

data["MUC5B_T_carrier"] = (
    (data[GT_COLUMN].astype(int) == 1)
    | (data[TT_COLUMN].astype(int) == 1)
).astype(int)


# Arrays MUST preserve this row order to reproduce train_test_split.
indices = np.arange(len(data))

y_all = data["target"].to_numpy(dtype=int)
z_all = data["predicted_risk"].to_numpy(dtype=float)
carrier_all = data["MUC5B_T_carrier"].to_numpy(dtype=int)
observed_target_all = data["target_was_observed"].to_numpy(dtype=bool)


# ============================================================
# Reconstruct each split and calculate the missing comparator
# ============================================================

rows = []

for repeat, repeat_rows in results.groupby("repeat", sort=True):
    split_seeds = repeat_rows["split_seed"].dropna().unique()

    if len(split_seeds) != 1:
        raise ValueError(
            f"Repeat {repeat} has {len(split_seeds)} split seeds; expected exactly 1."
        )

    split_seed = int(split_seeds[0])

    _, test_idx = train_test_split(
        indices,
        train_size=TRAIN_SIZE,
        stratify=y_all,
        random_state=split_seed,
    )

    y_test = y_all[test_idx]
    z_test = z_all[test_idx]
    c_test = carrier_all[test_idx]
    observed_test = observed_target_all[test_idx]

    # "Explicitly adjudicated negative" means:
    # observed FILD/FILA target AND observed value == 0.
    explicit_negative = (y_test == 0) & observed_test

    for _, r in repeat_rows.iterrows():
        t_main = float(r["main_threshold"])
        t_lower = float(r["lower_threshold"])
        t_matched = float(r["matched_zebra_threshold"])

        baseline_positive = z_test >= t_main

        rescue_band = (
            (z_test >= t_lower)
            & (z_test < t_main)
        )
        rescued_positive = (
            baseline_positive
            | (
                rescue_band
                & (c_test == 1)
            )
        )

        # IMPORTANT:
        # This threshold was learned on TRAINING negatives in the
        # original rescue script. We merely apply it to the untouched
        # test set, including the explicitly adjudicated-negative subset.
        matched_zebra_positive = z_test >= t_matched

        baseline_fpr, baseline_fp, explicit_n = safe_fpr(
            baseline_positive,
            explicit_negative,
        )
        rescued_fpr, rescued_fp, _ = safe_fpr(
            rescued_positive,
            explicit_negative,
        )
        matched_fpr, matched_fp, _ = safe_fpr(
            matched_zebra_positive,
            explicit_negative,
        )

        # Sanity check: reconstruction should exactly recover the
        # baseline and rescue adjudicated-negative FPRs already saved
        # by the original script.
        saved_baseline = float(r["baseline_explicit_negative_FPR"])
        saved_rescued = float(r["rescued_explicit_negative_FPR"])

        if np.isfinite(saved_baseline) and np.isfinite(baseline_fpr):
            if abs(saved_baseline - baseline_fpr) > RECONSTRUCTION_TOL:
                raise RuntimeError(
                    f"Repeat {repeat}: baseline reconstruction mismatch: "
                    f"{saved_baseline} vs {baseline_fpr}. "
                    "Check package versions, row ordering, or source data."
                )

        if np.isfinite(saved_rescued) and np.isfinite(rescued_fpr):
            if abs(saved_rescued - rescued_fpr) > RECONSTRUCTION_TOL:
                raise RuntimeError(
                    f"Repeat {repeat}: rescue reconstruction mismatch: "
                    f"{saved_rescued} vs {rescued_fpr}. "
                    "Check package versions, row ordering, or source data."
                )

        rows.append({
            "repeat": int(repeat),
            "split_seed": split_seed,
            "main_target_fpr": float(r["main_target_fpr"]),
            "lower_target_fpr": float(r["lower_target_fpr"]),
            "main_threshold": t_main,
            "lower_threshold": t_lower,
            "matched_zebra_threshold": t_matched,

            "explicit_negative_N": explicit_n,

            "baseline_explicit_negative_FP": baseline_fp,
            "rescued_explicit_negative_FP": rescued_fp,
            "matched_zebra_explicit_negative_FP": matched_fp,

            "baseline_explicit_negative_FPR": baseline_fpr,
            "rescued_explicit_negative_FPR": rescued_fpr,
            "matched_zebra_explicit_negative_FPR": matched_fpr,

            "rescue_minus_baseline_explicit_negative_FPR":
                rescued_fpr - baseline_fpr,

            "matched_minus_baseline_explicit_negative_FPR":
                matched_fpr - baseline_fpr,

            # THIS IS THE MISSING FAIR COMPARATOR.
            "rescue_minus_matched_explicit_negative_FPR":
                rescued_fpr - matched_fpr,

            "rescue_minus_matched_explicit_negative_FP":
                rescued_fp - matched_fp,
        })


detail = pd.DataFrame(rows)

OUTDIR.mkdir(parents=True, exist_ok=True)
detail.to_csv(DETAIL_OUT, index=False)

summary = summarize(detail)
summary.to_csv(SUMMARY_OUT, index=False)


# ============================================================
# Report the manuscript-relevant quantities
# ============================================================

print()
print("=" * 110)
print("EXPLICITLY ADJUDICATED NEGATIVES: RESCUE VS FAIR MATCHED-ZEBRA COMPARATOR")
print("=" * 110)

display_cols = [
    "main_target_fpr",
    "lower_target_fpr",
    "n_repeats",

    "baseline_explicit_negative_FPR_mean",
    "rescued_explicit_negative_FPR_mean",
    "matched_zebra_explicit_negative_FPR_mean",

    "rescue_minus_baseline_explicit_negative_FPR_mean",
    "rescue_minus_matched_explicit_negative_FPR_mean",
    "rescue_minus_matched_explicit_negative_FPR_q025",
    "rescue_minus_matched_explicit_negative_FPR_q975",
    "rescue_minus_matched_explicit_negative_FPR_fraction_gt_zero",
]

display_cols = [c for c in display_cols if c in summary.columns]

print(summary[display_cols].to_string(index=False))

print()
print(f"Wrote: {DETAIL_OUT}")
print(f"Wrote: {SUMMARY_OUT}")
print()
print(
    "Primary manuscript quantity: "
    "rescue_minus_matched_explicit_negative_FPR_mean"
)
