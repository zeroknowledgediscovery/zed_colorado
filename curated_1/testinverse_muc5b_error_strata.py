#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
from sklearn.metrics import roc_auc_score


# ============================================================
# Configuration
# ============================================================

BASE_DATA_FILE = "ILD_TOP_DRIVERS_DATA.csv"
TARGET_FILE = "REPHENOTYPES FOR IC.csv"
PRED_FILE = "PREDICTIONS_104W_PRED_WINDOW.parquet"

TARGET_NAME = "FILD or FILA ADJUDICATED"

GG_COLUMN = "rs35705950.1_G_2"
GT_COLUMN = "rs35705950.1_G_1"
TT_COLUMN = "rs35705950.1_G_0"

# Define ZeBRA operating thresholds by the FPR achieved among
# notebook-32 target==0 patients.
TARGET_FPRS = [
    0.005,
    0.01,
    0.02,
    0.05,
    0.10,
]

OUTDIR = Path(
    "./RESULTS/MUC5B_ZEBRA_ERROR_STRATA"
)
OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Helpers
# ============================================================

def odds_ratio_ci(a, b, c, d):
    """
    2x2 table:
                     carrier  noncarrier
       group A          a         b
       group B          c         d

    Returns OR(A vs B) and Wald 95% CI.
    Applies Haldane-Anscombe +0.5 correction if needed.
    """
    cells = np.array(
        [a, b, c, d],
        dtype=float,
    )

    if np.any(cells == 0):
        cells += 0.5

    a, b, c, d = cells

    OR = (a * d) / (b * c)

    se = np.sqrt(
        1 / a
        + 1 / b
        + 1 / c
        + 1 / d
    )

    z = 1.959963984540054

    lo = np.exp(
        np.log(OR)
        - z * se
    )

    hi = np.exp(
        np.log(OR)
        + z * se
    )

    return (
        float(OR),
        float(lo),
        float(hi),
    )


def threshold_for_target_fpr(
    negative_scores,
    target_fpr,
):
    """
    Choose the lowest score threshold giving approximately the
    requested upper-tail FPR among negatives.

    Because many ZeBRA scores may tie at the top end, the achieved
    FPR may differ from the nominal target.
    """
    negative_scores = np.asarray(
        negative_scores,
        dtype=float,
    )

    q = 1.0 - float(
        target_fpr
    )

    threshold = float(
        np.quantile(
            negative_scores,
            q,
            method="higher",
        )
    )

    achieved_fpr = float(
        np.mean(
            negative_scores
            >= threshold
        )
    )

    return (
        threshold,
        achieved_fpr,
    )


def carrier_stats(
    d,
    mask_a,
    mask_b,
    group_a_name,
    group_b_name,
):
    """
    Compare MUC5B carrier prevalence between two mutually exclusive
    subgroups within the same phenotype stratum.
    """
    a = d.loc[
        mask_a,
        "MUC5B_T_carrier",
    ].astype(int)

    b = d.loc[
        mask_b,
        "MUC5B_T_carrier",
    ].astype(int)

    if len(a) == 0 or len(b) == 0:
        return {
            "group_A": group_a_name,
            "group_B": group_b_name,
            "N_A": len(a),
            "N_B": len(b),
            "carrier_N_A": int(a.sum()) if len(a) else 0,
            "carrier_N_B": int(b.sum()) if len(b) else 0,
            "carrier_prev_A": float(a.mean()) if len(a) else np.nan,
            "carrier_prev_B": float(b.mean()) if len(b) else np.nan,
            "carrier_prev_ratio_A_over_B": np.nan,
            "OR_A_vs_B": np.nan,
            "OR_ci_low": np.nan,
            "OR_ci_high": np.nan,
            "fisher_p": np.nan,
        }

    a_carrier = int(
        a.sum()
    )
    a_non = int(
        len(a) - a_carrier
    )

    b_carrier = int(
        b.sum()
    )
    b_non = int(
        len(b) - b_carrier
    )

    prev_a = float(
        a.mean()
    )

    prev_b = float(
        b.mean()
    )

    ratio = (
        prev_a / prev_b
        if prev_b > 0
        else np.nan
    )

    OR, OR_lo, OR_hi = odds_ratio_ci(
        a_carrier,
        a_non,
        b_carrier,
        b_non,
    )

    fisher_or, fisher_p = fisher_exact(
        [
            [a_carrier, a_non],
            [b_carrier, b_non],
        ]
    )

    return {
        "group_A":
            group_a_name,

        "group_B":
            group_b_name,

        "N_A":
            len(a),

        "N_B":
            len(b),

        "carrier_N_A":
            a_carrier,

        "carrier_N_B":
            b_carrier,

        "carrier_prev_A":
            prev_a,

        "carrier_prev_B":
            prev_b,

        "carrier_prev_ratio_A_over_B":
            ratio,

        "OR_A_vs_B":
            OR,

        "OR_ci_low":
            OR_lo,

        "OR_ci_high":
            OR_hi,

        "fisher_exact_OR":
            fisher_or,

        "fisher_p":
            fisher_p,
    }


def genotype_distribution(
    d,
    mask,
    confusion_group,
    threshold_name,
):
    sub = d.loc[
        mask
    ]

    n = len(
        sub
    )

    if n == 0:
        return {
            "threshold_name":
                threshold_name,
            "confusion_group":
                confusion_group,
            "N":
                0,
            "GG_N":
                0,
            "GT_N":
                0,
            "TT_N":
                0,
            "GG_fraction":
                np.nan,
            "GT_fraction":
                np.nan,
            "TT_fraction":
                np.nan,
            "T_carrier_fraction":
                np.nan,
            "mean_T_dosage":
                np.nan,
        }

    GG = int(
        sub[
            "MUC5B_GG"
        ].sum()
    )

    GT = int(
        sub[
            "MUC5B_GT"
        ].sum()
    )

    TT = int(
        sub[
            "MUC5B_TT"
        ].sum()
    )

    return {
        "threshold_name":
            threshold_name,

        "confusion_group":
            confusion_group,

        "N":
            n,

        "GG_N":
            GG,

        "GT_N":
            GT,

        "TT_N":
            TT,

        "GG_fraction":
            GG / n,

        "GT_fraction":
            GT / n,

        "TT_fraction":
            TT / n,

        "T_carrier_fraction":
            (GT + TT) / n,

        "mean_T_dosage":
            (GT + 2 * TT) / n,
    }


# ============================================================
# Load source data
# ============================================================

print(
    "Loading source data..."
)

BASE_DATA = pd.read_csv(
    BASE_DATA_FILE
).drop(
    columns=[
        "target"
    ],
    errors="ignore",
)

TARGETS = (
    pd.read_csv(
        TARGET_FILE
    )
    .replace({
        "N": 0,
        "Y": 1,
        "n": 0,
        "y": 1,
        "na": np.nan,
    })
    .rename(
        columns={
            "arb_person_id":
                "patient_id"
        }
    )
)

PREDS = pd.read_parquet(
    PRED_FILE
)


# ============================================================
# Exact notebook-32 cohort construction, while preserving
# whether the FILD/FILA adjudication was actually observed.
# ============================================================

DATA = (
    BASE_DATA
    .merge(
        TARGETS[
            [
                "patient_id",
                TARGET_NAME,
            ]
        ].rename(
            columns={
                TARGET_NAME:
                    "target_original"
            }
        ),
        on="patient_id",
        how="left",
    )
    .merge(
        PREDS[
            [
                "patient_id",
                "predicted_risk",
            ]
        ],
        on="patient_id",
        how="left",
    )
)

DATA[
    "target_was_observed"
] = DATA[
    "target_original"
].notnull()

DATA[
    "target"
] = (
    pd.to_numeric(
        DATA[
            "target_original"
        ],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
)

# EXACT notebook-32 availability filter
DATA = DATA[
    DATA[
        "predicted_risk"
    ].notnull()
].copy()


# ============================================================
# Decode MUC5B genotype.
#
# The source genotype was categorical 0/1/2 and then one-hot
# encoded. Therefore:
#
#   sum == 1 -> observed genotype
#   sum == 0 -> original genotype missing
#   sum > 1  -> invalid encoding
# ============================================================

required = [
    GG_COLUMN,
    GT_COLUMN,
    TT_COLUMN,
]

missing = [
    c
    for c in required
    if c not in DATA.columns
]

if missing:
    raise ValueError(
        "Missing expected MUC5B columns:\n"
        + "\n".join(
            missing
        )
    )

for c in required:
    DATA[
        c
    ] = pd.to_numeric(
        DATA[
            c
        ],
        errors="coerce",
    )

dummy_sum = DATA[
    required
].sum(
    axis=1
)

print()
print(
    "MUC5B one-hot state sums:"
)
print(
    dummy_sum.value_counts()
    .sort_index()
    .to_string()
)

if (
    dummy_sum > 1
).any():
    raise ValueError(
        "Some MUC5B rows have >1 active genotype state."
    )

DATA = DATA.loc[
    dummy_sum == 1
].copy()

DATA[
    "MUC5B_GG"
] = DATA[
    GG_COLUMN
].astype(int)

DATA[
    "MUC5B_GT"
] = DATA[
    GT_COLUMN
].astype(int)

DATA[
    "MUC5B_TT"
] = DATA[
    TT_COLUMN
].astype(int)

DATA[
    "MUC5B_T_carrier"
] = (
    (
        DATA[
            "MUC5B_GT"
        ] == 1
    )
    |
    (
        DATA[
            "MUC5B_TT"
        ] == 1
    )
).astype(int)

DATA[
    "MUC5B_T_dosage"
] = (
    DATA[
        "MUC5B_GT"
    ]
    + 2 * DATA[
        "MUC5B_TT"
    ]
)


# ============================================================
# Basic cohort report
# ============================================================

POS = DATA[
    DATA[
        "target"
    ] == 1
].copy()

NEG = DATA[
    DATA[
        "target"
    ] == 0
].copy()

EXPLICIT_NEG = DATA[
    (
        DATA[
            "target_was_observed"
        ]
    )
    &
    (
        DATA[
            "target"
        ] == 0
    )
].copy()

print()
print("=" * 90)
print("COHORT")
print("=" * 90)

print(
    f"Total with called MUC5B genotype: "
    f"{len(DATA):,}"
)

print(
    f"FILD/FILA positive: "
    f"{len(POS):,}"
)

print(
    f"Notebook-32 negative: "
    f"{len(NEG):,}"
)

print(
    f"Explicitly adjudicated negative: "
    f"{len(EXPLICIT_NEG):,}"
)

print()
print(
    "MUC5B T-carrier prevalence:"
)

print(
    f"  positives: "
    f"{POS['MUC5B_T_carrier'].mean():.3%}"
)

print(
    f"  notebook-32 negatives: "
    f"{NEG['MUC5B_T_carrier'].mean():.3%}"
)

print(
    f"  explicit negatives: "
    f"{EXPLICIT_NEG['MUC5B_T_carrier'].mean():.3%}"
)

print()
print(
    "Raw ZeBRA AUC for FILD/FILA = "
    f"{roc_auc_score(DATA['target'], DATA['predicted_risk']):.4f}"
)


# ============================================================
# Define operating thresholds
# ============================================================

THRESHOLDS = []

negative_scores = NEG[
    "predicted_risk"
].to_numpy(
    dtype=float
)

for target_fpr in TARGET_FPRS:
    threshold, achieved_fpr = (
        threshold_for_target_fpr(
            negative_scores,
            target_fpr,
        )
    )

    sensitivity = float(
        np.mean(
            POS[
                "predicted_risk"
            ].to_numpy(
                dtype=float
            )
            >= threshold
        )
    )

    THRESHOLDS.append({
        "threshold_name":
            f"FPR_{100*target_fpr:g}pct",

        "threshold_source":
            "Notebook-32 target=0",

        "nominal_target_fpr":
            target_fpr,

        "threshold":
            threshold,

        "achieved_fpr_notebook32_neg":
            achieved_fpr,

        "sensitivity_FILD_positive":
            sensitivity,
    })


# Add full-cohort top-1% ZeBRA threshold for continuity with
# the earlier enrichment analysis.
top1_threshold = float(
    np.quantile(
        DATA[
            "predicted_risk"
        ],
        0.99,
        method="higher",
    )
)

THRESHOLDS.append({
    "threshold_name":
        "FULL_COHORT_TOP_1PCT",

    "threshold_source":
        "Full notebook-32 cohort",

    "nominal_target_fpr":
        np.nan,

    "threshold":
        top1_threshold,

    "achieved_fpr_notebook32_neg":
        float(
            np.mean(
                NEG[
                    "predicted_risk"
                ]
                >= top1_threshold
            )
        ),

    "sensitivity_FILD_positive":
        float(
            np.mean(
                POS[
                    "predicted_risk"
                ]
                >= top1_threshold
            )
        ),
})

THRESHOLD_DF = pd.DataFrame(
    THRESHOLDS
)

THRESHOLD_DF.to_csv(
    OUTDIR
    / "ZEBRA_ERROR_ANALYSIS_THRESHOLDS.csv",
    index=False,
)

print()
print("=" * 90)
print("THRESHOLDS")
print("=" * 90)
print(
    THRESHOLD_DF.to_string(
        index=False
    )
)


# ============================================================
# False-negative / false-positive MUC5B analyses
# ============================================================

COMPARISONS = []
GENOTYPE_ROWS = []
CONFUSION_ROWS = []

for row in THRESHOLDS:

    threshold_name = row[
        "threshold_name"
    ]

    threshold = float(
        row[
            "threshold"
        ]
    )

    # --------------------------------------------------------
    # FILD/FILA POSITIVES:
    # TP vs FN
    # --------------------------------------------------------

    pos_pred = (
        POS[
            "predicted_risk"
        ]
        >= threshold
    )

    TP = (
        pos_pred
    )

    FN = (
        ~pos_pred
    )

    n_tp = int(
        TP.sum()
    )

    n_fn = int(
        FN.sum()
    )

    # Primary question:
    # are FALSE NEGATIVES genetically different from TRUE POSITIVES?
    fn_vs_tp = carrier_stats(
        POS,
        FN,
        TP,
        "FN",
        "TP",
    )

    fn_vs_tp.update({
        "threshold_name":
            threshold_name,

        "threshold":
            threshold,

        "phenotype_stratum":
            "FILD/FILA positive",

        "comparison":
            "FN vs TP",
    })

    COMPARISONS.append(
        fn_vs_tp
    )

    GENOTYPE_ROWS.append(
        genotype_distribution(
            POS,
            TP,
            "TP",
            threshold_name,
        )
    )

    GENOTYPE_ROWS[-1][
        "phenotype_stratum"
    ] = "FILD/FILA positive"

    GENOTYPE_ROWS.append(
        genotype_distribution(
            POS,
            FN,
            "FN",
            threshold_name,
        )
    )

    GENOTYPE_ROWS[-1][
        "phenotype_stratum"
    ] = "FILD/FILA positive"

    # --------------------------------------------------------
    # NOTEBOOK-32 NEGATIVES:
    # FP vs TN
    # --------------------------------------------------------

    neg_pred = (
        NEG[
            "predicted_risk"
        ]
        >= threshold
    )

    FP = (
        neg_pred
    )

    TN = (
        ~neg_pred
    )

    n_fp = int(
        FP.sum()
    )

    n_tn = int(
        TN.sum()
    )

    # Primary question:
    # are FALSE POSITIVES enriched for MUC5B relative to TRUE NEGATIVES?
    fp_vs_tn = carrier_stats(
        NEG,
        FP,
        TN,
        "FP",
        "TN",
    )

    fp_vs_tn.update({
        "threshold_name":
            threshold_name,

        "threshold":
            threshold,

        "phenotype_stratum":
            "Notebook-32 FILD/FILA negative",

        "comparison":
            "FP vs TN",
    })

    COMPARISONS.append(
        fp_vs_tn
    )

    GENOTYPE_ROWS.append(
        genotype_distribution(
            NEG,
            FP,
            "FP",
            threshold_name,
        )
    )

    GENOTYPE_ROWS[-1][
        "phenotype_stratum"
    ] = "Notebook-32 FILD/FILA negative"

    GENOTYPE_ROWS.append(
        genotype_distribution(
            NEG,
            TN,
            "TN",
            threshold_name,
        )
    )

    GENOTYPE_ROWS[-1][
        "phenotype_stratum"
    ] = "Notebook-32 FILD/FILA negative"

    # --------------------------------------------------------
    # EXPLICITLY ADJUDICATED NEGATIVES:
    # FP vs TN using the SAME threshold
    # --------------------------------------------------------

    explicit_pred = (
        EXPLICIT_NEG[
            "predicted_risk"
        ]
        >= threshold
    )

    EFP = explicit_pred
    ETN = ~explicit_pred

    explicit_fp_vs_tn = carrier_stats(
        EXPLICIT_NEG,
        EFP,
        ETN,
        "FP",
        "TN",
    )

    explicit_fp_vs_tn.update({
        "threshold_name":
            threshold_name,

        "threshold":
            threshold,

        "phenotype_stratum":
            "Explicitly adjudicated FILD/FILA negative",

        "comparison":
            "FP vs TN",
    })

    COMPARISONS.append(
        explicit_fp_vs_tn
    )

    GENOTYPE_ROWS.append(
        genotype_distribution(
            EXPLICIT_NEG,
            EFP,
            "FP",
            threshold_name,
        )
    )

    GENOTYPE_ROWS[-1][
        "phenotype_stratum"
    ] = "Explicitly adjudicated FILD/FILA negative"

    GENOTYPE_ROWS.append(
        genotype_distribution(
            EXPLICIT_NEG,
            ETN,
            "TN",
            threshold_name,
        )
    )

    GENOTYPE_ROWS[-1][
        "phenotype_stratum"
    ] = "Explicitly adjudicated FILD/FILA negative"

    # --------------------------------------------------------
    # Confusion counts
    # --------------------------------------------------------

    CONFUSION_ROWS.append({
        "threshold_name":
            threshold_name,

        "threshold":
            threshold,

        "TP":
            n_tp,

        "FN":
            n_fn,

        "FP":
            n_fp,

        "TN":
            n_tn,

        "sensitivity":
            n_tp / (
                n_tp + n_fn
            ),

        "FPR":
            n_fp / (
                n_fp + n_tn
            ),

        "PPV":
            n_tp / (
                n_tp + n_fp
            )
            if (
                n_tp + n_fp
            ) > 0
            else np.nan,
    })


COMPARISON_DF = pd.DataFrame(
    COMPARISONS
)

GENOTYPE_DF = pd.DataFrame(
    GENOTYPE_ROWS
)

CONFUSION_DF = pd.DataFrame(
    CONFUSION_ROWS
)

COMPARISON_DF.to_csv(
    OUTDIR
    / "MUC5B_FALSE_NEGATIVE_FALSE_POSITIVE_COMPARISONS.csv",
    index=False,
)

GENOTYPE_DF.to_csv(
    OUTDIR
    / "MUC5B_GENOTYPE_BY_CONFUSION_GROUP.csv",
    index=False,
)

CONFUSION_DF.to_csv(
    OUTDIR
    / "ZEBRA_CONFUSION_COUNTS_BY_THRESHOLD.csv",
    index=False,
)


# ============================================================
# Console summaries
# ============================================================

print()
print("=" * 90)
print("FALSE NEGATIVES AMONG FILD/FILA POSITIVES")
print("=" * 90)

fn_table = COMPARISON_DF[
    COMPARISON_DF[
        "comparison"
    ] == "FN vs TP"
][
    [
        "threshold_name",
        "N_A",
        "N_B",
        "carrier_prev_A",
        "carrier_prev_B",
        "carrier_prev_ratio_A_over_B",
        "OR_A_vs_B",
        "OR_ci_low",
        "OR_ci_high",
        "fisher_p",
    ]
].copy()

fn_table = fn_table.rename(
    columns={
        "N_A":
            "FN_N",
        "N_B":
            "TP_N",
        "carrier_prev_A":
            "FN_carrier_prev",
        "carrier_prev_B":
            "TP_carrier_prev",
        "carrier_prev_ratio_A_over_B":
            "FN_over_TP_prevalence_ratio",
        "OR_A_vs_B":
            "FN_vs_TP_OR",
    }
)

print(
    fn_table.to_string(
        index=False
    )
)


print()
print("=" * 90)
print("FALSE POSITIVES AMONG NOTEBOOK-32 FILD/FILA NEGATIVES")
print("=" * 90)

fp_table = COMPARISON_DF[
    (
        COMPARISON_DF[
            "comparison"
        ] == "FP vs TN"
    )
    &
    (
        COMPARISON_DF[
            "phenotype_stratum"
        ]
        == "Notebook-32 FILD/FILA negative"
    )
][
    [
        "threshold_name",
        "N_A",
        "N_B",
        "carrier_prev_A",
        "carrier_prev_B",
        "carrier_prev_ratio_A_over_B",
        "OR_A_vs_B",
        "OR_ci_low",
        "OR_ci_high",
        "fisher_p",
    ]
].copy()

fp_table = fp_table.rename(
    columns={
        "N_A":
            "FP_N",
        "N_B":
            "TN_N",
        "carrier_prev_A":
            "FP_carrier_prev",
        "carrier_prev_B":
            "TN_carrier_prev",
        "carrier_prev_ratio_A_over_B":
            "FP_over_TN_prevalence_ratio",
        "OR_A_vs_B":
            "FP_vs_TN_OR",
    }
)

print(
    fp_table.to_string(
        index=False
    )
)


print()
print("=" * 90)
print("FALSE POSITIVES AMONG EXPLICITLY ADJUDICATED NEGATIVES")
print("=" * 90)

efp_table = COMPARISON_DF[
    (
        COMPARISON_DF[
            "comparison"
        ] == "FP vs TN"
    )
    &
    (
        COMPARISON_DF[
            "phenotype_stratum"
        ]
        == "Explicitly adjudicated FILD/FILA negative"
    )
][
    [
        "threshold_name",
        "N_A",
        "N_B",
        "carrier_prev_A",
        "carrier_prev_B",
        "carrier_prev_ratio_A_over_B",
        "OR_A_vs_B",
        "OR_ci_low",
        "OR_ci_high",
        "fisher_p",
    ]
].copy()

efp_table = efp_table.rename(
    columns={
        "N_A":
            "FP_N",
        "N_B":
            "TN_N",
        "carrier_prev_A":
            "FP_carrier_prev",
        "carrier_prev_B":
            "TN_carrier_prev",
        "carrier_prev_ratio_A_over_B":
            "FP_over_TN_prevalence_ratio",
        "OR_A_vs_B":
            "FP_vs_TN_OR",
    }
)

print(
    efp_table.to_string(
        index=False
    )
)


# ============================================================
# Plots
# ============================================================

# ------------------------------------------------------------
# 1. MUC5B carrier prevalence among FN vs TP
# ------------------------------------------------------------

plot_fn = COMPARISON_DF[
    COMPARISON_DF[
        "comparison"
    ] == "FN vs TP"
].copy()

x = np.arange(
    len(
        plot_fn
    )
)

fig, ax = plt.subplots(
    figsize=(10, 5)
)

width = 0.38

ax.bar(
    x - width / 2,
    plot_fn[
        "carrier_prev_A"
    ],
    width=width,
    label="False negatives",
)

ax.bar(
    x + width / 2,
    plot_fn[
        "carrier_prev_B"
    ],
    width=width,
    label="True positives",
)

ax.set_xticks(
    x
)

ax.set_xticklabels(
    plot_fn[
        "threshold_name"
    ],
    rotation=30,
    ha="right",
)

ax.set_ylabel(
    "MUC5B T-carrier prevalence"
)

ax.set_title(
    "FILD/FILA positives: MUC5B carriage in false negatives vs true positives"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MUC5B_FN_vs_TP_carrier_prevalence.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ------------------------------------------------------------
# 2. MUC5B carrier prevalence among FP vs TN
#    notebook-32 negatives
# ------------------------------------------------------------

plot_fp = COMPARISON_DF[
    (
        COMPARISON_DF[
            "comparison"
        ] == "FP vs TN"
    )
    &
    (
        COMPARISON_DF[
            "phenotype_stratum"
        ]
        == "Notebook-32 FILD/FILA negative"
    )
].copy()

x = np.arange(
    len(
        plot_fp
    )
)

fig, ax = plt.subplots(
    figsize=(10, 5)
)

ax.bar(
    x - width / 2,
    plot_fp[
        "carrier_prev_A"
    ],
    width=width,
    label="False positives",
)

ax.bar(
    x + width / 2,
    plot_fp[
        "carrier_prev_B"
    ],
    width=width,
    label="True negatives",
)

ax.set_xticks(
    x
)

ax.set_xticklabels(
    plot_fp[
        "threshold_name"
    ],
    rotation=30,
    ha="right",
)

ax.set_ylabel(
    "MUC5B T-carrier prevalence"
)

ax.set_title(
    "Notebook-32 FILD/FILA negatives: MUC5B carriage in false positives vs true negatives"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MUC5B_FP_vs_TN_carrier_prevalence.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ------------------------------------------------------------
# 3. Enrichment ratios:
#    FN/TP and FP/TN
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(10, 5)
)

ax.plot(
    plot_fn[
        "threshold_name"
    ],
    plot_fn[
        "carrier_prev_ratio_A_over_B"
    ],
    marker="o",
    label="FN / TP carrier prevalence",
)

ax.plot(
    plot_fp[
        "threshold_name"
    ],
    plot_fp[
        "carrier_prev_ratio_A_over_B"
    ],
    marker="o",
    label="FP / TN carrier prevalence",
)

ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_ylabel(
    "MUC5B carrier prevalence ratio"
)

ax.set_title(
    "MUC5B enrichment in ZeBRA errors"
)

ax.tick_params(
    axis="x",
    labelrotation=30,
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MUC5B_error_enrichment_ratios.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Final note
# ============================================================

print()
print("=" * 90)
print("INTERPRETATION GUIDE")
print("=" * 90)

print(
    "If FN carrier prevalence < TP carrier prevalence:"
)
print(
    "  ZeBRA is preferentially detecting genetically driven FILD cases."
)

print()
print(
    "If FN carrier prevalence > TP carrier prevalence:"
)
print(
    "  ZeBRA may be missing a genetically enriched subtype."
)

print()
print(
    "If FP carrier prevalence > TN carrier prevalence:"
)
print(
    "  ZeBRA false positives may carry latent inherited fibrotic susceptibility."
)

print()
print(
    "If FP carrier prevalence ~= TN carrier prevalence:"
)
print(
    "  There is little evidence that ZeBRA false positives are genetically enriched."
)

print()
print(
    "Outputs saved to:"
)
print(
    OUTDIR.resolve()
)
