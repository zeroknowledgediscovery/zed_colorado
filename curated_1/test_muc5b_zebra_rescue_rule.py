#!/usr/bin/env python3

from pathlib import Path
import secrets

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
from sklearn.model_selection import train_test_split
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

# Repeated outer split, matching the logic used in notebook 32.
TRAIN_SIZE = 0.40
N_REPEATS = 100

# Main clinically stringent ZeBRA operating points.
MAIN_FPRS = [
    0.005,
    0.010,
]

# Candidate lower ZeBRA thresholds defining the rescue band.
# Must be less stringent (higher FPR) than MAIN_FPRS.
LOWER_FPRS = [
    0.020,
    0.050,
    0.100,
]

OUTDIR = Path(
    "./RESULTS/MUC5B_ZEBRA_RESCUE_RULE"
)
OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Helpers
# ============================================================

def threshold_for_target_fpr(
    negative_scores,
    target_fpr,
):
    """
    Derive a ZeBRA threshold from TRAINING negatives only.

    Uses the upper-tail quantile corresponding approximately to
    the requested FPR. Because ZeBRA can have ties near 1.0,
    achieved FPR may differ from nominal FPR.
    """
    negative_scores = np.asarray(
        negative_scores,
        dtype=float,
    )

    threshold = float(
        np.quantile(
            negative_scores,
            1.0 - target_fpr,
            method="higher",
        )
    )

    achieved = float(
        np.mean(
            negative_scores >= threshold
        )
    )

    return threshold, achieved


def confusion_metrics(
    y,
    positive,
):
    y = np.asarray(
        y,
        dtype=int,
    )
    positive = np.asarray(
        positive,
        dtype=bool,
    )

    tp = int(
        np.sum(
            (y == 1) & positive
        )
    )

    fn = int(
        np.sum(
            (y == 1) & (~positive)
        )
    )

    fp = int(
        np.sum(
            (y == 0) & positive
        )
    )

    tn = int(
        np.sum(
            (y == 0) & (~positive)
        )
    )

    sensitivity = (
        tp / (tp + fn)
        if (tp + fn) > 0
        else np.nan
    )

    specificity = (
        tn / (tn + fp)
        if (tn + fp) > 0
        else np.nan
    )

    fpr = (
        fp / (fp + tn)
        if (fp + tn) > 0
        else np.nan
    )

    ppv = (
        tp / (tp + fp)
        if (tp + fp) > 0
        else np.nan
    )

    npv = (
        tn / (tn + fn)
        if (tn + fn) > 0
        else np.nan
    )

    return {
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "FPR": fpr,
        "PPV": ppv,
        "NPV": npv,
    }


def band_genetic_association(
    y,
    carrier,
    in_band,
):
    """
    Within the intermediate ZeBRA band only:
    test whether MUC5B carrier status is associated with FILD/FILA.
    """
    yb = np.asarray(
        y,
        dtype=int,
    )[in_band]

    cb = np.asarray(
        carrier,
        dtype=int,
    )[in_band]

    n_band = len(
        yb
    )

    if (
        n_band == 0
        or np.unique(yb).size < 2
        or np.unique(cb).size < 2
    ):
        return {
            "band_N": n_band,
            "band_cases": int(yb.sum()) if n_band else 0,
            "band_carriers": int(cb.sum()) if n_band else 0,
            "case_prev_carrier": np.nan,
            "case_prev_noncarrier": np.nan,
            "carrier_case_OR": np.nan,
            "carrier_case_fisher_p": np.nan,
        }

    case_carrier = int(
        np.sum(
            (yb == 1)
            & (cb == 1)
        )
    )

    control_carrier = int(
        np.sum(
            (yb == 0)
            & (cb == 1)
        )
    )

    case_noncarrier = int(
        np.sum(
            (yb == 1)
            & (cb == 0)
        )
    )

    control_noncarrier = int(
        np.sum(
            (yb == 0)
            & (cb == 0)
        )
    )

    carrier_case_prev = (
        case_carrier
        / (
            case_carrier
            + control_carrier
        )
        if (
            case_carrier
            + control_carrier
        ) > 0
        else np.nan
    )

    noncarrier_case_prev = (
        case_noncarrier
        / (
            case_noncarrier
            + control_noncarrier
        )
        if (
            case_noncarrier
            + control_noncarrier
        ) > 0
        else np.nan
    )

    OR, p = fisher_exact(
        [
            [
                case_carrier,
                control_carrier,
            ],
            [
                case_noncarrier,
                control_noncarrier,
            ],
        ]
    )

    return {
        "band_N":
            n_band,

        "band_cases":
            int(
                yb.sum()
            ),

        "band_carriers":
            int(
                cb.sum()
            ),

        "case_prev_carrier":
            carrier_case_prev,

        "case_prev_noncarrier":
            noncarrier_case_prev,

        "carrier_case_OR":
            float(OR),

        "carrier_case_fisher_p":
            float(p),
    }


def summarize_repeated(
    df,
    value_columns,
    group_columns,
):
    rows = []

    for keys, sub in df.groupby(
        group_columns
    ):
        if not isinstance(
            keys,
            tuple,
        ):
            keys = (
                keys,
            )

        row = dict(
            zip(
                group_columns,
                keys,
            )
        )

        row[
            "n_repeats"
        ] = len(
            sub
        )

        for col in value_columns:
            vals = (
                pd.to_numeric(
                    sub[col],
                    errors="coerce",
                )
                .dropna()
                .to_numpy(
                    dtype=float
                )
            )

            if len(
                vals
            ) == 0:
                continue

            row[
                f"{col}_mean"
            ] = float(
                np.mean(
                    vals
                )
            )

            row[
                f"{col}_median"
            ] = float(
                np.median(
                    vals
                )
            )

            row[
                f"{col}_q025"
            ] = float(
                np.quantile(
                    vals,
                    0.025,
                )
            )

            row[
                f"{col}_q975"
            ] = float(
                np.quantile(
                    vals,
                    0.975,
                )
            )

            row[
                f"{col}_fraction_gt_zero"
            ] = float(
                np.mean(
                    vals > 0
                )
            )

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# Load data
# ============================================================

print(
    "Loading data..."
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
# Exact notebook-32 cohort construction
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

# Exact notebook-32 ZeBRA availability filter.
DATA = DATA[
    DATA[
        "predicted_risk"
    ].notnull()
].copy()


# ============================================================
# Decode MUC5B rs35705950 genotype
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
    "MUC5B dummy-state sums:"
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
        "Some MUC5B rows are multi-hot."
    )

# Restrict to observed/called MUC5B genotype.
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
    + 2
    * DATA[
        "MUC5B_TT"
    ]
)


# ============================================================
# Cohort report
# ============================================================

print()
print("=" * 96)
print("COHORT")
print("=" * 96)

print(
    f"N with ZeBRA + called MUC5B genotype = "
    f"{len(DATA):,}"
)

print(
    f"FILD/FILA positives = "
    f"{int(DATA['target'].sum()):,}"
)

print(
    f"Notebook-32 negatives = "
    f"{int((DATA['target'] == 0).sum()):,}"
)

print(
    f"Raw ZeBRA FILD/FILA AUC = "
    f"{roc_auc_score(DATA['target'], DATA['predicted_risk']):.4f}"
)

print(
    f"MUC5B T-carrier prevalence in positives = "
    f"{DATA.loc[DATA['target'] == 1, 'MUC5B_T_carrier'].mean():.3%}"
)

print(
    f"MUC5B T-carrier prevalence in negatives = "
    f"{DATA.loc[DATA['target'] == 0, 'MUC5B_T_carrier'].mean():.3%}"
)


# ============================================================
# Full-cohort descriptive rescue analysis
# ============================================================

print()
print("=" * 96)
print("FULL-COHORT DESCRIPTIVE RESCUE ANALYSIS")
print("=" * 96)

FULL_ROWS = []

y_full = DATA[
    "target"
].to_numpy(
    dtype=int
)

z_full = DATA[
    "predicted_risk"
].to_numpy(
    dtype=float
)

carrier_full = DATA[
    "MUC5B_T_carrier"
].to_numpy(
    dtype=int
)

neg_scores_full = z_full[
    y_full == 0
]

for main_fpr in MAIN_FPRS:

    t_main, achieved_main = (
        threshold_for_target_fpr(
            neg_scores_full,
            main_fpr,
        )
    )

    baseline_positive = (
        z_full >= t_main
    )

    baseline = confusion_metrics(
        y_full,
        baseline_positive,
    )

    for lower_fpr in LOWER_FPRS:

        if lower_fpr <= main_fpr:
            continue

        t_lower, achieved_lower = (
            threshold_for_target_fpr(
                neg_scores_full,
                lower_fpr,
            )
        )

        rescue_band = (
            (z_full >= t_lower)
            &
            (z_full < t_main)
        )

        rescue_selected = (
            rescue_band
            &
            (carrier_full == 1)
        )

        rescued_positive = (
            baseline_positive
            |
            rescue_selected
        )

        rescued = confusion_metrics(
            y_full,
            rescued_positive,
        )

        rescued_tp = (
            rescued[
                "TP"
            ]
            - baseline[
                "TP"
            ]
        )

        added_fp = (
            rescued[
                "FP"
            ]
            - baseline[
                "FP"
            ]
        )

        baseline_fn = baseline[
            "FN"
        ]

        fraction_fn_rescued = (
            rescued_tp
            / baseline_fn
            if baseline_fn > 0
            else np.nan
        )

        rescue_ppv = (
            rescued_tp
            / (
                rescued_tp
                + added_fp
            )
            if (
                rescued_tp
                + added_fp
            ) > 0
            else np.nan
        )

        fp_per_rescued_tp = (
            added_fp
            / rescued_tp
            if rescued_tp > 0
            else np.inf
        )

        band_assoc = (
            band_genetic_association(
                y_full,
                carrier_full,
                rescue_band,
            )
        )

        row = {
            "main_target_fpr":
                main_fpr,

            "lower_target_fpr":
                lower_fpr,

            "main_threshold":
                t_main,

            "lower_threshold":
                t_lower,

            "baseline_sensitivity":
                baseline[
                    "sensitivity"
                ],

            "baseline_FPR":
                baseline[
                    "FPR"
                ],

            "baseline_PPV":
                baseline[
                    "PPV"
                ],

            "rescued_sensitivity":
                rescued[
                    "sensitivity"
                ],

            "rescued_FPR":
                rescued[
                    "FPR"
                ],

            "rescued_PPV":
                rescued[
                    "PPV"
                ],

            "delta_sensitivity":
                (
                    rescued[
                        "sensitivity"
                    ]
                    - baseline[
                        "sensitivity"
                    ]
                ),

            "delta_FPR":
                (
                    rescued[
                        "FPR"
                    ]
                    - baseline[
                        "FPR"
                    ]
                ),

            "delta_PPV":
                (
                    rescued[
                        "PPV"
                    ]
                    - baseline[
                        "PPV"
                    ]
                ),

            "rescued_TP":
                rescued_tp,

            "added_FP":
                added_fp,

            "fraction_baseline_FN_rescued":
                fraction_fn_rescued,

            "rescue_subgroup_PPV":
                rescue_ppv,

            "FP_per_rescued_TP":
                fp_per_rescued_tp,

            **band_assoc,
        }

        FULL_ROWS.append(
            row
        )


FULL_DF = pd.DataFrame(
    FULL_ROWS
)

FULL_DF.to_csv(
    OUTDIR
    / "FULL_COHORT_MUC5B_RESCUE_RULE.csv",
    index=False,
)

print(
    FULL_DF[
        [
            "main_target_fpr",
            "lower_target_fpr",
            "baseline_sensitivity",
            "baseline_FPR",
            "rescued_sensitivity",
            "rescued_FPR",
            "delta_sensitivity",
            "delta_FPR",
            "rescued_TP",
            "added_FP",
            "fraction_baseline_FN_rescued",
            "rescue_subgroup_PPV",
            "FP_per_rescued_TP",
            "carrier_case_OR",
            "carrier_case_fisher_p",
        ]
    ].to_string(
        index=False
    )
)


# ============================================================
# Repeated split validation
#
# Thresholds are estimated ONLY from training negatives.
# Evaluation is ONLY on untouched outer-test patients.
# ============================================================

print()
print("=" * 96)
print(
    f"REPEATED SPLIT VALIDATION: {N_REPEATS} repeats"
)
print("=" * 96)

REPEAT_ROWS = []

indices = np.arange(
    len(DATA)
)

y_all = DATA[
    "target"
].to_numpy(
    dtype=int
)

z_all = DATA[
    "predicted_risk"
].to_numpy(
    dtype=float
)

carrier_all = DATA[
    "MUC5B_T_carrier"
].to_numpy(
    dtype=int
)

observed_target_all = DATA[
    "target_was_observed"
].to_numpy(
    dtype=bool
)

for repeat in range(
    1,
    N_REPEATS + 1,
):

    split_seed = secrets.randbits(
        31
    )

    train_idx, test_idx = train_test_split(
        indices,
        train_size=TRAIN_SIZE,
        stratify=y_all,
        random_state=split_seed,
    )

    y_train = y_all[
        train_idx
    ]

    y_test = y_all[
        test_idx
    ]

    z_train = z_all[
        train_idx
    ]

    z_test = z_all[
        test_idx
    ]

    c_test = carrier_all[
        test_idx
    ]

    observed_test = observed_target_all[
        test_idx
    ]

    train_negative_scores = (
        z_train[
            y_train == 0
        ]
    )

    for main_fpr in MAIN_FPRS:

        t_main, train_main_fpr = (
            threshold_for_target_fpr(
                train_negative_scores,
                main_fpr,
            )
        )

        baseline_positive = (
            z_test >= t_main
        )

        baseline = confusion_metrics(
            y_test,
            baseline_positive,
        )

        for lower_fpr in LOWER_FPRS:

            if lower_fpr <= main_fpr:
                continue

            t_lower, train_lower_fpr = (
                threshold_for_target_fpr(
                    train_negative_scores,
                    lower_fpr,
                )
            )

            rescue_band = (
                (z_test >= t_lower)
                &
                (z_test < t_main)
            )

            rescue_selected = (
                rescue_band
                &
                (c_test == 1)
            )

            rescued_positive = (
                baseline_positive
                |
                rescue_selected
            )

            rescued = confusion_metrics(
                y_test,
                rescued_positive,
            )

            rescued_tp = (
                rescued[
                    "TP"
                ]
                - baseline[
                    "TP"
                ]
            )

            added_fp = (
                rescued[
                    "FP"
                ]
                - baseline[
                    "FP"
                ]
            )

            baseline_fn = baseline[
                "FN"
            ]

            fraction_fn_rescued = (
                rescued_tp
                / baseline_fn
                if baseline_fn > 0
                else np.nan
            )

            rescue_ppv = (
                rescued_tp
                / (
                    rescued_tp
                    + added_fp
                )
                if (
                    rescued_tp
                    + added_fp
                ) > 0
                else np.nan
            )

            fp_per_rescued_tp = (
                added_fp
                / rescued_tp
                if rescued_tp > 0
                else np.nan
            )

            # ------------------------------------------------
            # Explicitly adjudicated negative FPR increment,
            # using the same test-set rule.
            # ------------------------------------------------

            explicit_negative = (
                (y_test == 0)
                &
                observed_test
            )

            if explicit_negative.sum() > 0:

                base_explicit_fp = int(
                    np.sum(
                        baseline_positive
                        &
                        explicit_negative
                    )
                )

                rescue_explicit_fp = int(
                    np.sum(
                        rescued_positive
                        &
                        explicit_negative
                    )
                )

                added_explicit_fp = (
                    rescue_explicit_fp
                    - base_explicit_fp
                )

                explicit_negative_N = int(
                    explicit_negative.sum()
                )

                baseline_explicit_fpr = (
                    base_explicit_fp
                    / explicit_negative_N
                )

                rescued_explicit_fpr = (
                    rescue_explicit_fp
                    / explicit_negative_N
                )

            else:
                added_explicit_fp = np.nan
                baseline_explicit_fpr = np.nan
                rescued_explicit_fpr = np.nan

            band_assoc = (
                band_genetic_association(
                    y_test,
                    c_test,
                    rescue_band,
                )
            )

            REPEAT_ROWS.append({
                "repeat":
                    repeat,

                "split_seed":
                    split_seed,

                "main_target_fpr":
                    main_fpr,

                "lower_target_fpr":
                    lower_fpr,

                "main_threshold":
                    t_main,

                "lower_threshold":
                    t_lower,

                "train_achieved_main_fpr":
                    train_main_fpr,

                "train_achieved_lower_fpr":
                    train_lower_fpr,

                "baseline_sensitivity":
                    baseline[
                        "sensitivity"
                    ],

                "baseline_FPR":
                    baseline[
                        "FPR"
                    ],

                "baseline_PPV":
                    baseline[
                        "PPV"
                    ],

                "rescued_sensitivity":
                    rescued[
                        "sensitivity"
                    ],

                "rescued_FPR":
                    rescued[
                        "FPR"
                    ],

                "rescued_PPV":
                    rescued[
                        "PPV"
                    ],

                "delta_sensitivity":
                    (
                        rescued[
                            "sensitivity"
                        ]
                        - baseline[
                            "sensitivity"
                        ]
                    ),

                "delta_FPR":
                    (
                        rescued[
                            "FPR"
                        ]
                        - baseline[
                            "FPR"
                        ]
                    ),

                "delta_PPV":
                    (
                        rescued[
                            "PPV"
                        ]
                        - baseline[
                            "PPV"
                        ]
                    ),

                "rescued_TP":
                    rescued_tp,

                "added_FP":
                    added_fp,

                "fraction_baseline_FN_rescued":
                    fraction_fn_rescued,

                "rescue_subgroup_PPV":
                    rescue_ppv,

                "FP_per_rescued_TP":
                    fp_per_rescued_tp,

                "baseline_explicit_negative_FPR":
                    baseline_explicit_fpr,

                "rescued_explicit_negative_FPR":
                    rescued_explicit_fpr,

                "delta_explicit_negative_FPR":
                    (
                        rescued_explicit_fpr
                        - baseline_explicit_fpr
                        if (
                            np.isfinite(
                                rescued_explicit_fpr
                            )
                            and np.isfinite(
                                baseline_explicit_fpr
                            )
                        )
                        else np.nan
                    ),

                "added_explicit_negative_FP":
                    added_explicit_fp,

                **band_assoc,
            })

    if (
        repeat % 10
    ) == 0:
        print(
            f"completed repeat "
            f"{repeat}/{N_REPEATS}"
        )


REPEAT_DF = pd.DataFrame(
    REPEAT_ROWS
)

REPEAT_DF.to_csv(
    OUTDIR
    / "REPEATED_SPLIT_MUC5B_RESCUE_RULE.csv",
    index=False,
)


# ============================================================
# Aggregate repeated split results
# ============================================================

SUMMARY_COLS = [
    "baseline_sensitivity",
    "baseline_FPR",
    "baseline_PPV",

    "rescued_sensitivity",
    "rescued_FPR",
    "rescued_PPV",

    "delta_sensitivity",
    "delta_FPR",
    "delta_PPV",

    "rescued_TP",
    "added_FP",

    "fraction_baseline_FN_rescued",
    "rescue_subgroup_PPV",
    "FP_per_rescued_TP",

    "baseline_explicit_negative_FPR",
    "rescued_explicit_negative_FPR",
    "delta_explicit_negative_FPR",

    "carrier_case_OR",
]

SUMMARY = summarize_repeated(
    REPEAT_DF,
    SUMMARY_COLS,
    [
        "main_target_fpr",
        "lower_target_fpr",
    ],
)

SUMMARY.to_csv(
    OUTDIR
    / "REPEATED_SPLIT_MUC5B_RESCUE_SUMMARY.csv",
    index=False,
)

print()
print("=" * 96)
print("REPEATED-SPLIT SUMMARY")
print("=" * 96)

display_cols = [
    "main_target_fpr",
    "lower_target_fpr",
    "n_repeats",

    "baseline_sensitivity_mean",
    "rescued_sensitivity_mean",
    "delta_sensitivity_mean",
    "delta_sensitivity_q025",
    "delta_sensitivity_q975",

    "baseline_FPR_mean",
    "rescued_FPR_mean",
    "delta_FPR_mean",

    "fraction_baseline_FN_rescued_mean",
    "rescue_subgroup_PPV_mean",
    "FP_per_rescued_TP_mean",

    "delta_explicit_negative_FPR_mean",

    "carrier_case_OR_mean",
]

available_display_cols = [
    c
    for c in display_cols
    if c in SUMMARY.columns
]

print(
    SUMMARY[
        available_display_cols
    ].to_string(
        index=False
    )
)


# ============================================================
# Rank rescue configurations
#
# A simple descriptive efficiency metric:
#
# sensitivity gained per 1 percentage-point FPR added
# ============================================================

SUMMARY[
    "sensitivity_gain_per_1pct_FPR"
] = (
    SUMMARY[
        "delta_sensitivity_mean"
    ]
    / SUMMARY[
        "delta_FPR_mean"
    ]
    * 0.01
)

SUMMARY[
    "rescue_efficiency_TP_per_FP"
] = (
    1.0
    / SUMMARY[
        "FP_per_rescued_TP_mean"
    ]
)

SUMMARY = SUMMARY.sort_values(
    [
        "sensitivity_gain_per_1pct_FPR",
        "rescue_subgroup_PPV_mean",
    ],
    ascending=False,
)

SUMMARY.to_csv(
    OUTDIR
    / "REPEATED_SPLIT_MUC5B_RESCUE_SUMMARY_RANKED.csv",
    index=False,
)

print()
print("=" * 96)
print("RANKED RESCUE CONFIGURATIONS")
print("=" * 96)

rank_cols = [
    "main_target_fpr",
    "lower_target_fpr",

    "delta_sensitivity_mean",
    "delta_FPR_mean",

    "fraction_baseline_FN_rescued_mean",
    "rescue_subgroup_PPV_mean",
    "FP_per_rescued_TP_mean",

    "sensitivity_gain_per_1pct_FPR",
]

print(
    SUMMARY[
        rank_cols
    ].to_string(
        index=False
    )
)


# ============================================================
# Plots
# ============================================================

# ------------------------------------------------------------
# 1. Incremental sensitivity vs incremental FPR
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(8, 6)
)

for _, row in SUMMARY.iterrows():

    ax.scatter(
        row[
            "delta_FPR_mean"
        ],
        row[
            "delta_sensitivity_mean"
        ],
        s=80,
    )

    ax.annotate(
        (
            f"main {100*row['main_target_fpr']:.1f}%"
            f" → lower {100*row['lower_target_fpr']:.1f}%"
        ),
        (
            row[
                "delta_FPR_mean"
            ],
            row[
                "delta_sensitivity_mean"
            ],
        ),
        xytext=(
            5,
            5,
        ),
        textcoords="offset points",
        fontsize=9,
    )

ax.set_xlabel(
    "Incremental FPR from MUC5B rescue"
)

ax.set_ylabel(
    "Incremental sensitivity from MUC5B rescue"
)

ax.set_title(
    "MUC5B rescue rule: sensitivity gained vs false positives added"
)

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MUC5B_RESCUE_delta_sensitivity_vs_delta_FPR.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ------------------------------------------------------------
# 2. Fraction of ZeBRA false negatives rescued
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(10, 5)
)

labels = [
    (
        f"{100*r.main_target_fpr:.1f}%→"
        f"{100*r.lower_target_fpr:.1f}% FPR"
    )
    for r in SUMMARY.itertuples()
]

x = np.arange(
    len(
        SUMMARY
    )
)

ax.bar(
    x,
    SUMMARY[
        "fraction_baseline_FN_rescued_mean"
    ],
)

ax.set_xticks(
    x
)

ax.set_xticklabels(
    labels,
    rotation=30,
    ha="right",
)

ax.set_ylabel(
    "Fraction of baseline ZeBRA false negatives rescued"
)

ax.set_title(
    "MUC5B-based rescue of otherwise missed FILD/FILA cases"
)

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MUC5B_RESCUE_fraction_false_negatives_rescued.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ------------------------------------------------------------
# 3. Added false positives per rescued true positive
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(10, 5)
)

ax.bar(
    x,
    SUMMARY[
        "FP_per_rescued_TP_mean"
    ],
)

ax.set_xticks(
    x
)

ax.set_xticklabels(
    labels,
    rotation=30,
    ha="right",
)

ax.set_ylabel(
    "Added false positives per rescued true positive"
)

ax.set_title(
    "Cost of MUC5B rescue"
)

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MUC5B_RESCUE_FP_per_rescued_TP.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Interpretation guide
# ============================================================

print()
print("=" * 96)
print("INTERPRETATION GUIDE")
print("=" * 96)

print(
    "The rescue rule is useful if, on untouched test sets:"
)
print(
    "  1. delta_sensitivity is consistently > 0,"
)
print(
    "  2. delta_FPR is small,"
)
print(
    "  3. rescue_subgroup_PPV is materially higher than the "
    "baseline prevalence among the intermediate-score band,"
)
print(
    "  4. FP_per_rescued_TP is acceptable for the intended screening use,"
)
print(
    "  5. the intermediate band shows carrier_case_OR > 1."
)

print()
print(
    "A clinically useful result would look like:"
)
print(
    "  'At a ZeBRA threshold corresponding to 1% FPR, adding MUC5B "
    "carrier status only for patients in the 1%-to-5% ZeBRA risk band "
    "recovers X% of otherwise missed cases while adding Y% absolute FPR.'"
)

print()
print(
    "Outputs saved to:"
)
print(
    OUTDIR.resolve()
)
