#!/usr/bin/env python3

from pathlib import Path
import secrets

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split


# ============================================================
# CONFIGURATION
# ============================================================

BASE_DATA_FILE = "ILD_TOP_DRIVERS_DATA.csv"
TARGET_FILE = "REPHENOTYPES FOR IC.csv"
PRED_FILE = "PREDICTIONS_104W_PRED_WINDOW.parquet"

TARGET_NAME = "FILD or FILA ADJUDICATED"

GG_COLUMN = "rs35705950.1_G_2"
GT_COLUMN = "rs35705950.1_G_1"
TT_COLUMN = "rs35705950.1_G_0"

TRAIN_SIZE = 0.40
N_REPEATS = 100

# ------------------------------------------------------------
# FIXED exploratory regions motivated by the previously observed
# local MUC5B-vs-ZeBRA curves.
#
# Within these regions, exact ZeBRA score is ignored and risk is
# estimated using MUC5B carrier status within that region.
#
# Everywhere else, calibrated ZeBRA is used.
# ------------------------------------------------------------

MUC5B_SWITCH_WINDOWS = [
    (0.0, 50.0, "P00_P50"),
    (65.0, 75.0, "P65_P75"),
    (85.0, 95.0, "P85_P95"),
]

EVAL_TARGET_FPRS = [
    0.005,
    0.010,
    0.020,
    0.050,
    0.100,
]

# Beta-binomial / Jeffreys smoothing for P(Y|region,G).
ALPHA = 0.5
BETA = 0.5

OUTDIR = Path(
    "./RESULTS/FIXED_MUC5B_SWITCH_REGIONS"
)
OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# HELPERS
# ============================================================

def decode_muc5b(df):
    required = [
        GG_COLUMN,
        GT_COLUMN,
        TT_COLUMN,
    ]

    missing = [
        c for c in required
        if c not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing expected MUC5B columns:\n"
            + "\n".join(missing)
        )

    d = df.copy()

    for c in required:
        d[c] = pd.to_numeric(
            d[c],
            errors="coerce",
        )

    onehot_sum = d[
        required
    ].sum(
        axis=1
    )

    print()
    print("MUC5B dummy-state sums:")
    print(
        onehot_sum
        .value_counts()
        .sort_index()
        .to_string()
    )

    if (
        onehot_sum > 1
    ).any():
        raise ValueError(
            "Some rs35705950 rows are multi-hot."
        )

    # all-zero = missing original genotype
    d = d.loc[
        onehot_sum == 1
    ].copy()

    d["MUC5B_GG"] = (
        d[GG_COLUMN]
        .astype(int)
    )

    d["MUC5B_GT"] = (
        d[GT_COLUMN]
        .astype(int)
    )

    d["MUC5B_TT"] = (
        d[TT_COLUMN]
        .astype(int)
    )

    d["MUC5B_T_carrier"] = (
        (
            d["MUC5B_GT"] == 1
        )
        |
        (
            d["MUC5B_TT"] == 1
        )
    ).astype(int)

    d["MUC5B_T_dosage"] = (
        d["MUC5B_GT"]
        + 2 * d["MUC5B_TT"]
    )

    return d


def empirical_percentile(
    reference_scores,
    query_scores,
):
    """
    Percentile of query scores relative to the TRAINING reference
    distribution.

    This avoids using the test distribution to define percentiles.
    """
    ref = np.sort(
        np.asarray(
            reference_scores,
            dtype=float,
        )
    )

    q = np.asarray(
        query_scores,
        dtype=float,
    )

    ranks = np.searchsorted(
        ref,
        q,
        side="right",
    )

    return (
        100.0
        * ranks
        / len(ref)
    )


def switch_region_id(
    percentiles,
):
    p = np.asarray(
        percentiles,
        dtype=float,
    )

    region = np.full(
        len(p),
        "",
        dtype=object,
    )

    for lo, hi, name in MUC5B_SWITCH_WINDOWS:
        mask = (
            (p >= lo)
            &
            (p < hi)
        )

        region[
            mask
        ] = name

    return region


def estimate_region_genotype_risks(
    y_train,
    train_percentile,
    g_train,
):
    """
    Estimate smoothed P(Y=1 | percentile region, MUC5B carrier)
    using TRAINING data only.
    """
    rows = []
    lookup = {}

    for lo, hi, name in MUC5B_SWITCH_WINDOWS:

        region_mask = (
            (train_percentile >= lo)
            &
            (train_percentile < hi)
        )

        for g in [
            0,
            1,
        ]:
            mask = (
                region_mask
                &
                (g_train == g)
            )

            n = int(
                mask.sum()
            )

            cases = int(
                y_train[
                    mask
                ].sum()
            )

            risk = (
                cases + ALPHA
            ) / (
                n
                + ALPHA
                + BETA
            )

            controls = (
                n - cases
            )

            rows.append({
                "region":
                    name,

                "percentile_low":
                    lo,

                "percentile_high":
                    hi,

                "MUC5B_carrier":
                    g,

                "N":
                    n,

                "cases":
                    cases,

                "controls":
                    controls,

                "smoothed_case_probability":
                    risk,
            })

            lookup[
                (
                    name,
                    g,
                )
            ] = float(
                risk
            )

    table = pd.DataFrame(
        rows
    )

    # --------------------------------------------------------
    # Add within-region MUC5B OR diagnostics.
    # --------------------------------------------------------

    or_rows = []

    for lo, hi, name in MUC5B_SWITCH_WINDOWS:

        region_mask = (
            (train_percentile >= lo)
            &
            (train_percentile < hi)
        )

        yr = y_train[
            region_mask
        ]

        gr = g_train[
            region_mask
        ]

        a = int(
            np.sum(
                (yr == 1)
                &
                (gr == 1)
            )
        )

        b = int(
            np.sum(
                (yr == 0)
                &
                (gr == 1)
            )
        )

        c = int(
            np.sum(
                (yr == 1)
                &
                (gr == 0)
            )
        )

        d = int(
            np.sum(
                (yr == 0)
                &
                (gr == 0)
            )
        )

        fisher_or, fisher_p = fisher_exact(
            [
                [a, b],
                [c, d],
            ]
        )

        # Haldane correction for finite descriptive OR.
        aa, bb, cc, dd = map(
            float,
            [
                a,
                b,
                c,
                d,
            ],
        )

        if min(
            aa,
            bb,
            cc,
            dd,
        ) == 0:
            aa += 0.5
            bb += 0.5
            cc += 0.5
            dd += 0.5

        corrected_or = (
            aa * dd
            / (
                bb * cc
            )
        )

        or_rows.append({
            "region":
                name,

            "N":
                int(
                    region_mask.sum()
                ),

            "cases":
                int(
                    yr.sum()
                ),

            "carrier_case_N":
                a,

            "carrier_control_N":
                b,

            "noncarrier_case_N":
                c,

            "noncarrier_control_N":
                d,

            "MUC5B_OR":
                corrected_or,

            "fisher_OR":
                fisher_or,

            "fisher_p":
                fisher_p,
        })

    return (
        lookup,
        table,
        pd.DataFrame(
            or_rows
        ),
    )


def fit_zebra_calibrator(
    z_train,
    y_train,
):
    """
    Training-only monotone calibration of ZeBRA to disease probability.
    """
    iso = IsotonicRegression(
        y_min=0.0,
        y_max=1.0,
        out_of_bounds="clip",
    )

    iso.fit(
        z_train,
        y_train,
    )

    return iso


def make_hybrid_score(
    z,
    percentile,
    g,
    zebra_calibrator,
    region_genotype_lookup,
):
    """
    Outside switch windows:
        score = calibrated ZeBRA P(Y|Z)

    Inside switch windows:
        score = training-estimated P(Y|region,MUC5B)

    Thus exact ZeBRA ranking is deliberately discarded inside
    the specified genomic-dominant windows.
    """
    z = np.asarray(
        z,
        dtype=float,
    )

    percentile = np.asarray(
        percentile,
        dtype=float,
    )

    g = np.asarray(
        g,
        dtype=int,
    )

    score = zebra_calibrator.predict(
        z
    ).astype(float)

    region = switch_region_id(
        percentile
    )

    for i in range(
        len(score)
    ):
        if region[i] != "":
            score[i] = (
                region_genotype_lookup[
                    (
                        region[i],
                        int(
                            g[i]
                        ),
                    )
                ]
            )

    return (
        score,
        region,
    )


def threshold_for_target_fpr(
    negative_scores,
    target_fpr,
):
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
            negative_scores
            >= threshold
        )
    )

    return (
        threshold,
        achieved,
    )


def threshold_matching_empirical_fpr(
    negative_scores,
    desired_fpr,
):
    """
    Select threshold on TRAINING negative scores whose empirical
    FPR is closest to desired_fpr.
    """
    x = np.asarray(
        negative_scores,
        dtype=float,
    )

    unique_scores = np.unique(
        x
    )

    candidate_fprs = np.array(
        [
            np.mean(
                x >= t
            )
            for t in unique_scores
        ],
        dtype=float,
    )

    idx = int(
        np.argmin(
            np.abs(
                candidate_fprs
                - desired_fpr
            )
        )
    )

    return (
        float(
            unique_scores[idx]
        ),
        float(
            candidate_fprs[idx]
        ),
    )


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
            (y == 1)
            &
            positive
        )
    )

    fn = int(
        np.sum(
            (y == 1)
            &
            (~positive)
        )
    )

    fp = int(
        np.sum(
            (y == 0)
            &
            positive
        )
    )

    tn = int(
        np.sum(
            (y == 0)
            &
            (~positive)
        )
    )

    sensitivity = (
        tp / (
            tp + fn
        )
        if (
            tp + fn
        ) > 0
        else np.nan
    )

    specificity = (
        tn / (
            tn + fp
        )
        if (
            tn + fp
        ) > 0
        else np.nan
    )

    fpr = (
        fp / (
            fp + tn
        )
        if (
            fp + tn
        ) > 0
        else np.nan
    )

    ppv = (
        tp / (
            tp + fp
        )
        if (
            tp + fp
        ) > 0
        else np.nan
    )

    npv = (
        tn / (
            tn + fn
        )
        if (
            tn + fn
        ) > 0
        else np.nan
    )

    lr_pos = (
        sensitivity / fpr
        if (
            np.isfinite(
                sensitivity
            )
            and np.isfinite(
                fpr
            )
            and fpr > 0
        )
        else np.inf
    )

    lr_neg = (
        (
            1.0
            - sensitivity
        )
        / specificity
        if (
            np.isfinite(
                sensitivity
            )
            and np.isfinite(
                specificity
            )
            and specificity > 0
        )
        else np.nan
    )

    dor = (
        lr_pos / lr_neg
        if (
            np.isfinite(
                lr_pos
            )
            and np.isfinite(
                lr_neg
            )
            and lr_neg > 0
        )
        else np.nan
    )

    return {
        "TP":
            tp,

        "FN":
            fn,

        "FP":
            fp,

        "TN":
            tn,

        "sensitivity":
            sensitivity,

        "specificity":
            specificity,

        "FPR":
            fpr,

        "PPV":
            ppv,

        "NPV":
            npv,

        "LR_positive":
            lr_pos,

        "LR_negative":
            lr_neg,

        "diagnostic_odds_ratio":
            dor,
    }


def summarize_repeats(
    df,
    group_cols,
    value_cols,
):
    rows = []

    for keys, sub in df.groupby(
        group_cols
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
                group_cols,
                keys,
            )
        )

        row[
            "n_repeats"
        ] = len(
            sub
        )

        for col in value_cols:

            vals = (
                pd.to_numeric(
                    sub[
                        col
                    ],
                    errors="coerce",
                )
                .replace(
                    [
                        np.inf,
                        -np.inf,
                    ],
                    np.nan,
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
# LOAD DATA
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

DATA = DATA[
    DATA[
        "predicted_risk"
    ].notnull()
].copy()

DATA = decode_muc5b(
    DATA
)


# ============================================================
# COHORT REPORT
# ============================================================

print()
print("=" * 96)
print("COHORT")
print("=" * 96)

print(
    f"N = {len(DATA):,}"
)

print(
    f"FILD/FILA positives = "
    f"{int(DATA['target'].sum()):,}"
)

print(
    f"MUC5B carriers = "
    f"{int(DATA['MUC5B_T_carrier'].sum()):,}"
)

print(
    f"Raw ZeBRA AUC = "
    f"{roc_auc_score(DATA['target'], DATA['predicted_risk']):.4f}"
)

print()
print(
    "Fixed MUC5B-switch percentile windows:"
)

for lo, hi, name in MUC5B_SWITCH_WINDOWS:
    print(
        f"  {name}: "
        f"{lo:g} <= percentile < {hi:g}"
    )


# ============================================================
# REPEATED OUTER-SPLIT TEST
# ============================================================

indices = np.arange(
    len(
        DATA
    )
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

g_all = DATA[
    "MUC5B_T_carrier"
].to_numpy(
    dtype=int
)

RESULT_ROWS = []
OR_ROWS = []
RISK_ROWS = []

print()
print("=" * 96)
print(
    f"RUNNING {N_REPEATS} REPEATED OUTER SPLITS"
)
print("=" * 96)

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

    z_train = z_all[
        train_idx
    ]

    g_train = g_all[
        train_idx
    ]

    y_test = y_all[
        test_idx
    ]

    z_test = z_all[
        test_idx
    ]

    g_test = g_all[
        test_idx
    ]

    # --------------------------------------------------------
    # Percentiles are ALWAYS defined from training distribution.
    # --------------------------------------------------------

    train_pct = empirical_percentile(
        z_train,
        z_train,
    )

    test_pct = empirical_percentile(
        z_train,
        z_test,
    )

    # --------------------------------------------------------
    # Fit ZeBRA calibration and MUC5B regional risks on training.
    # --------------------------------------------------------

    zebra_calibrator = fit_zebra_calibrator(
        z_train,
        y_train,
    )

    (
        region_genotype_lookup,
        region_risk_table,
        region_or_table,
    ) = estimate_region_genotype_risks(
        y_train,
        train_pct,
        g_train,
    )

    region_risk_table.insert(
        0,
        "repeat",
        repeat,
    )

    region_risk_table.insert(
        1,
        "split_seed",
        split_seed,
    )

    RISK_ROWS.append(
        region_risk_table
    )

    region_or_table.insert(
        0,
        "repeat",
        repeat,
    )

    region_or_table.insert(
        1,
        "split_seed",
        split_seed,
    )

    OR_ROWS.append(
        region_or_table
    )

    # --------------------------------------------------------
    # Construct hybrid score.
    # --------------------------------------------------------

    train_hybrid, train_region = make_hybrid_score(
        z_train,
        train_pct,
        g_train,
        zebra_calibrator,
        region_genotype_lookup,
    )

    test_hybrid, test_region = make_hybrid_score(
        z_test,
        test_pct,
        g_test,
        zebra_calibrator,
        region_genotype_lookup,
    )

    # Comparator 1: raw ZeBRA.
    train_zebra_score = z_train
    test_zebra_score = z_test

    # Comparator 2: calibrated ZeBRA only.
    train_zebra_cal = (
        zebra_calibrator.predict(
            z_train
        )
    )

    test_zebra_cal = (
        zebra_calibrator.predict(
            z_test
        )
    )

    auc_raw = roc_auc_score(
        y_test,
        test_zebra_score,
    )

    auc_cal = roc_auc_score(
        y_test,
        test_zebra_cal,
    )

    auc_hybrid = roc_auc_score(
        y_test,
        test_hybrid,
    )

    # --------------------------------------------------------
    # Operating-point tests.
    #
    # We choose a HYBRID threshold at nominal training FPR,
    # then choose a RAW-ZeBRA threshold whose TRAINING FPR
    # matches the achieved hybrid FPR as closely as possible.
    #
    # This gives a fair matched-FPR comparison.
    # --------------------------------------------------------

    train_neg = (
        y_train == 0
    )

    test_neg = (
        y_test == 0
    )

    for target_fpr in EVAL_TARGET_FPRS:

        hybrid_threshold, hybrid_train_fpr = (
            threshold_for_target_fpr(
                train_hybrid[
                    train_neg
                ],
                target_fpr,
            )
        )

        raw_threshold, raw_train_fpr = (
            threshold_matching_empirical_fpr(
                train_zebra_score[
                    train_neg
                ],
                hybrid_train_fpr,
            )
        )

        cal_threshold, cal_train_fpr = (
            threshold_matching_empirical_fpr(
                train_zebra_cal[
                    train_neg
                ],
                hybrid_train_fpr,
            )
        )

        hybrid_positive = (
            test_hybrid
            >= hybrid_threshold
        )

        raw_positive = (
            test_zebra_score
            >= raw_threshold
        )

        cal_positive = (
            test_zebra_cal
            >= cal_threshold
        )

        mh = confusion_metrics(
            y_test,
            hybrid_positive,
        )

        mr = confusion_metrics(
            y_test,
            raw_positive,
        )

        mc = confusion_metrics(
            y_test,
            cal_positive,
        )

        RESULT_ROWS.append({
            "repeat":
                repeat,

            "split_seed":
                split_seed,

            "target_fpr":
                target_fpr,

            "auc_zebra_raw":
                auc_raw,

            "auc_zebra_calibrated":
                auc_cal,

            "auc_hybrid":
                auc_hybrid,

            "delta_auc_hybrid_vs_raw":
                (
                    auc_hybrid
                    - auc_raw
                ),

            "hybrid_train_threshold":
                hybrid_threshold,

            "hybrid_train_achieved_FPR":
                hybrid_train_fpr,

            "raw_zebra_train_threshold":
                raw_threshold,

            "raw_zebra_train_achieved_FPR":
                raw_train_fpr,

            "cal_zebra_train_threshold":
                cal_threshold,

            "cal_zebra_train_achieved_FPR":
                cal_train_fpr,

            "hybrid_sensitivity":
                mh[
                    "sensitivity"
                ],

            "hybrid_FPR":
                mh[
                    "FPR"
                ],

            "hybrid_PPV":
                mh[
                    "PPV"
                ],

            "hybrid_LR_positive":
                mh[
                    "LR_positive"
                ],

            "hybrid_LR_negative":
                mh[
                    "LR_negative"
                ],

            "hybrid_DOR":
                mh[
                    "diagnostic_odds_ratio"
                ],

            "raw_zebra_sensitivity":
                mr[
                    "sensitivity"
                ],

            "raw_zebra_FPR":
                mr[
                    "FPR"
                ],

            "raw_zebra_PPV":
                mr[
                    "PPV"
                ],

            "raw_zebra_LR_positive":
                mr[
                    "LR_positive"
                ],

            "raw_zebra_LR_negative":
                mr[
                    "LR_negative"
                ],

            "raw_zebra_DOR":
                mr[
                    "diagnostic_odds_ratio"
                ],

            "cal_zebra_sensitivity":
                mc[
                    "sensitivity"
                ],

            "cal_zebra_FPR":
                mc[
                    "FPR"
                ],

            "cal_zebra_PPV":
                mc[
                    "PPV"
                ],

            "cal_zebra_LR_positive":
                mc[
                    "LR_positive"
                ],

            "cal_zebra_LR_negative":
                mc[
                    "LR_negative"
                ],

            "delta_sensitivity_hybrid_vs_raw":
                (
                    mh[
                        "sensitivity"
                    ]
                    - mr[
                        "sensitivity"
                    ]
                ),

            "delta_FPR_hybrid_vs_raw":
                (
                    mh[
                        "FPR"
                    ]
                    - mr[
                        "FPR"
                    ]
                ),

            "delta_PPV_hybrid_vs_raw":
                (
                    mh[
                        "PPV"
                    ]
                    - mr[
                        "PPV"
                    ]
                ),

            "delta_LRplus_hybrid_vs_raw":
                (
                    mh[
                        "LR_positive"
                    ]
                    - mr[
                        "LR_positive"
                    ]
                    if (
                        np.isfinite(
                            mh[
                                "LR_positive"
                            ]
                        )
                        and np.isfinite(
                            mr[
                                "LR_positive"
                            ]
                        )
                    )
                    else np.nan
                ),

            "delta_LRminus_hybrid_vs_raw":
                (
                    mh[
                        "LR_negative"
                    ]
                    - mr[
                        "LR_negative"
                    ]
                ),
        })

    if (
        repeat % 10
    ) == 0:
        print(
            f"completed repeat "
            f"{repeat}/{N_REPEATS}"
        )


# ============================================================
# SAVE RAW RESULTS
# ============================================================

RESULT_DF = pd.DataFrame(
    RESULT_ROWS
)

OR_DF = pd.concat(
    OR_ROWS,
    ignore_index=True,
)

RISK_DF = pd.concat(
    RISK_ROWS,
    ignore_index=True,
)

RESULT_DF.to_csv(
    OUTDIR
    / "REPEATED_SPLIT_FIXED_SWITCH_RESULTS.csv",
    index=False,
)

OR_DF.to_csv(
    OUTDIR
    / "TRAINING_LOCAL_MUC5B_OR_BY_SWITCH_REGION.csv",
    index=False,
)

RISK_DF.to_csv(
    OUTDIR
    / "TRAINING_REGION_MUC5B_RISK_ESTIMATES.csv",
    index=False,
)


# ============================================================
# SUMMARIZE PERFORMANCE
# ============================================================

VALUE_COLS = [
    "auc_zebra_raw",
    "auc_zebra_calibrated",
    "auc_hybrid",
    "delta_auc_hybrid_vs_raw",

    "hybrid_sensitivity",
    "hybrid_FPR",
    "hybrid_PPV",
    "hybrid_LR_positive",
    "hybrid_LR_negative",
    "hybrid_DOR",

    "raw_zebra_sensitivity",
    "raw_zebra_FPR",
    "raw_zebra_PPV",
    "raw_zebra_LR_positive",
    "raw_zebra_LR_negative",
    "raw_zebra_DOR",

    "cal_zebra_sensitivity",
    "cal_zebra_FPR",
    "cal_zebra_PPV",
    "cal_zebra_LR_positive",
    "cal_zebra_LR_negative",

    "delta_sensitivity_hybrid_vs_raw",
    "delta_FPR_hybrid_vs_raw",
    "delta_PPV_hybrid_vs_raw",
    "delta_LRplus_hybrid_vs_raw",
    "delta_LRminus_hybrid_vs_raw",
]

SUMMARY = summarize_repeats(
    RESULT_DF,
    [
        "target_fpr",
    ],
    VALUE_COLS,
)

SUMMARY.to_csv(
    OUTDIR
    / "FIXED_SWITCH_VS_ZEBRA_SUMMARY.csv",
    index=False,
)


# ============================================================
# SUMMARIZE LOCAL MUC5B OR
# ============================================================

OR_SUMMARY = (
    OR_DF
    .groupby(
        "region"
    )
    .agg(
        n_repeats=(
            "repeat",
            "nunique",
        ),
        mean_N=(
            "N",
            "mean",
        ),
        mean_cases=(
            "cases",
            "mean",
        ),
        mean_MUC5B_OR=(
            "MUC5B_OR",
            "mean",
        ),
        median_MUC5B_OR=(
            "MUC5B_OR",
            "median",
        ),
        q025_MUC5B_OR=(
            "MUC5B_OR",
            lambda x:
                np.quantile(
                    x,
                    0.025,
                ),
        ),
        q975_MUC5B_OR=(
            "MUC5B_OR",
            lambda x:
                np.quantile(
                    x,
                    0.975,
                ),
        ),
        fraction_fisher_p_lt_005=(
            "fisher_p",
            lambda x:
                np.mean(
                    np.asarray(
                        x
                    )
                    < 0.05
                ),
        ),
    )
    .reset_index()
)

OR_SUMMARY.to_csv(
    OUTDIR
    / "LOCAL_MUC5B_OR_SWITCH_REGION_SUMMARY.csv",
    index=False,
)


# ============================================================
# CONSOLE REPORT
# ============================================================

print()
print("=" * 96)
print("FIXED-SWITCH PERFORMANCE SUMMARY")
print("=" * 96)

display_cols = [
    "target_fpr",
    "n_repeats",

    "auc_zebra_raw_mean",
    "auc_hybrid_mean",
    "delta_auc_hybrid_vs_raw_mean",

    "raw_zebra_sensitivity_mean",
    "hybrid_sensitivity_mean",
    "delta_sensitivity_hybrid_vs_raw_mean",
    "delta_sensitivity_hybrid_vs_raw_q025",
    "delta_sensitivity_hybrid_vs_raw_q975",
    "delta_sensitivity_hybrid_vs_raw_fraction_gt_zero",

    "raw_zebra_FPR_mean",
    "hybrid_FPR_mean",
    "delta_FPR_hybrid_vs_raw_mean",

    "raw_zebra_LR_positive_mean",
    "hybrid_LR_positive_mean",
    "delta_LRplus_hybrid_vs_raw_mean",

    "raw_zebra_LR_negative_mean",
    "hybrid_LR_negative_mean",
    "delta_LRminus_hybrid_vs_raw_mean",

    "raw_zebra_PPV_mean",
    "hybrid_PPV_mean",
]

print(
    SUMMARY[
        [
            c
            for c in display_cols
            if c in SUMMARY.columns
        ]
    ].to_string(
        index=False
    )
)

print()
print("=" * 96)
print("MUC5B OR IN THE FIXED SWITCH REGIONS")
print("=" * 96)

print(
    OR_SUMMARY.to_string(
        index=False
    )
)


# ============================================================
# PLOTS
# ============================================================

# ------------------------------------------------------------
# Plot 1: sensitivity at matched FPR
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(9, 6)
)

ax.plot(
    SUMMARY[
        "target_fpr"
    ],
    SUMMARY[
        "raw_zebra_sensitivity_mean"
    ],
    marker="o",
    label="ZeBRA alone",
)

ax.plot(
    SUMMARY[
        "target_fpr"
    ],
    SUMMARY[
        "hybrid_sensitivity_mean"
    ],
    marker="o",
    label="Fixed MUC5B-switch hybrid",
)

ax.set_xlabel(
    "Target FPR"
)

ax.set_ylabel(
    "Mean test sensitivity"
)

ax.set_title(
    "Fixed-region MUC5B switch vs ZeBRA at matched FPR"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MATCHED_FPR_SENSITIVITY.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ------------------------------------------------------------
# Plot 2: LR+
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(9, 6)
)

ax.plot(
    SUMMARY[
        "target_fpr"
    ],
    SUMMARY[
        "raw_zebra_LR_positive_mean"
    ],
    marker="o",
    label="ZeBRA alone",
)

ax.plot(
    SUMMARY[
        "target_fpr"
    ],
    SUMMARY[
        "hybrid_LR_positive_mean"
    ],
    marker="o",
    label="Fixed MUC5B-switch hybrid",
)

ax.set_xlabel(
    "Target FPR"
)

ax.set_ylabel(
    "Mean LR+"
)

ax.set_title(
    "Positive likelihood ratio: fixed MUC5B-switch hybrid"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MATCHED_FPR_LRPLUS.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ------------------------------------------------------------
# Plot 3: LR-
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(9, 6)
)

ax.plot(
    SUMMARY[
        "target_fpr"
    ],
    SUMMARY[
        "raw_zebra_LR_negative_mean"
    ],
    marker="o",
    label="ZeBRA alone",
)

ax.plot(
    SUMMARY[
        "target_fpr"
    ],
    SUMMARY[
        "hybrid_LR_negative_mean"
    ],
    marker="o",
    label="Fixed MUC5B-switch hybrid",
)

ax.set_xlabel(
    "Target FPR"
)

ax.set_ylabel(
    "Mean LR-"
)

ax.set_title(
    "Negative likelihood ratio: fixed MUC5B-switch hybrid"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "MATCHED_FPR_LRMINUS.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ------------------------------------------------------------
# Plot 4: local MUC5B OR in specified switch regions
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(8, 5)
)

x = np.arange(
    len(
        OR_SUMMARY
    )
)

ax.bar(
    x,
    OR_SUMMARY[
        "median_MUC5B_OR"
    ],
)

ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_xticks(
    x
)

ax.set_xticklabels(
    OR_SUMMARY[
        "region"
    ],
)

ax.set_ylabel(
    "Median training-split MUC5B odds ratio"
)

ax.set_title(
    "MUC5B disease association in fixed ZeBRA-percentile switch regions"
)

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "LOCAL_MUC5B_OR_SWITCH_REGIONS.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# INTERPRETATION GUIDE
# ============================================================

print()
print("=" * 96)
print("INTERPRETATION GUIDE")
print("=" * 96)

print(
    "Hybrid rule:"
)

print(
    "  ZeBRA percentile <50:       ignore exact ZeBRA; use MUC5B regional risk"
)

print(
    "  ZeBRA percentile 50-65:     use calibrated ZeBRA"
)

print(
    "  ZeBRA percentile 65-75:     ignore exact ZeBRA; use MUC5B regional risk"
)

print(
    "  ZeBRA percentile 75-85:     use calibrated ZeBRA"
)

print(
    "  ZeBRA percentile 85-95:     ignore exact ZeBRA; use MUC5B regional risk"
)

print(
    "  ZeBRA percentile >=95:      use calibrated ZeBRA"
)

print()
print(
    "The decisive comparison is hybrid vs raw ZeBRA at the SAME "
    "training-derived achieved FPR."
)

print()
print(
    "Support for the proposed switch rule requires:"
)

print(
    "  1. delta_sensitivity_hybrid_vs_raw > 0,"
)

print(
    "  2. test FPRs remain closely matched,"
)

print(
    "  3. hybrid LR+ > ZeBRA LR+,"
)

print(
    "  4. hybrid LR- < ZeBRA LR-,"
)

print(
    "  5. MUC5B OR >1 in the proposed switch regions."
)

print()
print(
    "IMPORTANT: these percentile regions were selected after inspecting "
    "the full-cohort local-information plots. Therefore this is an "
    "exploratory/post-hoc validation, not an independent confirmatory test. "
    "A subsequent external cohort or a new pre-specified split is needed "
    "for confirmatory inference."
)

print()
print(
    "Outputs saved under:"
)

print(
    OUTDIR.resolve()
)
