#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import chi2, norm
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

OUTDIR = Path(
    "./RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES"
)
OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

# Moving windows along the empirical ZeBRA risk distribution.
#
# A 10-percentile-point window moved every 2.5 percentile points
# gives substantial overlap and a reasonably smooth descriptive curve.
WINDOW_HALF_WIDTH_PERCENTILE = 5.0

WINDOW_CENTERS = np.arange(
    5.0,
    97.6,
    2.5,
)

MIN_WINDOW_N = 150
MIN_CASES = 8
MIN_CONTROLS = 20

# Very small ridge used only for numerical stabilization of MLE.
RIDGE = 1e-8

# A practical threshold for flagging windows where both local
# partial-information statistics are weak.
WEAK_LRT_THRESHOLD = 1.0


# ============================================================
# Logistic MLE using Newton/IRLS
#
# Avoids needing statsmodels.
# ============================================================

def sigmoid(x):
    x = np.clip(
        x,
        -35,
        35,
    )
    return 1.0 / (
        1.0 + np.exp(-x)
    )


def logistic_fit(
    X,
    y,
    ridge=RIDGE,
    max_iter=200,
    tol=1e-9,
):
    """
    Logistic-regression maximum likelihood fit with an unpenalized
    intercept and tiny numerical ridge on non-intercept terms.

    Returns coefficients, covariance matrix, fitted probabilities,
    and the ordinary (unpenalized) Bernoulli log-likelihood.
    """
    X = np.asarray(
        X,
        dtype=float,
    )

    y = np.asarray(
        y,
        dtype=float,
    )

    n, p = X.shape

    beta = np.zeros(
        p,
        dtype=float,
    )

    penalty = np.eye(
        p,
        dtype=float,
    )

    # Do not penalize intercept.
    penalty[
        0,
        0,
    ] = 0.0

    for _ in range(
        max_iter
    ):
        eta = X @ beta

        prob = sigmoid(
            eta
        )

        w = np.clip(
            prob * (
                1.0 - prob
            ),
            1e-8,
            None,
        )

        gradient = (
            X.T @ (
                y - prob
            )
            - ridge
            * penalty
            @ beta
        )

        hessian_pos = (
            X.T
            @ (
                X
                * w[:, None]
            )
            + ridge
            * penalty
        )

        try:
            step = np.linalg.solve(
                hessian_pos,
                gradient,
            )
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(
                hessian_pos
            ) @ gradient

        beta_new = (
            beta
            + step
        )

        if np.max(
            np.abs(
                beta_new
                - beta
            )
        ) < tol:
            beta = beta_new
            break

        beta = beta_new

    eta = X @ beta

    prob = sigmoid(
        eta
    )

    w = np.clip(
        prob
        * (
            1.0 - prob
        ),
        1e-8,
        None,
    )

    fisher = (
        X.T
        @ (
            X
            * w[:, None]
        )
        + ridge
        * penalty
    )

    try:
        cov = np.linalg.inv(
            fisher
        )
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(
            fisher
        )

    prob_clip = np.clip(
        prob,
        1e-12,
        1.0 - 1e-12,
    )

    ll = float(
        np.sum(
            y
            * np.log(
                prob_clip
            )
            + (
                1.0 - y
            )
            * np.log(
                1.0 - prob_clip
            )
        )
    )

    return {
        "beta":
            beta,

        "cov":
            cov,

        "prob":
            prob,

        "ll":
            ll,
    }


def coef_summary(
    fit,
    idx,
):
    beta = float(
        fit[
            "beta"
        ][idx]
    )

    var = float(
        fit[
            "cov"
        ][
            idx,
            idx,
        ]
    )

    se = np.sqrt(
        max(
            var,
            0.0,
        )
    )

    if se > 0:
        z = beta / se

        p = float(
            2.0
            * norm.sf(
                abs(
                    z
                )
            )
        )

        lo = beta - 1.959963984540054 * se
        hi = beta + 1.959963984540054 * se

    else:
        z = np.nan
        p = np.nan
        lo = np.nan
        hi = np.nan

    return {
        "beta":
            beta,

        "se":
            se,

        "z":
            z,

        "p":
            p,

        "OR":
            float(
                np.exp(
                    beta
                )
            ),

        "OR_ci_low":
            float(
                np.exp(
                    lo
                )
            )
            if np.isfinite(
                lo
            )
            else np.nan,

        "OR_ci_high":
            float(
                np.exp(
                    hi
                )
            )
            if np.isfinite(
                hi
            )
            else np.nan,
    }


def lrt(
    ll_large,
    ll_small,
    df=1,
):
    stat = max(
        2.0
        * (
            ll_large
            - ll_small
        ),
        0.0,
    )

    p = float(
        chi2.sf(
            stat,
            df,
        )
    )

    return (
        float(
            stat
        ),
        p,
    )


# ============================================================
# Genotype decoding
# ============================================================

def decode_muc5b(
    df,
):
    required = [
        GG_COLUMN,
        GT_COLUMN,
        TT_COLUMN,
    ]

    missing = [
        c
        for c in required
        if c not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing expected MUC5B columns:\n"
            + "\n".join(
                missing
            )
        )

    d = df.copy()

    for c in required:
        d[
            c
        ] = pd.to_numeric(
            d[
                c
            ],
            errors="coerce",
        )

    dummy_sum = d[
        required
    ].sum(
        axis=1
    )

    print()
    print(
        "MUC5B dummy-state sums:"
    )
    print(
        dummy_sum
        .value_counts()
        .sort_index()
        .to_string()
    )

    if (
        dummy_sum > 1
    ).any():
        raise ValueError(
            "Some MUC5B rows have more than one active genotype state."
        )

    # sum==0 is missing original genotype after pd.get_dummies().
    d = d.loc[
        dummy_sum == 1
    ].copy()

    d[
        "MUC5B_GG"
    ] = d[
        GG_COLUMN
    ].astype(int)

    d[
        "MUC5B_GT"
    ] = d[
        GT_COLUMN
    ].astype(int)

    d[
        "MUC5B_TT"
    ] = d[
        TT_COLUMN
    ].astype(int)

    d[
        "MUC5B_T_carrier"
    ] = (
        (
            d[
                "MUC5B_GT"
            ] == 1
        )
        |
        (
            d[
                "MUC5B_TT"
            ] == 1
        )
    ).astype(int)

    d[
        "MUC5B_T_dosage"
    ] = (
        d[
            "MUC5B_GT"
        ]
        + 2
        * d[
            "MUC5B_TT"
        ]
    )

    return d


# ============================================================
# Local-window analysis
# ============================================================

def analyze_window(
    d,
    center,
    half_width,
):
    lo_pct = (
        center
        - half_width
    )

    hi_pct = (
        center
        + half_width
    )

    mask = (
        (
            d[
                "ZeBRA_percentile"
            ]
            >= lo_pct
        )
        &
        (
            d[
                "ZeBRA_percentile"
            ]
            <= hi_pct
        )
    )

    w = d.loc[
        mask
    ].copy()

    n = len(
        w
    )

    if n < MIN_WINDOW_N:
        return None

    y = (
        w[
            "target"
        ]
        .astype(int)
        .to_numpy()
    )

    z_raw = (
        w[
            "predicted_risk"
        ]
        .astype(float)
        .to_numpy()
    )

    g = (
        w[
            "MUC5B_T_carrier"
        ]
        .astype(int)
        .to_numpy()
    )

    cases = int(
        y.sum()
    )

    controls = int(
        len(y)
        - cases
    )

    if (
        cases < MIN_CASES
        or controls < MIN_CONTROLS
        or np.unique(
            g
        ).size < 2
    ):
        return None

    z_sd = float(
        np.std(
            z_raw,
            ddof=0,
        )
    )

    z_mean = float(
        np.mean(
            z_raw
        )
    )

    if z_sd > 1e-12:
        z = (
            z_raw
            - z_mean
        ) / z_sd
    else:
        # In score-saturation windows, ZeBRA has no local variation.
        z = np.zeros_like(
            z_raw,
            dtype=float,
        )

    intercept = np.ones(
        len(y),
        dtype=float,
    )

    X_null = np.column_stack(
        [
            intercept,
        ]
    )

    X_z = np.column_stack(
        [
            intercept,
            z,
        ]
    )

    X_g = np.column_stack(
        [
            intercept,
            g,
        ]
    )

    X_both = np.column_stack(
        [
            intercept,
            z,
            g,
        ]
    )

    X_interaction = np.column_stack(
        [
            intercept,
            z,
            g,
            z * g,
        ]
    )

    fit_null = logistic_fit(
        X_null,
        y,
    )

    fit_z = logistic_fit(
        X_z,
        y,
    )

    fit_g = logistic_fit(
        X_g,
        y,
    )

    fit_both = logistic_fit(
        X_both,
        y,
    )

    fit_interaction = logistic_fit(
        X_interaction,
        y,
    )

    # --------------------------------------------------------
    # Comparable local predictive-information measures
    #
    # Conditional information in ZeBRA given genotype:
    #    2*(LL[Z+G] - LL[G])
    #
    # Conditional information in genotype given ZeBRA:
    #    2*(LL[Z+G] - LL[Z])
    #
    # Both are 1-df LRT statistics.
    # --------------------------------------------------------

    lrt_z_given_g, p_z_given_g = lrt(
        fit_both[
            "ll"
        ],
        fit_g[
            "ll"
        ],
        df=1,
    )

    lrt_g_given_z, p_g_given_z = lrt(
        fit_both[
            "ll"
        ],
        fit_z[
            "ll"
        ],
        df=1,
    )

    lrt_z_marginal, p_z_marginal = lrt(
        fit_z[
            "ll"
        ],
        fit_null[
            "ll"
        ],
        df=1,
    )

    lrt_g_marginal, p_g_marginal = lrt(
        fit_g[
            "ll"
        ],
        fit_null[
            "ll"
        ],
        df=1,
    )

    lrt_interaction_stat, p_interaction = lrt(
        fit_interaction[
            "ll"
        ],
        fit_both[
            "ll"
        ],
        df=1,
    )

    z_effect = coef_summary(
        fit_both,
        1,
    )

    g_effect = coef_summary(
        fit_both,
        2,
    )

    interaction_effect = coef_summary(
        fit_interaction,
        3,
    )

    # --------------------------------------------------------
    # Local discrimination
    # --------------------------------------------------------

    if np.unique(
        y
    ).size == 2:

        auc_z = (
            roc_auc_score(
                y,
                z_raw,
            )
            if z_sd > 1e-12
            else 0.5
        )

        auc_g = roc_auc_score(
            y,
            g,
        )

        auc_both = roc_auc_score(
            y,
            fit_both[
                "prob"
            ],
        )

        auc_interaction = roc_auc_score(
            y,
            fit_interaction[
                "prob"
            ],
        )

    else:
        auc_z = np.nan
        auc_g = np.nan
        auc_both = np.nan
        auc_interaction = np.nan

    # --------------------------------------------------------
    # Risk summaries
    # --------------------------------------------------------

    case_prev = float(
        np.mean(
            y
        )
    )

    carrier_mask = (
        g == 1
    )

    noncarrier_mask = (
        g == 0
    )

    carrier_case_prev = float(
        np.mean(
            y[
                carrier_mask
            ]
        )
    )

    noncarrier_case_prev = float(
        np.mean(
            y[
                noncarrier_mask
            ]
        )
    )

    return {
        "window_center_percentile":
            center,

        "window_lower_percentile":
            lo_pct,

        "window_upper_percentile":
            hi_pct,

        "N":
            n,

        "cases":
            cases,

        "controls":
            controls,

        "case_prevalence":
            case_prev,

        "MUC5B_carrier_N":
            int(
                g.sum()
            ),

        "MUC5B_carrier_prevalence":
            float(
                g.mean()
            ),

        "mean_ZeBRA":
            z_mean,

        "median_ZeBRA":
            float(
                np.median(
                    z_raw
                )
            ),

        "sd_ZeBRA":
            z_sd,

        "unique_ZeBRA_values":
            int(
                np.unique(
                    z_raw
                ).size
            ),

        # Local marginal information
        "LRT_Z_marginal":
            lrt_z_marginal,

        "LRT_Z_marginal_p":
            p_z_marginal,

        "LRT_G_marginal":
            lrt_g_marginal,

        "LRT_G_marginal_p":
            p_g_marginal,

        # Local conditional information
        "LRT_Z_given_G":
            lrt_z_given_g,

        "LRT_Z_given_G_p":
            p_z_given_g,

        "LRT_G_given_Z":
            lrt_g_given_z,

        "LRT_G_given_Z_p":
            p_g_given_z,

        "G_minus_Z_partial_LRT":
            (
                lrt_g_given_z
                - lrt_z_given_g
            ),

        # Interaction
        "LRT_ZxG_interaction":
            lrt_interaction_stat,

        "LRT_ZxG_interaction_p":
            p_interaction,

        # Effect sizes from additive model
        "OR_Z_per_local_SD":
            z_effect[
                "OR"
            ],

        "OR_Z_ci_low":
            z_effect[
                "OR_ci_low"
            ],

        "OR_Z_ci_high":
            z_effect[
                "OR_ci_high"
            ],

        "OR_Z_p":
            z_effect[
                "p"
            ],

        "OR_MUC5B_carrier":
            g_effect[
                "OR"
            ],

        "OR_MUC5B_ci_low":
            g_effect[
                "OR_ci_low"
            ],

        "OR_MUC5B_ci_high":
            g_effect[
                "OR_ci_high"
            ],

        "OR_MUC5B_p":
            g_effect[
                "p"
            ],

        "OR_interaction":
            interaction_effect[
                "OR"
            ],

        "OR_interaction_p":
            interaction_effect[
                "p"
            ],

        # Local AUCs
        "AUC_Z":
            auc_z,

        "AUC_MUC5B":
            auc_g,

        "AUC_Z_plus_MUC5B":
            auc_both,

        "AUC_interaction":
            auc_interaction,

        # Local disease prevalence by genotype
        "case_prevalence_MUC5B_carrier":
            carrier_case_prev,

        "case_prevalence_MUC5B_noncarrier":
            noncarrier_case_prev,

        "carrier_minus_noncarrier_case_prevalence":
            (
                carrier_case_prev
                - noncarrier_case_prev
            ),
    }


# ============================================================
# FPR-defined clinically relevant bands
# ============================================================

def threshold_for_fpr(
    negative_scores,
    fpr,
):
    return float(
        np.quantile(
            negative_scores,
            1.0 - fpr,
            method="higher",
        )
    )


def analyze_operating_bands(
    d,
):
    negative_scores = (
        d.loc[
            d[
                "target"
            ] == 0,
            "predicted_risk",
        ]
        .to_numpy(
            dtype=float
        )
    )

    # Descending-risk intervals.
    fpr_edges = [
        0.0,
        0.005,
        0.01,
        0.02,
        0.05,
        0.10,
        0.20,
        0.40,
    ]

    thresholds = {}

    for fpr in fpr_edges[
        1:
    ]:
        thresholds[
            fpr
        ] = threshold_for_fpr(
            negative_scores,
            fpr,
        )

    rows = []

    # Highest band first: above the 0.5%-FPR threshold.
    band_defs = []

    band_defs.append(
        (
            "above_0.5pct_FPR_threshold",
            thresholds[
                0.005
            ],
            np.inf,
        )
    )

    pairs = [
        (
            0.005,
            0.01,
        ),
        (
            0.01,
            0.02,
        ),
        (
            0.02,
            0.05,
        ),
        (
            0.05,
            0.10,
        ),
        (
            0.10,
            0.20,
        ),
        (
            0.20,
            0.40,
        ),
    ]

    for hi_fpr, lo_fpr in pairs:
        upper_t = thresholds[
            hi_fpr
        ]

        lower_t = thresholds[
            lo_fpr
        ]

        band_defs.append(
            (
                (
                    f"{100*hi_fpr:g}_to_"
                    f"{100*lo_fpr:g}pct_FPR_band"
                ),
                lower_t,
                upper_t,
            )
        )

    for name, lower_t, upper_t in band_defs:

        if np.isinf(
            upper_t
        ):
            sub = d[
                d[
                    "predicted_risk"
                ]
                >= lower_t
            ].copy()
        else:
            sub = d[
                (
                    d[
                        "predicted_risk"
                    ]
                    >= lower_t
                )
                &
                (
                    d[
                        "predicted_risk"
                    ]
                    < upper_t
                )
            ].copy()

        if len(
            sub
        ) < MIN_WINDOW_N:
            continue

        y = sub[
            "target"
        ].astype(int).to_numpy()

        z_raw = sub[
            "predicted_risk"
        ].astype(float).to_numpy()

        g = sub[
            "MUC5B_T_carrier"
        ].astype(int).to_numpy()

        if (
            y.sum() < MIN_CASES
            or (
                len(y)
                - y.sum()
            ) < MIN_CONTROLS
            or np.unique(
                g
            ).size < 2
        ):
            continue

        z_sd = np.std(
            z_raw
        )

        if z_sd > 1e-12:
            z = (
                z_raw
                - z_raw.mean()
            ) / z_sd
        else:
            z = np.zeros_like(
                z_raw
            )

        I = np.ones(
            len(y)
        )

        fit_z = logistic_fit(
            np.column_stack(
                [
                    I,
                    z,
                ]
            ),
            y,
        )

        fit_g = logistic_fit(
            np.column_stack(
                [
                    I,
                    g,
                ]
            ),
            y,
        )

        fit_both = logistic_fit(
            np.column_stack(
                [
                    I,
                    z,
                    g,
                ]
            ),
            y,
        )

        z_info, z_p = lrt(
            fit_both[
                "ll"
            ],
            fit_g[
                "ll"
            ],
        )

        g_info, g_p = lrt(
            fit_both[
                "ll"
            ],
            fit_z[
                "ll"
            ],
        )

        rows.append({
            "band":
                name,

            "lower_ZeBRA":
                float(
                    lower_t
                ),

            "upper_ZeBRA":
                (
                    float(
                        upper_t
                    )
                    if np.isfinite(
                        upper_t
                    )
                    else np.inf
                ),

            "N":
                len(
                    sub
                ),

            "cases":
                int(
                    y.sum()
                ),

            "case_prevalence":
                float(
                    y.mean()
                ),

            "mean_ZeBRA":
                float(
                    z_raw.mean()
                ),

            "LRT_Z_given_G":
                z_info,

            "LRT_Z_given_G_p":
                z_p,

            "LRT_G_given_Z":
                g_info,

            "LRT_G_given_Z_p":
                g_p,

            "G_minus_Z_partial_LRT":
                (
                    g_info
                    - z_info
                ),

            "AUC_Z":
                (
                    roc_auc_score(
                        y,
                        z_raw,
                    )
                    if z_sd > 1e-12
                    else 0.5
                ),

            "AUC_MUC5B":
                roc_auc_score(
                    y,
                    g,
                ),
        })

    return pd.DataFrame(
        rows
    )


# ============================================================
# Load and construct exact notebook-32 cohort
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

# Exact notebook-32 availability filter.
DATA = DATA[
    DATA[
        "predicted_risk"
    ].notnull()
].copy()

DATA = decode_muc5b(
    DATA
)

# Average-rank percentile keeps tied ZeBRA scores together.
DATA[
    "ZeBRA_percentile"
] = (
    DATA[
        "predicted_risk"
    ]
    .rank(
        method="average",
        pct=True,
    )
    * 100.0
)

print()
print("=" * 92)
print("COHORT")
print("=" * 92)

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
    f"ZeBRA global AUC = "
    f"{roc_auc_score(DATA['target'], DATA['predicted_risk']):.4f}"
)


# ============================================================
# Moving-window local-information curves
# ============================================================

WINDOW_ROWS = []

for center in WINDOW_CENTERS:

    result = analyze_window(
        DATA,
        center,
        WINDOW_HALF_WIDTH_PERCENTILE,
    )

    if result is not None:
        WINDOW_ROWS.append(
            result
        )

WINDOW_DF = pd.DataFrame(
    WINDOW_ROWS
)

WINDOW_DF.to_csv(
    OUTDIR
    / "LOCAL_ZEBRA_MUC5B_INFORMATION_WINDOWS.csv",
    index=False,
)


# ============================================================
# Identify candidate crossing and low-information regions
# ============================================================

if len(
    WINDOW_DF
):
    WINDOW_DF[
        "genomics_locally_stronger"
    ] = (
        WINDOW_DF[
            "LRT_G_given_Z"
        ]
        > WINDOW_DF[
            "LRT_Z_given_G"
        ]
    )

    WINDOW_DF[
        "both_weak"
    ] = (
        (
            WINDOW_DF[
                "LRT_G_given_Z"
            ]
            < WEAK_LRT_THRESHOLD
        )
        &
        (
            WINDOW_DF[
                "LRT_Z_given_G"
            ]
            < WEAK_LRT_THRESHOLD
        )
    )

    WINDOW_DF.to_csv(
        OUTDIR
        / "LOCAL_ZEBRA_MUC5B_INFORMATION_WINDOWS.csv",
        index=False,
    )

    CROSSING = WINDOW_DF[
        WINDOW_DF[
            "genomics_locally_stronger"
        ]
    ].copy()

    BOTH_WEAK = WINDOW_DF[
        WINDOW_DF[
            "both_weak"
        ]
    ].copy()

    CROSSING.to_csv(
        OUTDIR
        / "CANDIDATE_GENOMICS_DOMINANT_WINDOWS.csv",
        index=False,
    )

    BOTH_WEAK.to_csv(
        OUTDIR
        / "CANDIDATE_BOTH_WEAK_WINDOWS.csv",
        index=False,
    )


# ============================================================
# Clinically relevant FPR bands
# ============================================================

BAND_DF = analyze_operating_bands(
    DATA
)

BAND_DF.to_csv(
    OUTDIR
    / "FPR_DEFINED_LOCAL_INFORMATION_BANDS.csv",
    index=False,
)


# ============================================================
# Plot 1:
# Partial predictive-information curves
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "LRT_Z_given_G"
    ],
    marker="o",
    label="ZeBRA information given MUC5B",
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "LRT_G_given_Z"
    ],
    marker="o",
    label="MUC5B information given ZeBRA",
)

ax.axhline(
    WEAK_LRT_THRESHOLD,
    linestyle="--",
    linewidth=1,
    label="Weak-information reference",
)

ax.set_xlabel(
    "ZeBRA risk percentile (moving-window center)"
)

ax.set_ylabel(
    "Local partial likelihood-ratio statistic"
)

ax.set_title(
    "Local conditional predictive information: ZeBRA vs MUC5B"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "LOCAL_PARTIAL_INFORMATION_CURVES.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Plot 2:
# Difference curve -- positive means MUC5B contributes more
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 5)
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "G_minus_Z_partial_LRT"
    ],
    marker="o",
)

ax.axhline(
    0.0,
    linestyle="--",
)

ax.set_xlabel(
    "ZeBRA risk percentile (moving-window center)"
)

ax.set_ylabel(
    "MUC5B partial information - ZeBRA partial information"
)

ax.set_title(
    "Local crossover of predictive information"
)

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "LOCAL_INFORMATION_CROSSOVER.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Plot 3:
# Local odds-ratio curves
#
# ZeBRA OR is per 1 local SD within the window.
# MUC5B OR is carrier vs non-carrier.
# These effect sizes are useful but are NOT as directly
# comparable as the 1-df partial-LRT curves above.
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "OR_Z_per_local_SD"
    ],
    marker="o",
    label="ZeBRA OR per local SD",
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "OR_MUC5B_carrier"
    ],
    marker="o",
    label="MUC5B carrier OR",
)

ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_xlabel(
    "ZeBRA risk percentile (moving-window center)"
)

ax.set_ylabel(
    "Local adjusted odds ratio"
)

ax.set_title(
    "Local disease association of ZeBRA and MUC5B"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "LOCAL_ODDS_RATIO_CURVES.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Plot 4:
# Local AUC curves
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "AUC_Z"
    ],
    marker="o",
    label="ZeBRA only",
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "AUC_MUC5B"
    ],
    marker="o",
    label="MUC5B only",
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "AUC_Z_plus_MUC5B"
    ],
    marker="o",
    label="ZeBRA + MUC5B",
)

ax.axhline(
    0.5,
    linestyle="--",
)

ax.set_xlabel(
    "ZeBRA risk percentile (moving-window center)"
)

ax.set_ylabel(
    "Local AUC"
)

ax.set_title(
    "Local discrimination along the ZeBRA risk axis"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "LOCAL_AUC_CURVES.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Plot 5:
# Local case prevalence by MUC5B carrier status
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "case_prevalence_MUC5B_noncarrier"
    ],
    marker="o",
    label="MUC5B non-carrier",
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "case_prevalence_MUC5B_carrier"
    ],
    marker="o",
    label="MUC5B T-carrier",
)

ax.set_xlabel(
    "ZeBRA risk percentile (moving-window center)"
)

ax.set_ylabel(
    "FILD/FILA prevalence"
)

ax.set_title(
    "Disease prevalence by MUC5B status along the ZeBRA risk axis"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "LOCAL_CASE_PREVALENCE_BY_MUC5B.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Plot 6:
# Interaction evidence
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 5)
)

ax.plot(
    WINDOW_DF[
        "window_center_percentile"
    ],
    WINDOW_DF[
        "LRT_ZxG_interaction"
    ],
    marker="o",
)

ax.axhline(
    chi2.ppf(
        0.95,
        1,
    ),
    linestyle="--",
    label="Nominal p=0.05 threshold",
)

ax.set_xlabel(
    "ZeBRA risk percentile (moving-window center)"
)

ax.set_ylabel(
    "1-df LRT for ZeBRA × MUC5B interaction"
)

ax.set_title(
    "Local interaction between ZeBRA and MUC5B"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "LOCAL_INTERACTION_CURVE.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# Plot 7:
# FPR-defined bands
# ============================================================

if len(
    BAND_DF
):

    x = np.arange(
        len(
            BAND_DF
        )
    )

    width = 0.38

    fig, ax = plt.subplots(
        figsize=(11, 6)
    )

    ax.bar(
        x - width / 2,
        BAND_DF[
            "LRT_Z_given_G"
        ],
        width=width,
        label="ZeBRA information given MUC5B",
    )

    ax.bar(
        x + width / 2,
        BAND_DF[
            "LRT_G_given_Z"
        ],
        width=width,
        label="MUC5B information given ZeBRA",
    )

    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        BAND_DF[
            "band"
        ],
        rotation=30,
        ha="right",
    )

    ax.set_ylabel(
        "Partial likelihood-ratio statistic"
    )

    ax.set_title(
        "Conditional predictive information in clinically relevant ZeBRA bands"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        OUTDIR
        / "FPR_DEFINED_INFORMATION_BANDS.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# ============================================================
# Console output
# ============================================================

print()
print("=" * 92)
print("MOVING-WINDOW LOCAL INFORMATION")
print("=" * 92)

print(
    WINDOW_DF[
        [
            "window_center_percentile",
            "N",
            "cases",
            "mean_ZeBRA",
            "LRT_Z_given_G",
            "LRT_Z_given_G_p",
            "LRT_G_given_Z",
            "LRT_G_given_Z_p",
            "G_minus_Z_partial_LRT",
            "AUC_Z",
            "AUC_MUC5B",
            "OR_Z_per_local_SD",
            "OR_MUC5B_carrier",
        ]
    ].to_string(
        index=False
    )
)

print()
print("=" * 92)
print("CANDIDATE LOCAL CROSSING WINDOWS")
print(
    "MUC5B conditional information > ZeBRA conditional information"
)
print("=" * 92)

if len(
    CROSSING
):
    print(
        CROSSING[
            [
                "window_center_percentile",
                "window_lower_percentile",
                "window_upper_percentile",
                "N",
                "cases",
                "mean_ZeBRA",
                "LRT_Z_given_G",
                "LRT_G_given_Z",
                "G_minus_Z_partial_LRT",
                "OR_Z_per_local_SD",
                "OR_MUC5B_carrier",
            ]
        ].to_string(
            index=False
        )
    )
else:
    print(
        "No candidate crossing windows."
    )

print()
print("=" * 92)
print("CANDIDATE BOTH-WEAK WINDOWS")
print("=" * 92)

if len(
    BOTH_WEAK
):
    print(
        BOTH_WEAK[
            [
                "window_center_percentile",
                "N",
                "cases",
                "mean_ZeBRA",
                "LRT_Z_given_G",
                "LRT_G_given_Z",
            ]
        ].to_string(
            index=False
        )
    )
else:
    print(
        "No windows met the current both-weak criterion."
    )

print()
print("=" * 92)
print("FPR-DEFINED CLINICAL BANDS")
print("=" * 92)

if len(
    BAND_DF
):
    print(
        BAND_DF.to_string(
            index=False
        )
    )

print()
print("=" * 92)
print("INTERPRETATION")
print("=" * 92)

print(
    "The primary plot is LOCAL_PARTIAL_INFORMATION_CURVES.png."
)
print(
    "Because both curves are 1-df partial likelihood-ratio statistics, "
    "they are directly comparable."
)
print()
print(
    "The proposed hypothesis predicts:"
)
print(
    "  1. At the extreme high-ZeBRA end: LRT_Z_given_G >> LRT_G_given_Z."
)
print(
    "  2. Immediately below that region: the curves approach or cross, "
    "with LRT_G_given_Z > LRT_Z_given_G in a limited band."
)
print(
    "  3. Farther down the ZeBRA distribution: both statistics become small."
)
print()
print(
    "Important caveat: the x-axis is ZeBRA risk, not directly observed "
    "'true clinical liability'. These curves test the distributional prediction "
    "of that latent-liability hypothesis using ZeBRA as the observable axis."
)

print()
print(
    "Outputs saved to:"
)
print(
    OUTDIR.resolve()
)
