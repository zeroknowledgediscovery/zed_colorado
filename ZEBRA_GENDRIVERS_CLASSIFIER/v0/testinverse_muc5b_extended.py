#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import (
    pearsonr,
    spearmanr,
    mannwhitneyu,
    fisher_exact,
    chi2,
)

from sklearn.preprocessing import (
    StandardScaler,
    SplineTransformer,
)
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_predict,
)
from sklearn.metrics import (
    roc_auc_score,
    confusion_matrix,
    log_loss,
)

try:
    from tqdm.auto import tqdm
except Exception:
    def tqdm(x, **kwargs):
        return x


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

PERCENTILES = [
    50,
    60,
    70,
    75,
    80,
    85,
    90,
    95,
    97.5,
    99,
]

N_PERMUTATIONS = 5000
PERMUTATION_SEED = 1321

N_SPLINE_KNOTS = 5
CV_FOLDS = 5
CV_SEED = 1321


# ============================================================
# Utility functions
# ============================================================

def odds_ratio_ci(a, b, c, d, alpha=0.05):
    """
    2x2 table:

                      Carrier   Non-carrier
        High risk       a          b
        Lower risk      c          d

    Returns OR and approximate Wald 95% CI.
    Uses a 0.5 continuity correction if any cell is zero.
    """
    cells = np.array([a, b, c, d], dtype=float)

    if np.any(cells == 0):
        cells = cells + 0.5

    a2, b2, c2, d2 = cells

    OR = (a2 * d2) / (b2 * c2)

    se_log_or = np.sqrt(
        1 / a2
        + 1 / b2
        + 1 / c2
        + 1 / d2
    )

    zcrit = 1.959963984540054

    lo = np.exp(
        np.log(OR)
        - zcrit * se_log_or
    )

    hi = np.exp(
        np.log(OR)
        + zcrit * se_log_or
    )

    return float(OR), float(lo), float(hi)


def log_likelihood_binary(y, p):
    p = np.clip(
        np.asarray(p, dtype=float),
        1e-12,
        1 - 1e-12,
    )
    y = np.asarray(y, dtype=int)

    return float(
        np.sum(
            y * np.log(p)
            + (1 - y) * np.log(1 - p)
        )
    )


def high_risk_table(
    zebra,
    y,
    percentiles,
):
    baseline_prevalence = float(
        np.mean(y)
    )

    rows = []
    masks = []

    for pct in percentiles:
        threshold = float(
            np.percentile(
                zebra,
                pct,
            )
        )

        high_risk = (
            zebra >= threshold
        )

        pred_binary = high_risk.astype(int)

        tn, fp, fn, tp = confusion_matrix(
            y,
            pred_binary,
            labels=[0, 1],
        ).ravel()

        n_high = int(
            high_risk.sum()
        )

        actual_high_fraction = (
            n_high / len(y)
        )

        carrier_prevalence_high = (
            float(
                y[high_risk].mean()
            )
            if n_high > 0
            else np.nan
        )

        enrichment = (
            carrier_prevalence_high
            / baseline_prevalence
            if baseline_prevalence > 0
            else np.nan
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

        ppv = (
            tp / (tp + fp)
            if (tp + fp) > 0
            else np.nan
        )

        binary_auc = roc_auc_score(
            y,
            pred_binary,
        )

        fisher_or, fisher_p = fisher_exact(
            [
                [tp, fp],
                [fn, tn],
            ]
        )

        OR, OR_lo, OR_hi = odds_ratio_ci(
            tp,
            fp,
            fn,
            tn,
        )

        rows.append({
            "ZeBRA_percentile_cutoff": pct,
            "ZeBRA_threshold": threshold,
            "N_high_risk": n_high,
            "actual_high_fraction": actual_high_fraction,
            "MUC5B_T_carriers_high_risk": int(tp),
            "MUC5B_T_carrier_prevalence_high_risk":
                carrier_prevalence_high,
            "MUC5B_T_carrier_prevalence_all":
                baseline_prevalence,
            "enrichment_fold":
                enrichment,
            "sensitivity":
                sensitivity,
            "specificity":
                specificity,
            "PPV":
                ppv,
            "binary_split_AUC":
                binary_auc,
            "odds_ratio":
                OR,
            "odds_ratio_ci_low":
                OR_lo,
            "odds_ratio_ci_high":
                OR_hi,
            "fisher_exact_odds_ratio":
                fisher_or,
            "fisher_p":
                fisher_p,
        })

        masks.append(
            high_risk.copy()
        )

    return (
        pd.DataFrame(rows),
        masks,
    )


# ============================================================
# Load data
# ============================================================

print("Loading data...")

BASE_DATA = pd.read_csv(
    BASE_DATA_FILE
).drop(
    columns=["target"],
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
# EXACT cohort construction used in notebook 32
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
                    "target"
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

DATA["target"] = (
    pd.to_numeric(
        DATA["target"],
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


# ============================================================
# Inspect MUC5B columns
# ============================================================

MUC5B_COLUMNS = [
    c
    for c in DATA.columns
    if "rs35705950" in str(c)
]

print()
print("MUC5B columns found:")
for c in MUC5B_COLUMNS:
    print(" ", c)

print()
print("MUC5B column frequencies:")

for c in MUC5B_COLUMNS:
    print()
    print(c)
    print(
        DATA[c].value_counts(
            dropna=False
        ).to_string()
    )


# ============================================================
# Validate expected one-hot genotype encoding
#
# Assumed interpretation:
#   G_2 = GG = 0 T-risk alleles
#   G_1 = GT = 1 T-risk allele
#   G_0 = TT = 2 T-risk alleles
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
        + "\n".join(missing)
    )

for c in required:
    DATA[c] = pd.to_numeric(
        DATA[c],
        errors="coerce",
    )

DATA = DATA[
    DATA[
        required
    ].notnull().all(
        axis=1
    )
].copy()

onehot_sum = (
    DATA[
        required
    ].sum(
        axis=1
    )
)

print()
print("MUC5B dummy-state sums BEFORE filtering:")
print(
    onehot_sum.value_counts()
    .sort_index()
    .to_string()
)

# IMPORTANT:
# ILD_TOP_DRIVERS_DATA.csv was created with pd.get_dummies() from a
# categorical genotype coded 0/1/2. If the ORIGINAL genotype was missing,
# pd.get_dummies(dummy_na=False) produces:
#
#     G_0 = 0, G_1 = 0, G_2 = 0
#
# Therefore:
#   sum == 1  -> valid called genotype
#   sum == 0  -> missing/unavailable genotype
#   sum > 1   -> genuinely inconsistent encoding
#
# Do NOT interpret an all-zero row as TT or as a T-risk carrier.

n_missing_genotype = int(
    (onehot_sum == 0).sum()
)

n_valid_genotype = int(
    (onehot_sum == 1).sum()
)

n_invalid_multihot = int(
    (onehot_sum > 1).sum()
)

print()
print("MUC5B genotype availability:")
print(f"  valid called genotype (sum=1): {n_valid_genotype:,}")
print(f"  missing genotype      (sum=0): {n_missing_genotype:,}")
print(f"  invalid multi-hot     (sum>1): {n_invalid_multihot:,}")

if n_invalid_multihot > 0:
    bad = DATA.loc[
        onehot_sum > 1,
        ["patient_id", *required]
    ].head(20)

    print()
    print("Example invalid multi-hot rows:")
    print(
        bad.to_string(
            index=False
        )
    )

    raise ValueError(
        f"{n_invalid_multihot} rs35705950 rows have more than one "
        "active genotype state. These require inspection."
    )

# Retain ONLY patients with an observed/called rs35705950 genotype.
DATA = DATA.loc[
    onehot_sum == 1
].copy()

print()
print(
    f"Proceeding with {len(DATA):,} patients "
    "with a called MUC5B rs35705950 genotype."
)


# ============================================================
# Define genotype endpoints
# ============================================================

muc5b_GG = (
    DATA[
        GG_COLUMN
    ].astype(int)
)

muc5b_GT = (
    DATA[
        GT_COLUMN
    ].astype(int)
)

muc5b_TT = (
    DATA[
        TT_COLUMN
    ].astype(int)
)

# Primary endpoint:
# any T-risk allele = GT or TT
muc5b_T_carrier = (
    (muc5b_GT == 1)
    |
    (muc5b_TT == 1)
).astype(int)

# Dosage:
# GG=0, GT=1, TT=2
muc5b_T_dosage = (
    muc5b_GT
    + 2 * muc5b_TT
).astype(int)

zebra = (
    DATA[
        "predicted_risk"
    ]
    .astype(float)
    .to_numpy()
)

y = (
    muc5b_T_carrier
    .to_numpy(
        dtype=int
    )
)

dosage = (
    muc5b_T_dosage
    .to_numpy(
        dtype=float
    )
)

X = zebra.reshape(
    -1,
    1,
)


# ============================================================
# Cohort summary
# ============================================================

print()
print("=" * 78)
print("MUC5B rs35705950 T-RISK CARRIER ANALYSIS")
print("=" * 78)

print(f"N                    = {len(y):,}")
print(f"GG / non-carrier     = {(y == 0).sum():,}")
print(f"GT or TT / T-carrier = {y.sum():,}")
print(f"T-carrier prevalence = {y.mean():.4%}")

print()
print("Genotype-state counts:")
print(f"  GG = {int(muc5b_GG.sum()):,}")
print(f"  GT = {int(muc5b_GT.sum()):,}")
print(f"  TT = {int(muc5b_TT.sum()):,}")


# ============================================================
# 1. Global association
# ============================================================

pearson_r, pearson_p = pearsonr(
    zebra,
    y,
)

spearman_r, spearman_p = spearmanr(
    zebra,
    y,
)

pearson_dosage_r, pearson_dosage_p = pearsonr(
    zebra,
    dosage,
)

spearman_dosage_r, spearman_dosage_p = spearmanr(
    zebra,
    dosage,
)

z_carrier = zebra[
    y == 1
]

z_noncarrier = zebra[
    y == 0
]

u, mw_p = mannwhitneyu(
    z_carrier,
    z_noncarrier,
    alternative="two-sided",
)

raw_auc = roc_auc_score(
    y,
    zebra,
)

direction_independent_auc = max(
    raw_auc,
    1 - raw_auc,
)

print()
print("Global association")
print("------------------")
print(
    f"Point-biserial r = {pearson_r:.4f}, "
    f"p={pearson_p:.3e}"
)
print(
    f"Spearman rho     = {spearman_r:.4f}, "
    f"p={spearman_p:.3e}"
)
print(
    f"T-dosage Pearson r  = {pearson_dosage_r:.4f}, "
    f"p={pearson_dosage_p:.3e}"
)
print(
    f"T-dosage Spearman rho = {spearman_dosage_r:.4f}, "
    f"p={spearman_dosage_p:.3e}"
)
print(
    f"Mann-Whitney p   = {mw_p:.3e}"
)
print(
    f"Raw continuous AUC = {raw_auc:.4f}"
)
print(
    f"Direction-independent AUC = "
    f"{direction_independent_auc:.4f}"
)


# ============================================================
# 2. Linear logistic model:
#    ONLY predictor = ZeBRA
# ============================================================

cv = StratifiedKFold(
    n_splits=CV_FOLDS,
    shuffle=True,
    random_state=CV_SEED,
)

linear_model = make_pipeline(
    StandardScaler(),
    LogisticRegression(
        C=1e6,
        solver="lbfgs",
        max_iter=5000,
    ),
)

linear_cv_prob = cross_val_predict(
    linear_model,
    X,
    y,
    cv=cv,
    method="predict_proba",
    n_jobs=-1,
)[:, 1]

linear_cv_auc = roc_auc_score(
    y,
    linear_cv_prob,
)

linear_cv_logloss = log_loss(
    y,
    linear_cv_prob,
)

# Effect-size fit using standardized ZeBRA
zebra_sd = (
    (
        zebra
        - zebra.mean()
    )
    / zebra.std()
).reshape(
    -1,
    1,
)

effect_model = LogisticRegression(
    C=1e6,
    solver="lbfgs",
    max_iter=5000,
)

effect_model.fit(
    zebra_sd,
    y,
)

beta = float(
    effect_model.coef_[0, 0]
)

OR_per_SD = float(
    np.exp(beta)
)

print()
print("Linear ZeBRA -> MUC5B carrier model")
print("------------------------------------")
print(
    f"5-fold CV AUC      = {linear_cv_auc:.4f}"
)
print(
    f"5-fold CV log-loss = {linear_cv_logloss:.6f}"
)
print(
    f"OR per 1-SD higher ZeBRA = {OR_per_SD:.3f}"
)


# ============================================================
# 3. NONLINEAR ZeBRA -> carrier model using cubic splines
# ============================================================

spline_model = make_pipeline(
    SplineTransformer(
        n_knots=N_SPLINE_KNOTS,
        degree=3,
        include_bias=False,
    ),
    StandardScaler(),
    LogisticRegression(
        C=1e6,
        solver="lbfgs",
        max_iter=10000,
    ),
)

spline_cv_prob = cross_val_predict(
    spline_model,
    X,
    y,
    cv=cv,
    method="predict_proba",
    n_jobs=-1,
)[:, 1]

spline_cv_auc = roc_auc_score(
    y,
    spline_cv_prob,
)

spline_cv_logloss = log_loss(
    y,
    spline_cv_prob,
)

print()
print("Nonlinear spline ZeBRA -> MUC5B carrier model")
print("----------------------------------------------")
print(
    f"5-fold CV spline AUC      = {spline_cv_auc:.4f}"
)
print(
    f"5-fold CV spline log-loss = {spline_cv_logloss:.6f}"
)
print(
    f"ΔAUC spline - linear      = "
    f"{spline_cv_auc - linear_cv_auc:+.4f}"
)
print(
    f"Δlog-loss spline - linear = "
    f"{spline_cv_logloss - linear_cv_logloss:+.6f}"
)


# ============================================================
# 4. Exploratory likelihood-ratio test:
#    linear vs spline on the complete cohort
#
# Uses very weak regularization (C=1e6) to approximate MLE.
# CV performance above remains the preferred predictive estimate.
# ============================================================

linear_model.fit(
    X,
    y,
)

spline_model.fit(
    X,
    y,
)

p_linear_full = (
    linear_model.predict_proba(
        X
    )[:, 1]
)

p_spline_full = (
    spline_model.predict_proba(
        X
    )[:, 1]
)

ll_linear = log_likelihood_binary(
    y,
    p_linear_full,
)

ll_spline = log_likelihood_binary(
    y,
    p_spline_full,
)

spline_transformer = (
    spline_model.named_steps[
        "splinetransformer"
    ]
)

n_spline_features = int(
    spline_transformer.n_features_out_
)

# Intercept appears in both models.
# Linear model has one ZeBRA term.
lrt_df = max(
    n_spline_features - 1,
    1,
)

lrt_stat = max(
    2 * (
        ll_spline
        - ll_linear
    ),
    0.0,
)

lrt_p = float(
    chi2.sf(
        lrt_stat,
        lrt_df,
    )
)

print()
print("Exploratory nonlinearity test")
print("-----------------------------")
print(
    f"LRT statistic = {lrt_stat:.3f}"
)
print(
    f"df            = {lrt_df}"
)
print(
    f"p             = {lrt_p:.3e}"
)


# ============================================================
# 5. High-ZeBRA-risk enrichment
# ============================================================

THRESHOLD_RESULTS, HIGH_RISK_MASKS = high_risk_table(
    zebra,
    y,
    PERCENTILES,
)

print()
print("=" * 78)
print("HIGH ZeBRA RISK -> MUC5B T-CARRIER ENRICHMENT")
print("=" * 78)

print(
    THRESHOLD_RESULTS.to_string(
        index=False
    )
)

THRESHOLD_RESULTS.to_csv(
    "ZEBRA_HIGH_RISK_MUC5B_T_CARRIER_ENRICHMENT.csv",
    index=False,
)


# ============================================================
# 6. Specific top-1% / score-saturation diagnostics
# ============================================================

top99 = THRESHOLD_RESULTS.loc[
    THRESHOLD_RESULTS[
        "ZeBRA_percentile_cutoff"
    ] == 99
].iloc[0]

print()
print("Top-1% threshold diagnostic")
print("---------------------------")
print(
    f"Requested percentile = 99"
)
print(
    f"Score threshold       = "
    f"{top99['ZeBRA_threshold']:.8f}"
)
print(
    f"Actual N above cutoff = "
    f"{int(top99['N_high_risk']):,}"
)
print(
    f"Actual fraction       = "
    f"{top99['actual_high_fraction']:.4%}"
)
print(
    f"Carrier prevalence    = "
    f"{top99['MUC5B_T_carrier_prevalence_high_risk']:.4%}"
)
print(
    f"Fold enrichment       = "
    f"{top99['enrichment_fold']:.3f}x"
)
print(
    f"OR                    = "
    f"{top99['odds_ratio']:.3f} "
    f"({top99['odds_ratio_ci_low']:.3f}, "
    f"{top99['odds_ratio_ci_high']:.3f})"
)
print(
    f"Fisher p              = "
    f"{top99['fisher_p']:.3e}"
)

# Explicit score saturation at exactly 1.0
sat_mask = np.isclose(
    zebra,
    1.0,
)

if sat_mask.any():
    sat_n = int(
        sat_mask.sum()
    )

    sat_carrier_prev = float(
        y[
            sat_mask
        ].mean()
    )

    baseline_prev = float(
        y.mean()
    )

    sat_enrichment = (
        sat_carrier_prev
        / baseline_prev
    )

    a = int(
        y[
            sat_mask
        ].sum()
    )
    b = int(
        sat_n - a
    )
    c = int(
        y[
            ~sat_mask
        ].sum()
    )
    d = int(
        (~sat_mask).sum()
        - c
    )

    sat_fisher_or, sat_fisher_p = (
        fisher_exact(
            [
                [a, b],
                [c, d],
            ]
        )
    )

    sat_or, sat_lo, sat_hi = (
        odds_ratio_ci(
            a,
            b,
            c,
            d,
        )
    )

    print()
    print("Exact ZeBRA score == 1.0 diagnostic")
    print("------------------------------------")
    print(
        f"N score==1.0       = {sat_n:,}"
    )
    print(
        f"Carrier prevalence = "
        f"{sat_carrier_prev:.4%}"
    )
    print(
        f"Fold enrichment    = "
        f"{sat_enrichment:.3f}x"
    )
    print(
        f"OR                 = "
        f"{sat_or:.3f} "
        f"({sat_lo:.3f}, {sat_hi:.3f})"
    )
    print(
        f"Fisher p           = "
        f"{sat_fisher_p:.3e}"
    )
else:
    sat_n = 0
    sat_carrier_prev = np.nan
    sat_enrichment = np.nan
    sat_or = np.nan
    sat_lo = np.nan
    sat_hi = np.nan
    sat_fisher_p = np.nan


# ============================================================
# 7. Permutation test for MAXIMUM enrichment over all
#    candidate ZeBRA percentile thresholds
#
# This accounts for choosing the "best-looking" threshold
# after examining the enrichment curve.
# ============================================================

observed_enrichments = (
    THRESHOLD_RESULTS[
        "enrichment_fold"
    ].to_numpy(
        dtype=float
    )
)

observed_max_enrichment = float(
    np.nanmax(
        observed_enrichments
    )
)

observed_max_index = int(
    np.nanargmax(
        observed_enrichments
    )
)

observed_best_percentile = (
    PERCENTILES[
        observed_max_index
    ]
)

rng = np.random.default_rng(
    PERMUTATION_SEED
)

baseline_prevalence = float(
    y.mean()
)

perm_max_enrichment = np.empty(
    N_PERMUTATIONS,
    dtype=float,
)

# Also retain the fixed 99th-percentile enrichment
# under permutation.
pct99_index = PERCENTILES.index(
    99
)

perm_99_enrichment = np.empty(
    N_PERMUTATIONS,
    dtype=float,
)

print()
print(
    f"Running {N_PERMUTATIONS:,} "
    "MUC5B-label permutations..."
)

for b in tqdm(
    range(
        N_PERMUTATIONS
    ),
    desc="Permutation test",
):

    yp = rng.permutation(
        y
    )

    enrichment_values = []

    for mask in HIGH_RISK_MASKS:
        if mask.sum() == 0:
            enrichment_values.append(
                np.nan
            )
            continue

        high_prev = float(
            yp[
                mask
            ].mean()
        )

        enrichment_values.append(
            high_prev
            / baseline_prevalence
        )

    enrichment_values = np.asarray(
        enrichment_values,
        dtype=float,
    )

    perm_max_enrichment[b] = (
        np.nanmax(
            enrichment_values
        )
    )

    perm_99_enrichment[b] = (
        enrichment_values[
            pct99_index
        ]
    )


max_enrichment_permutation_p = (
    1
    + np.sum(
        perm_max_enrichment
        >= observed_max_enrichment
    )
) / (
    N_PERMUTATIONS
    + 1
)

observed_99_enrichment = float(
    THRESHOLD_RESULTS.loc[
        THRESHOLD_RESULTS[
            "ZeBRA_percentile_cutoff"
        ] == 99,
        "enrichment_fold",
    ].iloc[0]
)

fixed99_permutation_p = (
    1
    + np.sum(
        perm_99_enrichment
        >= observed_99_enrichment
    )
) / (
    N_PERMUTATIONS
    + 1
)


PERMUTATION_SUMMARY = pd.DataFrame([
    {
        "N_permutations":
            N_PERMUTATIONS,

        "observed_max_enrichment":
            observed_max_enrichment,

        "observed_best_percentile":
            observed_best_percentile,

        "permutation_p_max_enrichment":
            max_enrichment_permutation_p,

        "null_max_enrichment_mean":
            float(
                perm_max_enrichment.mean()
            ),

        "null_max_enrichment_q95":
            float(
                np.quantile(
                    perm_max_enrichment,
                    0.95,
                )
            ),

        "null_max_enrichment_q99":
            float(
                np.quantile(
                    perm_max_enrichment,
                    0.99,
                )
            ),

        "observed_99_enrichment":
            observed_99_enrichment,

        "permutation_p_fixed_99":
            fixed99_permutation_p,
    }
])

PERMUTATION_SUMMARY.to_csv(
    "ZEBRA_MUC5B_ENRICHMENT_PERMUTATION_TEST.csv",
    index=False,
)

print()
print("Permutation results")
print("-------------------")
print(
    PERMUTATION_SUMMARY.to_string(
        index=False
    )
)


# ============================================================
# 8. Main analysis summary
# ============================================================

SUMMARY = pd.DataFrame([
    {
        "N":
            len(y),

        "MUC5B_genotype_missing_excluded_N":
            n_missing_genotype,

        "MUC5B_genotype_valid_N":
            n_valid_genotype,

        "MUC5B_GG_N":
            int(
                muc5b_GG.sum()
            ),

        "MUC5B_GT_N":
            int(
                muc5b_GT.sum()
            ),

        "MUC5B_TT_N":
            int(
                muc5b_TT.sum()
            ),

        "MUC5B_T_carrier_N":
            int(
                y.sum()
            ),

        "MUC5B_T_carrier_prevalence":
            float(
                y.mean()
            ),

        "point_biserial_r":
            pearson_r,

        "point_biserial_p":
            pearson_p,

        "spearman_rho":
            spearman_r,

        "spearman_p":
            spearman_p,

        "T_dosage_pearson_r":
            pearson_dosage_r,

        "T_dosage_pearson_p":
            pearson_dosage_p,

        "T_dosage_spearman_rho":
            spearman_dosage_r,

        "T_dosage_spearman_p":
            spearman_dosage_p,

        "mannwhitney_p":
            mw_p,

        "raw_auc":
            raw_auc,

        "direction_independent_auc":
            direction_independent_auc,

        "linear_cv_auc":
            linear_cv_auc,

        "linear_cv_logloss":
            linear_cv_logloss,

        "spline_cv_auc":
            spline_cv_auc,

        "spline_cv_logloss":
            spline_cv_logloss,

        "spline_minus_linear_auc":
            spline_cv_auc
            - linear_cv_auc,

        "spline_minus_linear_logloss":
            spline_cv_logloss
            - linear_cv_logloss,

        "nonlinearity_LRT_stat":
            lrt_stat,

        "nonlinearity_LRT_df":
            lrt_df,

        "nonlinearity_LRT_p":
            lrt_p,

        "OR_per_SD_zebra":
            OR_per_SD,

        "observed_max_tail_enrichment":
            observed_max_enrichment,

        "best_tail_percentile":
            observed_best_percentile,

        "permutation_p_max_tail_enrichment":
            max_enrichment_permutation_p,

        "top99_enrichment":
            observed_99_enrichment,

        "top99_permutation_p":
            fixed99_permutation_p,

        "score1_N":
            sat_n,

        "score1_enrichment":
            sat_enrichment,

        "score1_OR":
            sat_or,

        "score1_OR_ci_low":
            sat_lo,

        "score1_OR_ci_high":
            sat_hi,

        "score1_fisher_p":
            sat_fisher_p,
    }
])

SUMMARY.to_csv(
    "ZEBRA_MUC5B_T_CARRIER_ASSOCIATION_EXTENDED.csv",
    index=False,
)

print()
print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)
print(
    SUMMARY.to_string(
        index=False
    )
)


# ============================================================
# 9. Plot: carrier status
# ============================================================

fig, ax = plt.subplots(
    figsize=(7, 5)
)

ax.boxplot(
    [
        zebra[
            y == 0
        ],
        zebra[
            y == 1
        ],
    ],
    labels=[
        "GG / no T risk allele",
        "GT or TT / T carrier",
    ],
    showmeans=True,
)

ax.set_ylabel(
    "ZeBRA predicted_risk"
)

ax.set_title(
    "ZeBRA score by MUC5B rs35705950 T-risk carrier status"
)

plt.tight_layout()

plt.savefig(
    "ZEBRA_BY_MUC5B_T_CARRIER.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()


# ============================================================
# 10. Plot: T-allele dosage
# ============================================================

fig, ax = plt.subplots(
    figsize=(7, 5)
)

ax.boxplot(
    [
        zebra[
            dosage == 0
        ],
        zebra[
            dosage == 1
        ],
        zebra[
            dosage == 2
        ],
    ],
    labels=[
        "GG (0 T)",
        "GT (1 T)",
        "TT (2 T)",
    ],
    showmeans=True,
)

ax.set_ylabel(
    "ZeBRA predicted_risk"
)

ax.set_title(
    "ZeBRA score by MUC5B rs35705950 T-allele dosage"
)

plt.tight_layout()

plt.savefig(
    "ZEBRA_BY_MUC5B_T_DOSAGE.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()


# ============================================================
# 11. Plot: tail enrichment
# ============================================================

fig, ax = plt.subplots(
    figsize=(7, 5)
)

ax.plot(
    THRESHOLD_RESULTS[
        "ZeBRA_percentile_cutoff"
    ],
    THRESHOLD_RESULTS[
        "enrichment_fold"
    ],
    marker="o",
)

ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_xlabel(
    "ZeBRA percentile threshold"
)

ax.set_ylabel(
    "Fold enrichment of MUC5B T carriers"
)

ax.set_title(
    "MUC5B T-risk carrier enrichment at high ZeBRA risk"
)

plt.tight_layout()

plt.savefig(
    "ZEBRA_HIGH_RISK_MUC5B_T_CARRIER_ENRICHMENT.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()


# ============================================================
# 12. Plot: nonlinear relationship
#
# Shows empirical carrier rate across ZeBRA quantile bins plus
# the fitted spline probability curve.
# ============================================================

plot_df = pd.DataFrame({
    "zebra": zebra,
    "carrier": y,
})

plot_df["bin"] = pd.qcut(
    plot_df[
        "zebra"
    ].rank(
        method="first"
    ),
    q=20,
    labels=False,
    duplicates="drop",
)

empirical = (
    plot_df
    .groupby(
        "bin"
    )
    .agg(
        mean_zebra=(
            "zebra",
            "mean",
        ),
        carrier_rate=(
            "carrier",
            "mean",
        ),
        n=(
            "carrier",
            "size",
        ),
    )
    .reset_index()
)

grid = np.linspace(
    float(
        zebra.min()
    ),
    float(
        zebra.max()
    ),
    500,
).reshape(
    -1,
    1,
)

spline_prob_grid = (
    spline_model.predict_proba(
        grid
    )[:, 1]
)

linear_prob_grid = (
    linear_model.predict_proba(
        grid
    )[:, 1]
)

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.scatter(
    empirical[
        "mean_zebra"
    ],
    empirical[
        "carrier_rate"
    ],
    s=35,
    label="Observed carrier rate (20 quantile bins)",
)

ax.plot(
    grid.ravel(),
    linear_prob_grid,
    linewidth=2,
    label="Linear logistic",
)

ax.plot(
    grid.ravel(),
    spline_prob_grid,
    linewidth=2,
    label="Spline logistic",
)

ax.set_xlabel(
    "ZeBRA predicted_risk"
)

ax.set_ylabel(
    "P(MUC5B T-risk carrier)"
)

ax.set_title(
    "Nonlinear association between ZeBRA risk and MUC5B carrier status"
)

ax.legend()

plt.tight_layout()

plt.savefig(
    "ZEBRA_MUC5B_NONLINEAR_ASSOCIATION.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()


# ============================================================
# 13. Plot: permutation null distribution of max enrichment
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.hist(
    perm_max_enrichment,
    bins=40,
)

ax.axvline(
    observed_max_enrichment,
    linestyle="--",
    linewidth=2,
    label=(
        f"Observed max = "
        f"{observed_max_enrichment:.3f}x"
    ),
)

ax.set_xlabel(
    "Maximum enrichment across tested ZeBRA thresholds"
)

ax.set_ylabel(
    "Permutation count"
)

ax.set_title(
    "Permutation null: maximum MUC5B enrichment over ZeBRA thresholds"
)

ax.legend()

plt.tight_layout()

plt.savefig(
    "ZEBRA_MUC5B_MAX_ENRICHMENT_PERMUTATION.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()


print()
print("Analysis complete.")
print("Saved outputs:")
for filename in [
    "ZEBRA_MUC5B_T_CARRIER_ASSOCIATION_EXTENDED.csv",
    "ZEBRA_HIGH_RISK_MUC5B_T_CARRIER_ENRICHMENT.csv",
    "ZEBRA_MUC5B_ENRICHMENT_PERMUTATION_TEST.csv",
    "ZEBRA_BY_MUC5B_T_CARRIER.png",
    "ZEBRA_BY_MUC5B_T_DOSAGE.png",
    "ZEBRA_HIGH_RISK_MUC5B_T_CARRIER_ENRICHMENT.png",
    "ZEBRA_MUC5B_NONLINEAR_ASSOCIATION.png",
    "ZEBRA_MUC5B_MAX_ENRICHMENT_PERMUTATION.png",
]:
    print(" ", filename)
