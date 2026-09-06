#!/usr/bin/env python3

"""
Compare:
    1. ZeBRA alone
    2. Fixed MUC5B-switch hybrid
    3. Upper convex hull of the UNION of both ROC curves

The hybrid is generated out-of-fold so that every subject receives a
prediction from a model fitted without that subject.

The two empirical ROC dataframes are then concatenated:

    roc_union = pd.concat([roc_zebra, roc_hybrid])

and the union is passed to zedstat.processRoc followed by:

    zt_union.smooth(..., convexify=True)

zedstat therefore computes the upper ROC convex hull across BOTH
classifiers.

Important:
The hull is not a third ordinary scalar classifier.  It is the best
operating envelope achievable by selecting between the two classifiers
at different operating regions (and, along hull line segments, by
randomization/mixture in the usual ROC-convex-hull interpretation).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_curve
from sklearn.model_selection import StratifiedKFold

from zedstat import zedstat


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

# Fixed regions from the previous exploratory analysis.
MUC5B_SWITCH_WINDOWS = [
    (0.0, 50.0, "P00_P50"),
    (65.0, 75.0, "P65_P75"),
    (85.0, 95.0, "P85_P95"),
]

N_FOLDS = 5
RANDOM_STATE = 1321

# Jeffreys smoothing for P(Y | region, genotype).
ALPHA = 0.5
BETA = 0.5

# zedstat ROC interpolation grid.
ZEDSTAT_STEP = 0.0005

# processRoc settings.
ZEDSTAT_ORDER = 3
ZEDSTAT_ALPHA = 0.05

OUTDIR = Path(
    "./RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL"
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
            "Some rs35705950 rows have >1 active genotype state."
        )

    # All-zero corresponds to missing original genotype.
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

    return d


def empirical_percentile(
    reference_scores,
    query_scores,
):
    """
    Percentile of query scores relative to TRAINING scores only.
    """
    ref = np.sort(
        np.asarray(
            reference_scores,
            dtype=float,
        )
    )

    query = np.asarray(
        query_scores,
        dtype=float,
    )

    rank = np.searchsorted(
        ref,
        query,
        side="right",
    )

    return (
        100.0
        * rank
        / len(ref)
    )


def fit_zebra_calibrator(
    z_train,
    y_train,
):
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


def estimate_region_genotype_risks(
    y_train,
    train_percentile,
    g_train,
):
    """
    Smoothed P(Y=1 | switch-region, MUC5B carrier status)
    estimated from TRAINING data only.
    """
    lookup = {}

    rows = []

    for lo, hi, name in MUC5B_SWITCH_WINDOWS:

        in_region = (
            (train_percentile >= lo)
            &
            (train_percentile < hi)
        )

        for g in [
            0,
            1,
        ]:
            mask = (
                in_region
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

            probability = (
                cases + ALPHA
            ) / (
                n
                + ALPHA
                + BETA
            )

            lookup[
                (
                    name,
                    g,
                )
            ] = float(
                probability
            )

            rows.append({
                "region":
                    name,

                "MUC5B_carrier":
                    g,

                "N":
                    n,

                "cases":
                    cases,

                "probability":
                    probability,
            })

    return (
        lookup,
        pd.DataFrame(
            rows
        ),
    )


def region_for_percentile(
    percentile,
):
    p = float(
        percentile
    )

    for lo, hi, name in MUC5B_SWITCH_WINDOWS:
        if (
            p >= lo
            and p < hi
        ):
            return name

    return None


def make_hybrid_score(
    z,
    percentiles,
    g,
    zebra_calibrator,
    region_lookup,
):
    """
    Outside selected regions:
        isotonic-calibrated ZeBRA probability.

    Inside selected regions:
        replace exact ZeBRA ranking by
        P(Y | region, MUC5B carrier status).
    """
    z = np.asarray(
        z,
        dtype=float,
    )

    percentiles = np.asarray(
        percentiles,
        dtype=float,
    )

    g = np.asarray(
        g,
        dtype=int,
    )

    score = np.asarray(
        zebra_calibrator.predict(
            z
        ),
        dtype=float,
    )

    used_region = np.full(
        len(score),
        "",
        dtype=object,
    )

    for i in range(
        len(score)
    ):
        region = region_for_percentile(
            percentiles[i]
        )

        if region is not None:
            score[i] = (
                region_lookup[
                    (
                        region,
                        int(
                            g[i]
                        ),
                    )
                ]
            )

            used_region[i] = region

    return (
        score,
        used_region,
    )


def make_roc_df(
    y,
    score,
    model_name,
):
    fpr, tpr, threshold = roc_curve(
        y,
        score,
        pos_label=1,
        drop_intermediate=False,
    )

    return pd.DataFrame({
        "fpr":
            fpr,

        "tpr":
            tpr,

        "threshold":
            threshold,

        "model":
            model_name,
    })


def process_with_zedstat(
    roc_df,
    total_samples,
    positive_samples,
    prevalence,
    convexify,
):
    """
    zedstat AUC on an empirical ROC dataframe.

    For the union ROC we deliberately omit model-specific thresholds,
    because a point can originate from either classifier.
    """
    zdf = roc_df[
        [
            "fpr",
            "tpr",
        ]
    ].copy()

    zt = zedstat.processRoc(
        df=zdf,
        order=ZEDSTAT_ORDER,
        total_samples=total_samples,
        positive_samples=positive_samples,
        alpha=ZEDSTAT_ALPHA,
        prevalence=prevalence,
    )

    zt.smooth(
        STEP=ZEDSTAT_STEP,
        interpolate=True,
        convexify=convexify,
    )

    # Current zedstat.auc() returns a tuple:
    #     (nominal_auc, lower_bound, upper_bound)
    # rather than a scalar.  We want the nominal AUC of the
    # currently processed ROC curve here.
    auc_result = zt.auc()

    if isinstance(
        auc_result,
        (tuple, list, np.ndarray),
    ):
        auc_value = float(
            auc_result[0]
        )
    else:
        auc_value = float(
            auc_result
        )

    processed = (
        zt.get()
        .reset_index()
    )

    # Current zedstat normally returns fpr as the index.
    if (
        "fpr" not in processed.columns
        and processed.columns[0] != "fpr"
    ):
        processed = processed.rename(
            columns={
                processed.columns[0]:
                    "fpr"
            }
        )

    return (
        zt,
        processed,
        auc_value,
    )


def upper_envelope_at_fpr(
    roc_a,
    roc_b,
    grid,
):
    """
    Pointwise max of the two interpolated ROC curves.

    This is NOT necessarily the convex hull; it is included only
    as a diagnostic comparison.
    """
    def prep(df):
        d = (
            df[
                [
                    "fpr",
                    "tpr",
                ]
            ]
            .sort_values(
                [
                    "fpr",
                    "tpr",
                ],
                ascending=[
                    True,
                    False,
                ],
            )
            .groupby(
                "fpr",
                as_index=False,
            )
            .first()
            .sort_values(
                "fpr"
            )
        )

        return (
            d["fpr"].to_numpy(
                dtype=float
            ),
            d["tpr"].to_numpy(
                dtype=float
            ),
        )

    xa, ya = prep(
        roc_a
    )

    xb, yb = prep(
        roc_b
    )

    ta = np.interp(
        grid,
        xa,
        ya,
    )

    tb = np.interp(
        grid,
        xb,
        yb,
    )

    return np.maximum(
        ta,
        tb,
    )


# ============================================================
# LOAD AND BUILD COHORT
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

DATA = DATA.reset_index(
    drop=True
)

print()
print("=" * 96)
print("COHORT")
print("=" * 96)

print(
    f"N = {len(DATA):,}"
)

print(
    f"Cases = "
    f"{int(DATA['target'].sum()):,}"
)

print(
    f"Controls = "
    f"{int((DATA['target'] == 0).sum()):,}"
)

print(
    f"MUC5B carriers = "
    f"{int(DATA['MUC5B_T_carrier'].sum()):,}"
)


# ============================================================
# 5-FOLD CROSS-FITTED HYBRID PREDICTIONS
# ============================================================

y = DATA[
    "target"
].to_numpy(
    dtype=int
)

z = DATA[
    "predicted_risk"
].to_numpy(
    dtype=float
)

g = DATA[
    "MUC5B_T_carrier"
].to_numpy(
    dtype=int
)

oof_hybrid = np.full(
    len(DATA),
    np.nan,
    dtype=float,
)

oof_zebra_calibrated = np.full(
    len(DATA),
    np.nan,
    dtype=float,
)

oof_region = np.full(
    len(DATA),
    "",
    dtype=object,
)

fold_id = np.full(
    len(DATA),
    -1,
    dtype=int,
)

RISK_TABLES = []

cv = StratifiedKFold(
    n_splits=N_FOLDS,
    shuffle=True,
    random_state=RANDOM_STATE,
)

print()
print("=" * 96)
print(
    f"GENERATING {N_FOLDS}-FOLD OUT-OF-FOLD HYBRID SCORES"
)
print("=" * 96)

for fold, (
    train_idx,
    test_idx,
) in enumerate(
    cv.split(
        z.reshape(
            -1,
            1,
        ),
        y,
    ),
    start=1,
):

    print(
        f"fold {fold}/{N_FOLDS}: "
        f"train={len(train_idx):,}, "
        f"test={len(test_idx):,}"
    )

    y_train = y[
        train_idx
    ]

    z_train = z[
        train_idx
    ]

    g_train = g[
        train_idx
    ]

    z_test = z[
        test_idx
    ]

    g_test = g[
        test_idx
    ]

    train_pct = empirical_percentile(
        z_train,
        z_train,
    )

    test_pct = empirical_percentile(
        z_train,
        z_test,
    )

    calibrator = fit_zebra_calibrator(
        z_train,
        y_train,
    )

    (
        region_lookup,
        risk_table,
    ) = estimate_region_genotype_risks(
        y_train,
        train_pct,
        g_train,
    )

    risk_table.insert(
        0,
        "fold",
        fold,
    )

    RISK_TABLES.append(
        risk_table
    )

    hybrid_test, region_test = make_hybrid_score(
        z_test,
        test_pct,
        g_test,
        calibrator,
        region_lookup,
    )

    oof_hybrid[
        test_idx
    ] = hybrid_test

    oof_zebra_calibrated[
        test_idx
    ] = calibrator.predict(
        z_test
    )

    oof_region[
        test_idx
    ] = region_test

    fold_id[
        test_idx
    ] = fold


if np.isnan(
    oof_hybrid
).any():
    raise RuntimeError(
        "Some subjects did not receive an OOF hybrid prediction."
    )


OOF = DATA[
    [
        "patient_id",
        "target",
        "predicted_risk",
        "MUC5B_T_carrier",
    ]
].copy()

OOF[
    "fold"
] = fold_id

OOF[
    "zebra_calibrated_oof"
] = oof_zebra_calibrated

OOF[
    "hybrid_oof"
] = oof_hybrid

OOF[
    "switch_region"
] = oof_region

OOF.to_csv(
    OUTDIR
    / "OOF_ZEBRA_HYBRID_SCORES.csv",
    index=False,
)

pd.concat(
    RISK_TABLES,
    ignore_index=True,
).to_csv(
    OUTDIR
    / "OOF_TRAINING_REGION_RISK_TABLES.csv",
    index=False,
)


# ============================================================
# GENERATE THE TWO ROC CURVES
# ============================================================

roc_zebra = make_roc_df(
    y,
    z,
    "ZeBRA",
)

roc_hybrid = make_roc_df(
    y,
    oof_hybrid,
    "Fixed_MUC5B_switch_hybrid",
)

roc_zebra.to_csv(
    OUTDIR
    / "ROC_ZEBRA.csv",
    index=False,
)

roc_hybrid.to_csv(
    OUTDIR
    / "ROC_HYBRID.csv",
    index=False,
)


# ============================================================
# INDIVIDUAL zedstat AUCs
#
# convexify=False:
# preserve each empirical ROC shape while interpolating.
# ============================================================

n_total = len(
    y
)

n_positive = int(
    y.sum()
)

prevalence = float(
    y.mean()
)

(
    zt_zebra,
    zed_zebra_df,
    auc_zebra,
) = process_with_zedstat(
    roc_zebra,
    total_samples=n_total,
    positive_samples=n_positive,
    prevalence=prevalence,
    convexify=False,
)

(
    zt_hybrid,
    zed_hybrid_df,
    auc_hybrid,
) = process_with_zedstat(
    roc_hybrid,
    total_samples=n_total,
    positive_samples=n_positive,
    prevalence=prevalence,
    convexify=False,
)


# ============================================================
# CONCATENATE THE TWO ROC DATAFRAMES
#
# This is the central requested operation.
# ============================================================

roc_union = pd.concat(
    [
        roc_zebra[
            [
                "fpr",
                "tpr",
                "model",
            ]
        ],
        roc_hybrid[
            [
                "fpr",
                "tpr",
                "model",
            ]
        ],
    ],
    ignore_index=True,
)

roc_union.to_csv(
    OUTDIR
    / "ROC_UNION_RAW_POINTS.csv",
    index=False,
)


# ============================================================
# zedstat UPPER CONVEX HULL OF THE UNION
#
# processRoc only needs fpr/tpr here.
#
# zedstat's preparation step retains the highest TPR if the two
# curves contain the same FPR, then smooth(convexify=True)
# constructs the upper ROC hull.
# ============================================================

(
    zt_hull,
    hull_df,
    auc_hull,
) = process_with_zedstat(
    roc_union,
    total_samples=n_total,
    positive_samples=n_positive,
    prevalence=prevalence,
    convexify=True,
)

hull_df.to_csv(
    OUTDIR
    / "ROC_UNION_ZEDSTAT_CONVEX_HULL.csv",
    index=False,
)


# ============================================================
# ALSO COMPUTE INDIVIDUAL CONVEXIFIED AUCs
#
# Useful diagnostic: union-hull improvement should not simply be
# caused by convexification of one curve by itself.
# ============================================================

(
    _,
    zebra_hull_df,
    auc_zebra_hull,
) = process_with_zedstat(
    roc_zebra,
    total_samples=n_total,
    positive_samples=n_positive,
    prevalence=prevalence,
    convexify=True,
)

(
    _,
    hybrid_hull_df,
    auc_hybrid_hull,
) = process_with_zedstat(
    roc_hybrid,
    total_samples=n_total,
    positive_samples=n_positive,
    prevalence=prevalence,
    convexify=True,
)


# ============================================================
# POINTWISE UPPER ENVELOPE DIAGNOSTIC
#
# Different from convex hull; included so we can see how much
# additional area comes specifically from hull convexification.
# ============================================================

grid = np.arange(
    0.0,
    1.0 + ZEDSTAT_STEP / 2.0,
    ZEDSTAT_STEP,
)

envelope_tpr = upper_envelope_at_fpr(
    roc_zebra,
    roc_hybrid,
    grid,
)

envelope_df = pd.DataFrame({
    "fpr":
        grid,

    "tpr":
        envelope_tpr,
})

envelope_df.to_csv(
    OUTDIR
    / "ROC_POINTWISE_UPPER_ENVELOPE.csv",
    index=False,
)

# Use zedstat here too, but WITHOUT convexifying.
(
    _,
    envelope_processed,
    auc_envelope,
) = process_with_zedstat(
    envelope_df,
    total_samples=n_total,
    positive_samples=n_positive,
    prevalence=prevalence,
    convexify=False,
)


# ============================================================
# SUMMARY
# ============================================================

SUMMARY = pd.DataFrame([
    {
        "curve":
            "ZeBRA empirical/interpolated",

        "zedstat_auc":
            auc_zebra,
    },
    {
        "curve":
            "Hybrid empirical/interpolated",

        "zedstat_auc":
            auc_hybrid,
    },
    {
        "curve":
            "ZeBRA own convex hull",

        "zedstat_auc":
            auc_zebra_hull,
    },
    {
        "curve":
            "Hybrid own convex hull",

        "zedstat_auc":
            auc_hybrid_hull,
    },
    {
        "curve":
            "Pointwise max of ZeBRA and hybrid",

        "zedstat_auc":
            auc_envelope,
    },
    {
        "curve":
            "Convex hull of concatenated ZeBRA + hybrid ROC points",

        "zedstat_auc":
            auc_hull,
    },
])

SUMMARY[
    "gain_vs_ZeBRA_empirical"
] = (
    SUMMARY[
        "zedstat_auc"
    ]
    - auc_zebra
)

SUMMARY[
    "gain_vs_best_individual_empirical"
] = (
    SUMMARY[
        "zedstat_auc"
    ]
    - max(
        auc_zebra,
        auc_hybrid,
    )
)

SUMMARY.to_csv(
    OUTDIR
    / "ZEDSTAT_ROC_HULL_AUC_SUMMARY.csv",
    index=False,
)


# ============================================================
# FIND WHICH ORIGINAL CURVE IS BETTER AT EACH FPR
#
# This helps identify the switching regions along the ROC.
# ============================================================

def interpolate_roc(
    roc_df,
    grid,
):
    d = (
        roc_df[
            [
                "fpr",
                "tpr",
            ]
        ]
        .sort_values(
            [
                "fpr",
                "tpr",
            ],
            ascending=[
                True,
                False,
            ],
        )
        .groupby(
            "fpr",
            as_index=False,
        )
        .first()
        .sort_values(
            "fpr"
        )
    )

    return np.interp(
        grid,
        d[
            "fpr"
        ].to_numpy(
            dtype=float
        ),
        d[
            "tpr"
        ].to_numpy(
            dtype=float
        ),
    )


zebra_grid_tpr = interpolate_roc(
    roc_zebra,
    grid,
)

hybrid_grid_tpr = interpolate_roc(
    roc_hybrid,
    grid,
)

winner = np.where(
    hybrid_grid_tpr > zebra_grid_tpr,
    "Hybrid",
    np.where(
        zebra_grid_tpr > hybrid_grid_tpr,
        "ZeBRA",
        "Tie",
    ),
)

SWITCH_GRID = pd.DataFrame({
    "fpr":
        grid,

    "zebra_tpr":
        zebra_grid_tpr,

    "hybrid_tpr":
        hybrid_grid_tpr,

    "winner":
        winner,
})

SWITCH_GRID.to_csv(
    OUTDIR
    / "ROC_MODEL_WINNER_BY_FPR.csv",
    index=False,
)


# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 8)
)

ax.plot(
    roc_zebra[
        "fpr"
    ],
    roc_zebra[
        "tpr"
    ],
    linewidth=2,
    label=(
        f"ZeBRA "
        f"(AUC={auc_zebra:.4f})"
    ),
)

ax.plot(
    roc_hybrid[
        "fpr"
    ],
    roc_hybrid[
        "tpr"
    ],
    linewidth=2,
    label=(
        f"Fixed MUC5B-switch hybrid "
        f"(AUC={auc_hybrid:.4f})"
    ),
)

ax.plot(
    hull_df[
        "fpr"
    ],
    hull_df[
        "tpr"
    ],
    linewidth=3,
    linestyle="--",
    label=(
        f"zedstat convex hull of union "
        f"(AUC={auc_hull:.4f})"
    ),
)

ax.plot(
    [
        0,
        1,
    ],
    [
        0,
        1,
    ],
    linestyle=":",
    linewidth=1,
)

ax.set_xlim(
    0,
    1,
)

ax.set_ylim(
    0,
    1,
)

ax.set_xlabel(
    "False positive rate"
)

ax.set_ylabel(
    "True positive rate"
)

ax.set_title(
    "ZeBRA vs MUC5B-switch hybrid\n"
    "zedstat upper convex hull of concatenated ROC curves"
)

ax.legend(
    loc="lower right"
)

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "ROC_ZEBRA_HYBRID_ZEDSTAT_CONVEX_HULL.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# LOW-FPR ZOOM
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 7)
)

ax.plot(
    roc_zebra[
        "fpr"
    ],
    roc_zebra[
        "tpr"
    ],
    linewidth=2,
    label="ZeBRA",
)

ax.plot(
    roc_hybrid[
        "fpr"
    ],
    roc_hybrid[
        "tpr"
    ],
    linewidth=2,
    label="Fixed MUC5B-switch hybrid",
)

ax.plot(
    hull_df[
        "fpr"
    ],
    hull_df[
        "tpr"
    ],
    linewidth=3,
    linestyle="--",
    label="zedstat union convex hull",
)

ax.set_xlim(
    0,
    0.10,
)

# Let matplotlib choose the useful y-range from the points.
low_mask = (
    hull_df[
        "fpr"
    ]
    <= 0.10
)

if low_mask.any():
    ymax = min(
        1.0,
        float(
            hull_df.loc[
                low_mask,
                "tpr",
            ].max()
        )
        + 0.05,
    )

    ax.set_ylim(
        0,
        ymax,
    )

ax.set_xlabel(
    "False positive rate"
)

ax.set_ylabel(
    "True positive rate"
)

ax.set_title(
    "Low-FPR ROC comparison"
)

ax.legend(
    loc="lower right"
)

fig.tight_layout()

fig.savefig(
    OUTDIR
    / "ROC_ZEBRA_HYBRID_ZEDSTAT_CONVEX_HULL_LOW_FPR.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print("=" * 96)
print("ZEDSTAT ROC AUC RESULTS")
print("=" * 96)

print(
    SUMMARY.to_string(
        index=False
    )
)

print()
print(
    f"ZeBRA AUC                      = "
    f"{auc_zebra:.6f}"
)

print(
    f"Hybrid AUC                     = "
    f"{auc_hybrid:.6f}"
)

print(
    f"Pointwise upper-envelope AUC   = "
    f"{auc_envelope:.6f}"
)

print(
    f"UNION CONVEX-HULL AUC          = "
    f"{auc_hull:.6f}"
)

print()
print(
    "Convex-hull gain over best individual empirical ROC = "
    f"{auc_hull - max(auc_zebra, auc_hybrid):.6f}"
)

print()
print("=" * 96)
print("IMPORTANT INTERPRETATION")
print("=" * 96)

print(
    "The union convex hull answers:"
)

print(
    "  If we are allowed to choose whichever of the two classifiers "
    "gives the superior ROC operating point, what is the upper "
    "achievable ROC envelope?"
)

print()
print(
    "It is NOT the ROC of a single scalar prediction rule unless an "
    "explicit decision policy is constructed that implements the "
    "switching represented by the hull."
)

print()
print(
    "Outputs:"
)

for filename in [
    "OOF_ZEBRA_HYBRID_SCORES.csv",
    "ROC_ZEBRA.csv",
    "ROC_HYBRID.csv",
    "ROC_UNION_RAW_POINTS.csv",
    "ROC_UNION_ZEDSTAT_CONVEX_HULL.csv",
    "ROC_POINTWISE_UPPER_ENVELOPE.csv",
    "ROC_MODEL_WINNER_BY_FPR.csv",
    "ZEDSTAT_ROC_HULL_AUC_SUMMARY.csv",
    "ROC_ZEBRA_HYBRID_ZEDSTAT_CONVEX_HULL.png",
    "ROC_ZEBRA_HYBRID_ZEDSTAT_CONVEX_HULL_LOW_FPR.png",
]:
    print(
        "  ",
        OUTDIR
        / filename,
    )
