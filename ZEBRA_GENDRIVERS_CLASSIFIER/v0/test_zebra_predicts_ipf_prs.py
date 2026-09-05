#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr

from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_predict
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.feature_selection import mutual_info_regression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


# ============================================================
# Configuration
# ============================================================

BASE_DATA_FILE = "ILD_TOP_DRIVERS_DATA.csv"
TARGET_FILE = "REPHENOTYPES FOR IC.csv"
PRED_FILE = "PREDICTIONS_104W_PRED_WINDOW.parquet"

TARGET_NAME = "FILD or FILA ADJUDICATED"

OUTDIR = Path("./RESULTS/ZEBRA_PREDICTS_IPF_PRS")
OUTDIR.mkdir(parents=True, exist_ok=True)

CV_FOLDS = 5
CV_SEED = 1321

# ------------------------------------------------------------
# Literature-weighted IPF genetic risk score
#
# ORs from published IPF meta-analysis:
#
# rs35705950 / MUC5B:
#   T vs reference, OR = 3.85
#
# rs2736100 / TERT:
#   C vs A, OR = 0.70
#   Therefore A vs C risk OR = 1 / 0.70
#
# rs2609255 / FAM13A:
#   G vs T, OR = 1.37
#
# rs2076295 / DSP:
#   G vs T, OR = 1.31
#
# rs12610495 / DPP9:
#   G vs A, OR = 1.29
#
# The PRS weight is log(OR) per risk allele.
#
# "base" is the allele whose dosage is represented by the
# _0/_1/_2 dummy columns in ILD_TOP_DRIVERS_DATA.csv.
# ------------------------------------------------------------

LOCI = [
    {
        "rsid": "rs35705950",
        "gene": "MUC5B",
        "prefix": "rs35705950.1_G",
        "base": "G",
        "risk_allele": "T",
        "risk_or": 3.85,
    },
    {
        "rsid": "rs2736100",
        "gene": "TERT",
        "prefix": "rs2736100_A",
        "base": "A",
        "risk_allele": "A",
        "risk_or": 1.0 / 0.70,
    },
    {
        "rsid": "rs2609255",
        "gene": "FAM13A",
        "prefix": "rs2609255_T",
        "base": "T",
        "risk_allele": "G",
        "risk_or": 1.37,
    },
    {
        "rsid": "rs2076295",
        "gene": "DSP",
        "prefix": "rs2076295_T",
        "base": "T",
        "risk_allele": "G",
        "risk_or": 1.31,
    },
    {
        "rsid": "rs12610495",
        "gene": "DPP9",
        "prefix": "rs12610495_A",
        "base": "A",
        "risk_allele": "G",
        "risk_or": 1.29,
    },
]


# ============================================================
# Helpers
# ============================================================

def decode_base_dosage(df, prefix):
    """
    Decode the original 0/1/2 allele dosage from three dummy
    columns prefix_0, prefix_1, prefix_2.

    pd.get_dummies() encodes a missing original genotype as
    0,0,0. Those rows are returned as NaN.
    """
    cols = [
        f"{prefix}_0",
        f"{prefix}_1",
        f"{prefix}_2",
    ]

    missing_cols = [
        c for c in cols
        if c not in df.columns
    ]

    if missing_cols:
        raise ValueError(
            f"Missing genotype columns for {prefix}: "
            + ", ".join(missing_cols)
        )

    x = (
        df[cols]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
    )

    row_sum = x.sum(
        axis=1,
        min_count=1,
    )

    invalid = (
        row_sum > 1
    )

    if invalid.any():
        raise ValueError(
            f"{prefix}: {int(invalid.sum())} rows "
            "have >1 active dummy state."
        )

    dosage = (
        0 * x[cols[0]]
        + 1 * x[cols[1]]
        + 2 * x[cols[2]]
    ).astype(float)

    # all-zero dummy state = missing original genotype
    dosage[
        row_sum == 0
    ] = np.nan

    return dosage


def construct_prs(df):
    out = df.copy()

    weight_rows = []

    for locus in LOCI:
        base_dosage = decode_base_dosage(
            out,
            locus["prefix"],
        )

        if (
            locus["risk_allele"]
            == locus["base"]
        ):
            risk_dosage = (
                base_dosage
            )
            conversion = (
                f"{locus['base']} dosage"
            )
        else:
            # biallelic locus:
            # dosage(other allele) = 2 - dosage(base allele)
            risk_dosage = (
                2.0
                - base_dosage
            )
            conversion = (
                f"2 - {locus['base']} dosage"
            )

        beta = float(
            np.log(
                locus["risk_or"]
            )
        )

        dosage_col = (
            f"{locus['rsid']}_"
            f"{locus['gene']}_risk_dosage"
        )

        weighted_col = (
            f"{locus['rsid']}_"
            f"{locus['gene']}_weighted"
        )

        out[
            dosage_col
        ] = risk_dosage

        out[
            weighted_col
        ] = (
            risk_dosage
            * beta
        )

        weight_rows.append({
            **locus,
            "beta_log_OR":
                beta,
            "dosage_conversion":
                conversion,
            "N_called":
                int(
                    risk_dosage
                    .notnull()
                    .sum()
                ),
            "N_missing":
                int(
                    risk_dosage
                    .isnull()
                    .sum()
                ),
        })

    weights = pd.DataFrame(
        weight_rows
    )

    weighted_cols = [
        f"{x['rsid']}_{x['gene']}_weighted"
        for x in LOCI
    ]

    no_muc5b_cols = [
        c
        for c in weighted_cols
        if "MUC5B" not in c
    ]

    dosage_cols = [
        f"{x['rsid']}_{x['gene']}_risk_dosage"
        for x in LOCI
    ]

    no_muc5b_dosage_cols = [
        c
        for c in dosage_cols
        if "MUC5B" not in c
    ]

    # --------------------------------------------------------
    # Complete-case scores.
    #
    # Do not silently impute missing genotypes for this first
    # analysis. Missingness is retained and complete cases are
    # identified explicitly.
    # --------------------------------------------------------

    out[
        "PRS5_complete"
    ] = out[
        weighted_cols
    ].notnull().all(
        axis=1
    )

    out[
        "PRS4_noMUC5B_complete"
    ] = out[
        no_muc5b_cols
    ].notnull().all(
        axis=1
    )

    out[
        "PRS5"
    ] = out[
        weighted_cols
    ].sum(
        axis=1,
        min_count=len(
            weighted_cols
        ),
    )

    out[
        "PRS4_noMUC5B"
    ] = out[
        no_muc5b_cols
    ].sum(
        axis=1,
        min_count=len(
            no_muc5b_cols
        ),
    )

    # Unweighted risk-allele counts are retained as a sensitivity
    # analysis.
    out[
        "risk_allele_count_5"
    ] = out[
        dosage_cols
    ].sum(
        axis=1,
        min_count=len(
            dosage_cols
        ),
    )

    out[
        "risk_allele_count_4_noMUC5B"
    ] = out[
        no_muc5b_dosage_cols
    ].sum(
        axis=1,
        min_count=len(
            no_muc5b_dosage_cols
        ),
    )

    # Z-score each PRS using only complete observations.
    for score in [
        "PRS5",
        "PRS4_noMUC5B",
    ]:
        valid = out[
            score
        ].notnull()

        mu = out.loc[
            valid,
            score,
        ].mean()

        sd = out.loc[
            valid,
            score,
        ].std()

        out[
            score + "_z"
        ] = (
            out[
                score
            ]
            - mu
        ) / sd

    return (
        out,
        weights,
    )


def association_summary(
    df,
    score_col,
    group_name,
):
    work = df[
        [
            "predicted_risk",
            score_col,
        ]
    ].dropna().copy()

    zebra = (
        work[
            "predicted_risk"
        ]
        .astype(float)
        .to_numpy()
    )

    prs = (
        work[
            score_col
        ]
        .astype(float)
        .to_numpy()
    )

    pearson_r, pearson_p = pearsonr(
        zebra,
        prs,
    )

    spearman_r, spearman_p = spearmanr(
        zebra,
        prs,
    )

    # --------------------------------------------------------
    # Cross-validated continuous prediction:
    # PRS ~ ZeBRA
    # --------------------------------------------------------

    X = zebra.reshape(
        -1,
        1,
    )

    cv = KFold(
        n_splits=CV_FOLDS,
        shuffle=True,
        random_state=CV_SEED,
    )

    p_prs = cross_val_predict(
        LinearRegression(),
        X,
        prs,
        cv=cv,
        n_jobs=-1,
    )

    cv_r2 = r2_score(
        prs,
        p_prs,
    )

    # standardized slope
    z_zebra = (
        (
            zebra
            - zebra.mean()
        )
        / zebra.std()
    ).reshape(
        -1,
        1,
    )

    z_prs = (
        prs
        - prs.mean()
    ) / prs.std()

    lm = LinearRegression().fit(
        z_zebra,
        z_prs,
    )

    standardized_beta = float(
        lm.coef_[0]
    )

    # Non-parametric dependence estimate
    mi = float(
        mutual_info_regression(
            X,
            prs,
            random_state=CV_SEED,
        )[0]
    )

    # --------------------------------------------------------
    # Can ZeBRA identify individuals in the top 10% of PRS?
    # --------------------------------------------------------

    high_prs_threshold = float(
        np.quantile(
            prs,
            0.90,
        )
    )

    high_prs = (
        prs >= high_prs_threshold
    ).astype(int)

    auc_high_prs_raw = roc_auc_score(
        high_prs,
        zebra,
    )

    cv_strat = StratifiedKFold(
        n_splits=CV_FOLDS,
        shuffle=True,
        random_state=CV_SEED,
    )

    classifier = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=5000,
        ),
    )

    p_high = cross_val_predict(
        classifier,
        X,
        high_prs,
        cv=cv_strat,
        method="predict_proba",
        n_jobs=-1,
    )[:, 1]

    auc_high_prs_cv = roc_auc_score(
        high_prs,
        p_high,
    )

    return {
        "group":
            group_name,

        "PRS":
            score_col,

        "N":
            len(work),

        "pearson_r":
            pearson_r,

        "pearson_p":
            pearson_p,

        "spearman_rho":
            spearman_r,

        "spearman_p":
            spearman_p,

        "standardized_linear_beta":
            standardized_beta,

        "cv_linear_R2":
            cv_r2,

        "mutual_information":
            mi,

        "top10_PRS_threshold":
            high_prs_threshold,

        "top10_PRS_N":
            int(
                high_prs.sum()
            ),

        "raw_ZeBRA_AUC_for_top10_PRS":
            auc_high_prs_raw,

        "cv_ZeBRA_AUC_for_top10_PRS":
            auc_high_prs_cv,
    }


def zebra_quantile_table(
    df,
    score_col,
    group_name,
):
    work = df[
        [
            "predicted_risk",
            score_col,
        ]
    ].dropna().copy()

    work[
        "ZeBRA_decile"
    ] = pd.qcut(
        work[
            "predicted_risk"
        ].rank(
            method="first"
        ),
        q=10,
        labels=False,
    ) + 1

    table = (
        work
        .groupby(
            "ZeBRA_decile"
        )
        .agg(
            N=(
                score_col,
                "size",
            ),
            mean_ZeBRA=(
                "predicted_risk",
                "mean",
            ),
            mean_PRS=(
                score_col,
                "mean",
            ),
            median_PRS=(
                score_col,
                "median",
            ),
            sd_PRS=(
                score_col,
                "std",
            ),
        )
        .reset_index()
    )

    table.insert(
        0,
        "group",
        group_name,
    )

    table.insert(
        1,
        "PRS",
        score_col,
    )

    return table


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

DATA = DATA[
    DATA[
        "predicted_risk"
    ].notnull()
].copy()

print(
    f"Notebook-32 cohort N = {len(DATA):,}"
)
print(
    f"FILD/FILA positives = "
    f"{int(DATA['target'].sum()):,}"
)


# ============================================================
# Construct externally weighted PRSs
# ============================================================

DATA, WEIGHTS = construct_prs(
    DATA
)

WEIGHTS.to_csv(
    OUTDIR
    / "IPF_PRS_weights_used.csv",
    index=False,
)

print()
print(
    "PRS weights used:"
)
print(
    WEIGHTS[
        [
            "rsid",
            "gene",
            "prefix",
            "base",
            "risk_allele",
            "risk_or",
            "beta_log_OR",
            "dosage_conversion",
            "N_called",
            "N_missing",
        ]
    ].to_string(
        index=False
    )
)

print()
print(
    "Complete PRS counts:"
)
print(
    "  PRS5:",
    int(
        DATA[
            "PRS5"
        ].notnull().sum()
    ),
)
print(
    "  PRS4_noMUC5B:",
    int(
        DATA[
            "PRS4_noMUC5B"
        ].notnull().sum()
    ),
)


# ============================================================
# Save patient-level scores
# ============================================================

patient_cols = [
    "patient_id",
    "target",
    "target_was_observed",
    "predicted_risk",
    "PRS5",
    "PRS5_z",
    "PRS4_noMUC5B",
    "PRS4_noMUC5B_z",
    "risk_allele_count_5",
    "risk_allele_count_4_noMUC5B",
]

for locus in LOCI:
    patient_cols.append(
        f"{locus['rsid']}_"
        f"{locus['gene']}_risk_dosage"
    )

DATA[
    patient_cols
].to_csv(
    OUTDIR
    / "patient_level_ZeBRA_and_IPF_PRS.csv",
    index=False,
)


# ============================================================
# Define analysis groups
# ============================================================

GROUPS = {
    "Full notebook-32 cohort":
        DATA,

    "Notebook-32 target=0":
        DATA[
            DATA[
                "target"
            ] == 0
        ].copy(),

    "Explicitly adjudicated FILD/FILA negative":
        DATA[
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
        ].copy(),

    "FILD/FILA positive":
        DATA[
            DATA[
                "target"
            ] == 1
        ].copy(),
}


# ============================================================
# Does ZeBRA predict the PRS?
# ============================================================

SUMMARY_ROWS = []
DECILE_TABLES = []

for group_name, group_df in GROUPS.items():

    print()
    print(
        "=" * 88
    )
    print(
        group_name
    )
    print(
        "=" * 88
    )

    for score_col in [
        "PRS5_z",
        "PRS4_noMUC5B_z",
    ]:
        result = association_summary(
            group_df,
            score_col,
            group_name,
        )

        SUMMARY_ROWS.append(
            result
        )

        deciles = zebra_quantile_table(
            group_df,
            score_col,
            group_name,
        )

        DECILE_TABLES.append(
            deciles
        )

        print()
        print(
            score_col
        )
        print(
            f"  N = {result['N']:,}"
        )
        print(
            f"  Pearson r = "
            f"{result['pearson_r']:.4f}, "
            f"p={result['pearson_p']:.3e}"
        )
        print(
            f"  Spearman rho = "
            f"{result['spearman_rho']:.4f}, "
            f"p={result['spearman_p']:.3e}"
        )
        print(
            f"  CV linear R^2 = "
            f"{result['cv_linear_R2']:.5f}"
        )
        print(
            f"  Mutual information = "
            f"{result['mutual_information']:.5f}"
        )
        print(
            f"  ZeBRA CV AUC for top-10% PRS = "
            f"{result['cv_ZeBRA_AUC_for_top10_PRS']:.4f}"
        )


SUMMARY = pd.DataFrame(
    SUMMARY_ROWS
)

SUMMARY.to_csv(
    OUTDIR
    / "ZEBRA_PREDICTS_IPF_PRS_SUMMARY.csv",
    index=False,
)

DECILES = pd.concat(
    DECILE_TABLES,
    ignore_index=True,
)

DECILES.to_csv(
    OUTDIR
    / "PRS_BY_ZEBRA_DECILE.csv",
    index=False,
)


# ============================================================
# Primary plots: full cohort
# ============================================================

for score_col, title in [
    (
        "PRS5_z",
        "5-locus IPF PRS",
    ),
    (
        "PRS4_noMUC5B_z",
        "4-locus IPF PRS excluding MUC5B",
    ),
]:

    work = DATA[
        [
            "predicted_risk",
            score_col,
        ]
    ].dropna()

    # Scatter + linear relationship
    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    ax.scatter(
        work[
            "predicted_risk"
        ],
        work[
            score_col
        ],
        s=8,
        alpha=0.15,
    )

    x_grid = np.linspace(
        work[
            "predicted_risk"
        ].min(),
        work[
            "predicted_risk"
        ].max(),
        200,
    )

    lm = LinearRegression().fit(
        work[
            [
                "predicted_risk"
            ]
        ],
        work[
            score_col
        ],
    )

    ax.plot(
        x_grid,
        lm.predict(
            x_grid.reshape(
                -1,
                1,
            )
        ),
        linewidth=2,
    )

    ax.set_xlabel(
        "ZeBRA predicted_risk"
    )
    ax.set_ylabel(
        "Standardized IPF genetic risk score"
    )
    ax.set_title(
        f"ZeBRA vs {title}"
    )

    fig.tight_layout()

    fig.savefig(
        OUTDIR
        / (
            score_col
            + "_vs_ZeBRA_scatter.png"
        ),
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )

    # Mean PRS by ZeBRA decile
    dec = DECILES[
        (
            DECILES[
                "group"
            ]
            == "Full notebook-32 cohort"
        )
        &
        (
            DECILES[
                "PRS"
            ]
            == score_col
        )
    ]

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    ax.plot(
        dec[
            "ZeBRA_decile"
        ],
        dec[
            "mean_PRS"
        ],
        marker="o",
    )

    ax.axhline(
        0.0,
        linestyle="--",
    )

    ax.set_xlabel(
        "ZeBRA risk decile"
    )
    ax.set_ylabel(
        "Mean standardized IPF genetic risk score"
    )
    ax.set_title(
        f"{title} across ZeBRA risk deciles"
    )

    fig.tight_layout()

    fig.savefig(
        OUTDIR
        / (
            score_col
            + "_by_ZeBRA_decile.png"
        ),
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# ============================================================
# Final report
# ============================================================

print()
print(
    "=" * 88
)
print(
    "FINAL SUMMARY"
)
print(
    "=" * 88
)

print(
    SUMMARY.to_string(
        index=False
    )
)

print()
print(
    "Primary interpretation:"
)
print(
    "  PRS5_z tests whether ZeBRA tracks a compact IPF genetic-risk score "
    "including MUC5B."
)
print(
    "  PRS4_noMUC5B_z tests whether any relationship persists after removing "
    "the dominant MUC5B locus."
)
print(
    "  The strongest evidence that ZeBRA carries distributed genomic information "
    "would be reproducible Pearson/Spearman association, positive out-of-sample "
    "R^2, and AUC >0.5 for identifying the top PRS decile."
)

print()
print(
    "Outputs:"
)
for filename in [
    "IPF_PRS_weights_used.csv",
    "patient_level_ZeBRA_and_IPF_PRS.csv",
    "ZEBRA_PREDICTS_IPF_PRS_SUMMARY.csv",
    "PRS_BY_ZEBRA_DECILE.csv",
    "PRS5_z_vs_ZeBRA_scatter.png",
    "PRS4_noMUC5B_z_vs_ZeBRA_scatter.png",
    "PRS5_z_by_ZeBRA_decile.png",
    "PRS4_noMUC5B_z_by_ZeBRA_decile.png",
]:
    print(
        " ",
        OUTDIR
        / filename
    )
