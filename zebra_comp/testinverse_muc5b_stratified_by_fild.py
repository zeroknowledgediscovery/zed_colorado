#!/usr/bin/env python3
from pathlib import Path

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import (
    pearsonr,
    spearmanr,
    fisher_exact,
)

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_predict,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score,
    confusion_matrix,
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

CV_FOLDS = 5
CV_SEED = 1321

OUTDIR = Path(
    "./RESULTS/MUC5B_STRATIFIED_BY_FILD"
)
OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Helpers
# ============================================================

def safe_name(s):
    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(s),
    ).strip("_")


def odds_ratio_ci(a, b, c, d):
    cells = np.array(
        [a, b, c, d],
        dtype=float,
    )

    if np.any(cells == 0):
        cells += 0.5

    a, b, c, d = cells

    OR = (
        a * d
        / (b * c)
    )

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


def build_threshold_table(
    zebra,
    y,
    percentiles,
):
    baseline = float(
        y.mean()
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

        mask = (
            zebra >= threshold
        )

        n_high = int(
            mask.sum()
        )

        # Extremely small groups are too unstable to interpret.
        if n_high < 10:
            continue

        pred = mask.astype(int)

        tn, fp, fn, tp = confusion_matrix(
            y,
            pred,
            labels=[0, 1],
        ).ravel()

        high_prev = float(
            y[mask].mean()
        )

        enrichment = (
            high_prev / baseline
            if baseline > 0
            else np.nan
        )

        OR, OR_lo, OR_hi = odds_ratio_ci(
            tp,
            fp,
            fn,
            tn,
        )

        fisher_or, fisher_p = fisher_exact(
            [
                [tp, fp],
                [fn, tn],
            ]
        )

        rows.append({
            "ZeBRA_percentile_cutoff":
                pct,

            "ZeBRA_threshold":
                threshold,

            "N_high_risk":
                n_high,

            "actual_high_fraction":
                n_high / len(y),

            "MUC5B_T_carriers_high_risk":
                int(tp),

            "MUC5B_T_carrier_prevalence_high_risk":
                high_prev,

            "MUC5B_T_carrier_prevalence_all":
                baseline,

            "enrichment_fold":
                enrichment,

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

        masks.append({
            "percentile":
                pct,
            "mask":
                mask.copy(),
        })

    return (
        pd.DataFrame(rows),
        masks,
    )


def permutation_max_enrichment(
    y,
    threshold_masks,
    observed_table,
    n_permutations,
    seed,
):
    if len(observed_table) == 0:
        return {
            "observed_max_enrichment":
                np.nan,
            "observed_best_percentile":
                np.nan,
            "permutation_p_max_enrichment":
                np.nan,
            "null_max_enrichment_mean":
                np.nan,
            "null_max_enrichment_q95":
                np.nan,
            "null_max_enrichment_q99":
                np.nan,
            "observed_99_enrichment":
                np.nan,
            "permutation_p_fixed_99":
                np.nan,
        }

    baseline = float(
        y.mean()
    )

    observed_values = (
        observed_table[
            "enrichment_fold"
        ].to_numpy(
            dtype=float
        )
    )

    max_idx = int(
        np.nanargmax(
            observed_values
        )
    )

    observed_max = float(
        observed_values[
            max_idx
        ]
    )

    best_pct = float(
        observed_table.iloc[
            max_idx
        ][
            "ZeBRA_percentile_cutoff"
        ]
    )

    # Locate a 99th-percentile row if it exists.
    p99_rows = observed_table[
        observed_table[
            "ZeBRA_percentile_cutoff"
        ] == 99
    ]

    if len(p99_rows):
        observed_99 = float(
            p99_rows.iloc[0][
                "enrichment_fold"
            ]
        )
    else:
        observed_99 = np.nan

    rng = np.random.default_rng(
        seed
    )

    perm_max = np.empty(
        n_permutations,
        dtype=float,
    )

    perm_99 = np.full(
        n_permutations,
        np.nan,
        dtype=float,
    )

    mask_by_pct = {
        float(x["percentile"]):
            x["mask"]
        for x in threshold_masks
    }

    print(
        f"  running {n_permutations:,} "
        "label permutations..."
    )

    for b in tqdm(
        range(
            n_permutations
        ),
        desc="  permutations",
        leave=False,
    ):
        yp = rng.permutation(
            y
        )

        vals = []

        for item in threshold_masks:
            mask = item[
                "mask"
            ]

            vals.append(
                float(
                    yp[
                        mask
                    ].mean()
                ) / baseline
            )

        vals = np.asarray(
            vals,
            dtype=float,
        )

        perm_max[b] = (
            np.nanmax(
                vals
            )
        )

        if 99.0 in mask_by_pct:
            perm_99[b] = (
                float(
                    yp[
                        mask_by_pct[
                            99.0
                        ]
                    ].mean()
                )
                / baseline
            )

    p_max = (
        1
        + np.sum(
            perm_max
            >= observed_max
        )
    ) / (
        n_permutations
        + 1
    )

    if np.isfinite(
        observed_99
    ):
        p_99 = (
            1
            + np.sum(
                perm_99
                >= observed_99
            )
        ) / (
            np.sum(
                np.isfinite(
                    perm_99
                )
            )
            + 1
        )
    else:
        p_99 = np.nan

    return {
        "observed_max_enrichment":
            observed_max,

        "observed_best_percentile":
            best_pct,

        "permutation_p_max_enrichment":
            p_max,

        "null_max_enrichment_mean":
            float(
                perm_max.mean()
            ),

        "null_max_enrichment_q95":
            float(
                np.quantile(
                    perm_max,
                    0.95,
                )
            ),

        "null_max_enrichment_q99":
            float(
                np.quantile(
                    perm_max,
                    0.99,
                )
            ),

        "observed_99_enrichment":
            observed_99,

        "permutation_p_fixed_99":
            p_99,
    }


def analyze_group(
    group_name,
    d,
    permutation_seed,
):
    print()
    print("=" * 88)
    print(group_name)
    print("=" * 88)

    if len(d) < 50:
        print(
            "Too few patients; skipping."
        )
        return None

    zebra = (
        d[
            "predicted_risk"
        ]
        .astype(float)
        .to_numpy()
    )

    y = (
        d[
            "MUC5B_T_carrier"
        ]
        .astype(int)
        .to_numpy()
    )

    if np.unique(
        y
    ).size < 2:
        print(
            "Only one MUC5B carrier class; skipping."
        )
        return None

    n_carrier = int(
        y.sum()
    )

    n_noncarrier = int(
        (y == 0).sum()
    )

    prevalence = float(
        y.mean()
    )

    print(
        f"N={len(y):,}; "
        f"carriers={n_carrier:,}; "
        f"non-carriers={n_noncarrier:,}; "
        f"carrier prevalence={prevalence:.3%}"
    )

    # --------------------------------------------------------
    # Global continuous association
    # --------------------------------------------------------

    point_r, point_p = pearsonr(
        zebra,
        y,
    )

    spear_r, spear_p = spearmanr(
        zebra,
        y,
    )

    raw_auc = roc_auc_score(
        y,
        zebra,
    )

    # --------------------------------------------------------
    # CV logistic ZeBRA -> MUC5B carrier
    # --------------------------------------------------------

    minority = min(
        n_carrier,
        n_noncarrier,
    )

    n_splits = min(
        CV_FOLDS,
        minority,
    )

    if n_splits >= 2:
        cv = StratifiedKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=CV_SEED,
        )

        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                C=1e6,
                solver="lbfgs",
                max_iter=5000,
            ),
        )

        cv_prob = cross_val_predict(
            model,
            zebra.reshape(-1, 1),
            y,
            cv=cv,
            method="predict_proba",
            n_jobs=-1,
        )[:, 1]

        cv_auc = roc_auc_score(
            y,
            cv_prob,
        )
    else:
        cv_auc = np.nan

    # --------------------------------------------------------
    # High-risk enrichment
    # --------------------------------------------------------

    threshold_table, masks = build_threshold_table(
        zebra,
        y,
        PERCENTILES,
    )

    if len(
        threshold_table
    ):
        threshold_table.insert(
            0,
            "group",
            group_name,
        )

    perm = permutation_max_enrichment(
        y,
        masks,
        threshold_table,
        N_PERMUTATIONS,
        permutation_seed,
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary = {
        "group":
            group_name,

        "N":
            len(y),

        "MUC5B_T_carrier_N":
            n_carrier,

        "MUC5B_T_carrier_prevalence":
            prevalence,

        "point_biserial_r":
            point_r,

        "point_biserial_p":
            point_p,

        "spearman_rho":
            spear_r,

        "spearman_p":
            spear_p,

        "continuous_auc":
            raw_auc,

        "cv_logistic_auc":
            cv_auc,

        **perm,
    }

    print()
    print(
        f"Continuous AUC: {raw_auc:.4f}"
    )
    print(
        f"CV logistic AUC: {cv_auc:.4f}"
    )
    print(
        f"Point-biserial r={point_r:.4f}, "
        f"p={point_p:.3e}"
    )

    if np.isfinite(
        perm[
            "observed_max_enrichment"
        ]
    ):
        print(
            "Max tail enrichment: "
            f"{perm['observed_max_enrichment']:.3f}x "
            f"at percentile "
            f"{perm['observed_best_percentile']:g}; "
            f"permutation p="
            f"{perm['permutation_p_max_enrichment']:.4g}"
        )

    # --------------------------------------------------------
    # Plot enrichment curve for this group
    # --------------------------------------------------------

    if len(
        threshold_table
    ):
        fig, ax = plt.subplots(
            figsize=(7, 5)
        )

        ax.plot(
            threshold_table[
                "ZeBRA_percentile_cutoff"
            ],
            threshold_table[
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
            group_name
            + "\nMUC5B carrier enrichment at high ZeBRA risk"
        )

        fig.tight_layout()

        fig.savefig(
            OUTDIR
            / (
                safe_name(
                    group_name
                )
                + "_MUC5B_enrichment.png"
            ),
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(
            fig
        )

    return (
        summary,
        threshold_table,
    )


# ============================================================
# Load data
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
# Construct notebook-32 cohort while preserving whether the
# adjudication value was actually observed before fillna(0)
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

# EXACT notebook-32 ZeBRA availability filter
DATA = DATA[
    DATA[
        "predicted_risk"
    ].notnull()
].copy()


# ============================================================
# Decode MUC5B one-hot genotype and drop missing genotype rows
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

n_zero = int(
    (
        dummy_sum == 0
    ).sum()
)

n_one = int(
    (
        dummy_sum == 1
    ).sum()
)

n_gt1 = int(
    (
        dummy_sum > 1
    ).sum()
)

print()
print(
    "MUC5B genotype dummy-state sums:"
)
print(
    dummy_sum.value_counts()
    .sort_index()
    .to_string()
)
print()
print(
    f"valid called genotype: {n_one:,}"
)
print(
    f"missing genotype:      {n_zero:,}"
)
print(
    f"invalid multi-hot:     {n_gt1:,}"
)

if n_gt1 > 0:
    raise ValueError(
        "Some MUC5B rows are multi-hot."
    )

# Only called genotype
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


# ============================================================
# Define the four scientifically useful strata
# ============================================================

GROUPS = {
    # Reference result
    "Full notebook-32 cohort":
        DATA,

    # This is the direct current notebook-32 negative class.
    # IMPORTANT: it includes patients whose adjudication was missing
    # and was converted to 0 by notebook 32.
    "Notebook-32 target=0 (includes unadjudicated)":
        DATA[
            DATA[
                "target"
            ] == 0
        ].copy(),

    # Stronger test: only patients with an actually observed
    # adjudication explicitly coded negative.
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

    # Positive comparison group.
    "FILD/FILA positive":
        DATA[
            DATA[
                "target"
            ] == 1
        ].copy(),
}


# ============================================================
# Show strata sizes before running
# ============================================================

print()
print("=" * 88)
print("STRATA")
print("=" * 88)

for name, d in GROUPS.items():
    print(
        f"{name}: "
        f"N={len(d):,}, "
        f"MUC5B carriers="
        f"{int(d['MUC5B_T_carrier'].sum()):,}, "
        f"FILD positives="
        f"{int(d['target'].sum()):,}"
    )


# ============================================================
# Analyze each stratum
# ============================================================

ALL_SUMMARY = []
ALL_THRESHOLDS = []

for idx, (
    name,
    d,
) in enumerate(
    GROUPS.items(),
    start=1,
):

    result = analyze_group(
        name,
        d,
        permutation_seed=(
            PERMUTATION_SEED
            + idx * 1000
        ),
    )

    if result is None:
        continue

    summary, threshold_table = result

    ALL_SUMMARY.append(
        summary
    )

    if len(
        threshold_table
    ):
        ALL_THRESHOLDS.append(
            threshold_table
        )


# ============================================================
# Save results
# ============================================================

SUMMARY_DF = pd.DataFrame(
    ALL_SUMMARY
)

SUMMARY_DF.to_csv(
    OUTDIR
    / "MUC5B_ZEBRA_STRATIFIED_SUMMARY.csv",
    index=False,
)

if ALL_THRESHOLDS:
    THRESHOLD_DF = pd.concat(
        ALL_THRESHOLDS,
        ignore_index=True,
    )

    THRESHOLD_DF.to_csv(
        OUTDIR
        / "MUC5B_ZEBRA_STRATIFIED_THRESHOLDS.csv",
        index=False,
    )
else:
    THRESHOLD_DF = pd.DataFrame()


# ============================================================
# Compact comparison plot:
# maximum tail enrichment across groups
# ============================================================

plot_summary = SUMMARY_DF[
    SUMMARY_DF[
        "observed_max_enrichment"
    ].notnull()
].copy()

if len(
    plot_summary
):
    fig, ax = plt.subplots(
        figsize=(10, 5)
    )

    x = np.arange(
        len(
            plot_summary
        )
    )

    ax.bar(
        x,
        plot_summary[
            "observed_max_enrichment"
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
        plot_summary[
            "group"
        ],
        rotation=25,
        ha="right",
    )

    ax.set_ylabel(
        "Maximum MUC5B T-carrier enrichment"
    )

    ax.set_title(
        "Extreme ZeBRA-tail MUC5B enrichment by FILD stratum"
    )

    fig.tight_layout()

    fig.savefig(
        OUTDIR
        / "MUC5B_enrichment_by_FILD_stratum.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# ============================================================
# Final console report
# ============================================================

print()
print("=" * 88)
print("FINAL STRATIFIED SUMMARY")
print("=" * 88)

print(
    SUMMARY_DF.to_string(
        index=False
    )
)

print()
print(
    "Primary question:"
)
print(
    "Does the MUC5B tail enrichment persist when FILD/FILA-positive "
    "patients are removed?"
)
print()
print(
    "Focus especially on these rows:"
)
print(
    "  1. Notebook-32 target=0 (includes unadjudicated)"
)
print(
    "  2. Explicitly adjudicated FILD/FILA negative"
)
print()
print(
    "and on:"
)
print(
    "  observed_max_enrichment"
)
print(
    "  observed_best_percentile"
)
print(
    "  permutation_p_max_enrichment"
)
print(
    "  observed_99_enrichment"
)
print(
    "  permutation_p_fixed_99"
)

print()
print(
    "Outputs saved under:"
)
print(
    OUTDIR.resolve()
)
