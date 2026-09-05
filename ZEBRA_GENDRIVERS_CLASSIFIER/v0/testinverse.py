import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import (
    pearsonr,
    spearmanr,
    mannwhitneyu,
    fisher_exact,
)

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_predict,
)
from sklearn.metrics import (
    roc_auc_score,
    confusion_matrix,
)


# ============================================================
# Load data
# ============================================================

BASE_DATA = pd.read_csv(
    "ILD_TOP_DRIVERS_DATA.csv"
).drop(
    columns=["target"],
    errors="ignore",
)

TARGETS = (
    pd.read_csv(
        "REPHENOTYPES FOR IC.csv"
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
            "arb_person_id": "patient_id"
        }
    )
)

PREDS = pd.read_parquet(
    "PREDICTIONS_104W_PRED_WINDOW.parquet"
)

TARGET_NAME = "FILD or FILA ADJUDICATED"


# ============================================================
# EXACT cohort construction used in notebook 32
# ============================================================

DATA = (
    BASE_DATA
    .merge(
        TARGETS[
            ["patient_id", TARGET_NAME]
        ].rename(
            columns={
                TARGET_NAME: "target"
            }
        ),
        on="patient_id",
        how="left",
    )
    .merge(
        PREDS[
            ["patient_id", "predicted_risk"]
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
    DATA["predicted_risk"].notnull()
].copy()


# ============================================================
# Inspect MUC5B rs35705950 genotype columns
# ============================================================

MUC5B_COLUMNS = [
    c for c in DATA.columns
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
        )
    )


# ============================================================
# Expected one-hot encoding
#
# Assuming:
#
#   rs35705950.1_G_2 = GG
#   rs35705950.1_G_1 = GT
#   rs35705950.1_G_0 = TT
#
# Since T is the canonical MUC5B risk allele:
#
#   risk carrier = GT or TT
#                = NOT GG
#
# ============================================================

GG_COLUMN = "rs35705950.1_G_2"
GT_COLUMN = "rs35705950.1_G_1"
TT_COLUMN = "rs35705950.1_G_0"


required = [
    GG_COLUMN,
    GT_COLUMN,
    TT_COLUMN,
]

missing = [
    c for c in required
    if c not in DATA.columns
]

if missing:
    raise ValueError(
        "Missing expected MUC5B columns:\n"
        + "\n".join(missing)
    )


# ============================================================
# Convert genotype-state columns to numeric
# ============================================================

for c in required:
    DATA[c] = pd.to_numeric(
        DATA[c],
        errors="coerce",
    )


# Keep subjects with interpretable MUC5B genotype
DATA = DATA[
    DATA[
        required
    ].notnull().all(
        axis=1
    )
].copy()


# ============================================================
# Sanity check one-hot encoding
# ============================================================

onehot_sum = (
    DATA[
        required
    ].sum(
        axis=1
    )
)

print()
print("One-hot genotype-state sums:")
print(
    onehot_sum.value_counts().sort_index()
)

if not np.all(
    np.isclose(
        onehot_sum,
        1
    )
):
    print()
    print(
        "WARNING: MUC5B genotype columns are not "
        "strictly one-hot for every patient."
    )


# ============================================================
# Define biologically relevant endpoints
# ============================================================

# GG = no T risk allele
muc5b_GG = (
    DATA[
        GG_COLUMN
    ].astype(int)
)

# GT = one T risk allele
muc5b_GT = (
    DATA[
        GT_COLUMN
    ].astype(int)
)

# TT = two T risk alleles
muc5b_TT = (
    DATA[
        TT_COLUMN
    ].astype(int)
)

# PRIMARY endpoint:
# any T risk allele = GT or TT
muc5b_T_carrier = (
    (
        muc5b_GT == 1
    )
    |
    (
        muc5b_TT == 1
    )
).astype(int)


# T dosage:
#
# GG -> 0
# GT -> 1
# TT -> 2
muc5b_T_dosage = (
    muc5b_GT
    + 2 * muc5b_TT
).astype(int)


# ============================================================
# ZeBRA predictor
# ============================================================

zebra = DATA[
    "predicted_risk"
].astype(float).to_numpy()

X = zebra.reshape(
    -1,
    1,
)

y = muc5b_T_carrier.to_numpy(
    dtype=int
)


# ============================================================
# Cohort summary
# ============================================================

print()
print("=" * 70)
print("MUC5B T-RISK CARRIER ANALYSIS")
print("=" * 70)

print(f"N = {len(y):,}")

print(
    f"GG / non-carrier = "
    f"{int((y == 0).sum()):,}"
)

print(
    f"GT or TT / T-carrier = "
    f"{int(y.sum()):,}"
)

print(
    f"T-carrier prevalence = "
    f"{y.mean():.4%}"
)

print()

print("Genotype states:")

print(
    "  GG:",
    int(muc5b_GG.sum()),
)

print(
    "  GT:",
    int(muc5b_GT.sum()),
)

print(
    "  TT:",
    int(muc5b_TT.sum()),
)


# ============================================================
# 1. Correlation:
# ZeBRA vs MUC5B T-carrier
# ============================================================

pearson_r, pearson_p = pearsonr(
    zebra,
    y,
)

spearman_r, spearman_p = spearmanr(
    zebra,
    y,
)

print()
print("Association of ZeBRA with MUC5B T-carrier status")
print("-------------------------------------------------")

print(
    f"Pearson / point-biserial r = "
    f"{pearson_r:.4f}, "
    f"p = {pearson_p:.3e}"
)

print(
    f"Spearman rho               = "
    f"{spearman_r:.4f}, "
    f"p = {spearman_p:.3e}"
)


# ============================================================
# 2. Correlation with T-allele dosage 0/1/2
# ============================================================

dosage = (
    muc5b_T_dosage
    .to_numpy(
        dtype=float
    )
)

pearson_dosage_r, pearson_dosage_p = (
    pearsonr(
        zebra,
        dosage,
    )
)

spearman_dosage_r, spearman_dosage_p = (
    spearmanr(
        zebra,
        dosage,
    )
)

print()
print("Association of ZeBRA with MUC5B T-allele dosage")
print("------------------------------------------------")

print(
    f"Pearson r = "
    f"{pearson_dosage_r:.4f}, "
    f"p = {pearson_dosage_p:.3e}"
)

print(
    f"Spearman rho = "
    f"{spearman_dosage_r:.4f}, "
    f"p = {spearman_dosage_p:.3e}"
)


# ============================================================
# 3. Distribution difference:
# T carriers vs non-carriers
# ============================================================

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

print()
print("ZeBRA distribution by MUC5B T-carrier status")
print("---------------------------------------------")

print(
    f"Mean ZeBRA, T-carrier     = "
    f"{z_carrier.mean():.4f}"
)

print(
    f"Mean ZeBRA, non-carrier   = "
    f"{z_noncarrier.mean():.4f}"
)

print(
    f"Median ZeBRA, T-carrier   = "
    f"{np.median(z_carrier):.4f}"
)

print(
    f"Median ZeBRA, non-carrier = "
    f"{np.median(z_noncarrier):.4f}"
)

print(
    f"Mann-Whitney p            = "
    f"{mw_p:.3e}"
)


# ============================================================
# 4. Raw continuous ZeBRA AUC
#
# Can ZeBRA rank T carriers above non-carriers?
# ============================================================

raw_auc = roc_auc_score(
    y,
    zebra,
)

direction_independent_auc = max(
    raw_auc,
    1 - raw_auc,
)

print()
print("Continuous ZeBRA discrimination")
print("--------------------------------")

print(
    f"Raw ZeBRA -> MUC5B T-carrier AUC = "
    f"{raw_auc:.4f}"
)

print(
    f"Direction-independent AUC         = "
    f"{direction_independent_auc:.4f}"
)


# ============================================================
# 5. 5-fold CV logistic prediction
#
# ONLY predictor = ZeBRA
# ============================================================

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=1321,
)

model = make_pipeline(
    StandardScaler(),
    LogisticRegression(
        max_iter=5000,
    ),
)

p = cross_val_predict(
    model,
    X,
    y,
    cv=cv,
    method="predict_proba",
)[:, 1]

cv_auc = roc_auc_score(
    y,
    p,
)

print()
print(
    f"5-fold CV AUC: "
    f"ZeBRA -> MUC5B T-carrier = "
    f"{cv_auc:.4f}"
)


# ============================================================
# 6. Effect size:
# odds ratio per 1-SD increase in ZeBRA
# ============================================================

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
    penalty=None,
    solver="lbfgs",
    max_iter=5000,
)

effect_model.fit(
    zebra_sd,
    y,
)

beta = (
    effect_model
    .coef_[0, 0]
)

OR = np.exp(
    beta
)

print()
print(
    f"Odds ratio for MUC5B T-carrier "
    f"per 1-SD higher ZeBRA = "
    f"{OR:.3f}"
)


# ============================================================
# 7. High-ZeBRA-risk enrichment analysis
#
# Does restricting to progressively higher ZeBRA risk
# enrich for MUC5B T carriers?
# ============================================================

percentiles = [
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

baseline_prevalence = y.mean()

threshold_rows = []

for pct in percentiles:

    threshold = np.percentile(
        zebra,
        pct,
    )

    high_risk = (
        zebra >= threshold
    ).astype(int)

    tn, fp, fn, tp = confusion_matrix(
        y,
        high_risk,
        labels=[0, 1],
    ).ravel()

    n_high = int(
        high_risk.sum()
    )

    carrier_prevalence_high = (
        y[
            high_risk == 1
        ].mean()
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
        if tp + fn > 0
        else np.nan
    )

    specificity = (
        tn / (tn + fp)
        if tn + fp > 0
        else np.nan
    )

    ppv = (
        tp / (tp + fp)
        if tp + fp > 0
        else np.nan
    )

    binary_auc = roc_auc_score(
        y,
        high_risk,
    )

    odds_ratio, fisher_p = fisher_exact(
        [
            [tp, fp],
            [fn, tn],
        ]
    )

    threshold_rows.append({
        "ZeBRA_percentile_cutoff": pct,
        "ZeBRA_threshold": threshold,
        "N_high_risk": n_high,
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
            odds_ratio,

        "fisher_p":
            fisher_p,
    })


THRESHOLD_RESULTS = pd.DataFrame(
    threshold_rows
)

print()
print("=" * 70)
print("HIGH ZeBRA RISK -> MUC5B T-CARRIER ENRICHMENT")
print("=" * 70)

print(
    THRESHOLD_RESULTS
)

THRESHOLD_RESULTS.to_csv(
    "ZEBRA_HIGH_RISK_MUC5B_T_CARRIER_ENRICHMENT.csv",
    index=False,
)


# ============================================================
# 8. Main summary
# ============================================================

SUMMARY = pd.DataFrame([
    {
        "N": len(y),

        "MUC5B_GG_N":
            int(muc5b_GG.sum()),

        "MUC5B_GT_N":
            int(muc5b_GT.sum()),

        "MUC5B_TT_N":
            int(muc5b_TT.sum()),

        "MUC5B_T_carrier_N":
            int(y.sum()),

        "MUC5B_T_carrier_prevalence":
            y.mean(),

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

        "mean_zebra_T_carrier":
            z_carrier.mean(),

        "mean_zebra_noncarrier":
            z_noncarrier.mean(),

        "median_zebra_T_carrier":
            np.median(z_carrier),

        "median_zebra_noncarrier":
            np.median(z_noncarrier),

        "mannwhitney_p":
            mw_p,

        "raw_auc":
            raw_auc,

        "direction_independent_auc":
            direction_independent_auc,

        "cv_auc":
            cv_auc,

        "OR_per_SD_zebra":
            OR,
    }
])

SUMMARY.to_csv(
    "ZEBRA_MUC5B_T_CARRIER_ASSOCIATION.csv",
    index=False,
)

print(
    SUMMARY
)


# ============================================================
# 9. Plot:
# ZeBRA distribution by carrier status
# ============================================================

fig, ax = plt.subplots(
    figsize=(7, 5)
)

ax.boxplot(
    [
        zebra[y == 0],
        zebra[y == 1],
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
# 10. Plot:
# ZeBRA by genotype dosage
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
# 11. Plot:
# MUC5B T-carrier enrichment at high ZeBRA risk
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
