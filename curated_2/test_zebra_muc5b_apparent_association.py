#!/usr/bin/env python3
"""Pooled ZeBRA–MUC5B association plots for manuscript framing.

This analysis is intentionally descriptive. It asks whether ZeBRA risk appears
associated with MUC5B rs35705950 T-carrier status in the pooled notebook-32
cohort. It does NOT establish that ZeBRA reconstructs genotype. The companion
`testinverse_muc5b_stratified_by_fild.py` analysis tests whether this pooled
association persists after stratifying on FILD status.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler


BASE_DATA_FILE = "ILD_TOP_DRIVERS_DATA.csv"
TARGET_FILE = "REPHENOTYPES FOR IC.csv"
PRED_FILE = "PREDICTIONS_104W_PRED_WINDOW.parquet"
TARGET_NAME = "FILD or FILA ADJUDICATED"

GG_COLUMN = "rs35705950.1_G_2"
GT_COLUMN = "rs35705950.1_G_1"
TT_COLUMN = "rs35705950.1_G_0"

CV_FOLDS = 5
CV_SEED = 1321
N_SPLINE_KNOTS = 5

OUTDIR = Path("./RESULTS/ZEBRA_MUC5B_APPARENT_ASSOCIATION")
OUTDIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Load the exact notebook-32 cohort sources.
# ---------------------------------------------------------------------

base = pd.read_csv(BASE_DATA_FILE).drop(columns=["target"], errors="ignore")

targets = (
    pd.read_csv(TARGET_FILE)
    .replace({"N": 0, "Y": 1, "n": 0, "y": 1, "na": np.nan})
    .rename(columns={"arb_person_id": "patient_id"})
)

preds = pd.read_parquet(PRED_FILE)

data = (
    base
    .merge(
        targets[["patient_id", TARGET_NAME]].rename(columns={TARGET_NAME: "target"}),
        on="patient_id",
        how="left",
    )
    .merge(
        preds[["patient_id", "predicted_risk"]],
        on="patient_id",
        how="left",
    )
)

data["target"] = (
    pd.to_numeric(data["target"], errors="coerce")
    .fillna(0)
    .astype(int)
)
data = data[data["predicted_risk"].notnull()].copy()


# ---------------------------------------------------------------------
# Decode MUC5B rs35705950. The one-hot encoding used historically is:
#   G_2 = GG, G_1 = GT, G_0 = TT.
# An all-zero row represents an unavailable genotype and is excluded.
# ---------------------------------------------------------------------

required = [GG_COLUMN, GT_COLUMN, TT_COLUMN]
missing = [c for c in required if c not in data.columns]
if missing:
    raise ValueError("Missing required MUC5B columns: " + ", ".join(missing))

for c in required:
    data[c] = pd.to_numeric(data[c], errors="coerce")

data = data[data[required].notnull().all(axis=1)].copy()
onehot_sum = data[required].sum(axis=1)

if (onehot_sum > 1).any():
    raise ValueError("MUC5B one-hot encoding contains multi-hot rows.")

data = data.loc[onehot_sum == 1].copy()
data["MUC5B_T_carrier"] = (
    (data[GT_COLUMN].astype(int) == 1)
    | (data[TT_COLUMN].astype(int) == 1)
).astype(int)

zebra = data["predicted_risk"].astype(float).to_numpy()
y = data["MUC5B_T_carrier"].astype(int).to_numpy()
X = zebra.reshape(-1, 1)


# ---------------------------------------------------------------------
# Compact pooled association statistics. These quantify the *apparent*
# association only; the FILD-stratified companion analysis is required for
# interpretation.
# ---------------------------------------------------------------------

pearson_r, pearson_p = pearsonr(zebra, y)
spearman_r, spearman_p = spearmanr(zebra, y)
raw_auc = roc_auc_score(y, zebra)
direction_independent_auc = max(raw_auc, 1.0 - raw_auc)

cv = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=CV_SEED)

linear_model = make_pipeline(
    StandardScaler(),
    LogisticRegression(C=1e6, solver="lbfgs", max_iter=5000),
)
linear_cv_prob = cross_val_predict(
    linear_model, X, y, cv=cv, method="predict_proba", n_jobs=-1
)[:, 1]
linear_cv_auc = roc_auc_score(y, linear_cv_prob)
linear_cv_logloss = log_loss(y, linear_cv_prob)

spline_model = make_pipeline(
    SplineTransformer(
        n_knots=N_SPLINE_KNOTS,
        degree=3,
        include_bias=False,
    ),
    StandardScaler(),
    LogisticRegression(C=1e6, solver="lbfgs", max_iter=10000),
)
spline_cv_prob = cross_val_predict(
    spline_model, X, y, cv=cv, method="predict_proba", n_jobs=-1
)[:, 1]
spline_cv_auc = roc_auc_score(y, spline_cv_prob)
spline_cv_logloss = log_loss(y, spline_cv_prob)

summary = pd.DataFrame([
    {
        "N_with_called_MUC5B": len(data),
        "MUC5B_T_carrier_N": int(y.sum()),
        "MUC5B_T_carrier_prevalence": float(y.mean()),
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "spearman_rho": float(spearman_r),
        "spearman_p": float(spearman_p),
        "raw_ZeBRA_to_carrier_AUC": float(raw_auc),
        "direction_independent_AUC": float(direction_independent_auc),
        "linear_logistic_CV_AUC": float(linear_cv_auc),
        "linear_logistic_CV_logloss": float(linear_cv_logloss),
        "spline_logistic_CV_AUC": float(spline_cv_auc),
        "spline_logistic_CV_logloss": float(spline_cv_logloss),
        "CV_folds": CV_FOLDS,
        "CV_seed": CV_SEED,
    }
])
summary.to_csv(
    OUTDIR / "ZEBRA_MUC5B_APPARENT_ASSOCIATION_SUMMARY.csv",
    index=False,
)

print("Pooled apparent ZeBRA–MUC5B association")
print(summary.to_string(index=False))
print(
    "\nInterpret only together with testinverse_muc5b_stratified_by_fild.py; "
    "this pooled analysis does not demonstrate genotype reconstruction."
)


# ---------------------------------------------------------------------
# Figure 1: pooled ZeBRA distribution by MUC5B carrier status.
# ---------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot(
    [zebra[y == 0], zebra[y == 1]],
    labels=["GG / no T-risk allele", "GT or TT / T carrier"],
    showmeans=True,
)
ax.set_ylabel("ZeBRA predicted_risk")
ax.set_title("Pooled apparent association: ZeBRA score by MUC5B status")
plt.tight_layout()
plt.savefig(
    OUTDIR / "ZEBRA_BY_MUC5B_T_CARRIER.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close(fig)


# ---------------------------------------------------------------------
# Figure 2: apparent nonlinear pooled relationship. Empirical carrier rates
# are shown across ZeBRA quantile bins, with linear and spline fits. Again,
# the disease-stratified companion analysis is needed to interpret this plot.
# ---------------------------------------------------------------------

linear_model.fit(X, y)
spline_model.fit(X, y)

plot_df = pd.DataFrame({"zebra": zebra, "carrier": y})
plot_df["bin"] = pd.qcut(
    plot_df["zebra"].rank(method="first"),
    q=20,
    labels=False,
    duplicates="drop",
)
empirical = (
    plot_df.groupby("bin")
    .agg(
        mean_zebra=("zebra", "mean"),
        carrier_rate=("carrier", "mean"),
        n=("carrier", "size"),
    )
    .reset_index()
)

grid = np.linspace(float(zebra.min()), float(zebra.max()), 500).reshape(-1, 1)
linear_prob_grid = linear_model.predict_proba(grid)[:, 1]
spline_prob_grid = spline_model.predict_proba(grid)[:, 1]

fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(
    empirical["mean_zebra"],
    empirical["carrier_rate"],
    s=35,
    label="Observed carrier rate (20 quantile bins)",
)
ax.plot(grid.ravel(), linear_prob_grid, linewidth=2, label="Linear logistic")
ax.plot(grid.ravel(), spline_prob_grid, linewidth=2, label="Spline logistic")
ax.set_xlabel("ZeBRA predicted_risk")
ax.set_ylabel("P(MUC5B T-risk carrier)")
ax.set_title("Pooled apparent ZeBRA–MUC5B relationship")
ax.legend()
plt.tight_layout()
plt.savefig(
    OUTDIR / "ZEBRA_MUC5B_NONLINEAR_ASSOCIATION.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close(fig)

print(f"\nSaved pooled apparent-association outputs to {OUTDIR}")
