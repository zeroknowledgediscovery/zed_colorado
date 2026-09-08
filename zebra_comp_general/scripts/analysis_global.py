#!/usr/bin/env python3
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler

from common import metric_summary


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def safe_auc(y: np.ndarray, score: np.ndarray) -> float:
    return float(roc_auc_score(y, score)) if np.unique(y).size == 2 else np.nan


def safe_ap(y: np.ndarray, score: np.ndarray) -> float:
    return float(average_precision_score(y, score)) if np.unique(y).size == 2 else np.nan


def probs_metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    eps = 1e-12
    p = np.clip(np.asarray(p, float), eps, 1 - eps)
    return {
        "auc": safe_auc(y, p),
        "average_precision": safe_ap(y, p),
        "brier": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
    }


def _lightgbm_or_logistic(seed: int):
    try:
        from lightgbm import LGBMClassifier
        return LGBMClassifier(
            n_estimators=250,
            learning_rate=0.04,
            num_leaves=15,
            subsample=0.9,
            colsample_bytree=0.8,
            reg_lambda=1.0,
            random_state=seed,
            verbosity=-1,
        )
    except Exception:
        return LogisticRegression(C=1.0, solver="lbfgs", max_iter=5000, random_state=seed)


def _fit_predict_matrix(Xtr: pd.DataFrame, ytr: np.ndarray, Xte: pd.DataFrame, seed: int, nonlinear: bool = True) -> np.ndarray:
    imputer = SimpleImputer(strategy="median")
    xtr = imputer.fit_transform(Xtr)
    xte = imputer.transform(Xte)
    if nonlinear:
        model = _lightgbm_or_logistic(seed)
        if isinstance(model, LogisticRegression):
            scaler = StandardScaler()
            xtr = scaler.fit_transform(xtr)
            xte = scaler.transform(xte)
    else:
        scaler = StandardScaler()
        xtr = scaler.fit_transform(xtr)
        xte = scaler.transform(xte)
        model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=5000, random_state=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.fit(xtr, ytr)
    return model.predict_proba(xte)[:, 1]


def _operating_point(y_train: np.ndarray, s_train: np.ndarray, y_test: np.ndarray, s_test: np.ndarray, fpr: float) -> dict[str, float]:
    neg = np.asarray(s_train)[np.asarray(y_train) == 0]
    if len(neg) == 0:
        return {"threshold": np.nan, "test_fpr": np.nan, "test_sensitivity": np.nan, "lr_plus": np.nan}
    thr = float(np.quantile(neg, 1 - fpr, method="higher"))
    pred = np.asarray(s_test) >= thr
    y = np.asarray(y_test)
    controls = y == 0
    cases = y == 1
    tfpr = float(pred[controls].mean()) if controls.any() else np.nan
    sens = float(pred[cases].mean()) if cases.any() else np.nan
    lr = sens / tfpr if np.isfinite(tfpr) and tfpr > 0 else np.inf if sens > 0 else np.nan
    return {"threshold": thr, "test_fpr": tfpr, "test_sensitivity": sens, "lr_plus": lr}


def run_global_comparison(d: pd.DataFrame, feature_cols: list[str], outroot: Path, params: dict[str, Any], cfg: dict[str, Any]) -> None:
    out = ensure_dir(outroot / "01_GLOBAL_COMPARISON")
    nrep = int(params["global_repeats"])
    train_size = float(params["train_size"])
    seed = int(params["random_state"])
    fprs = cfg.get("analysis", {}).get("operating_fprs", [0.005, 0.01, 0.02, 0.05, 0.10])
    y = d.target.to_numpy(int)
    splitter = StratifiedShuffleSplit(n_splits=nrep, train_size=train_size, random_state=seed)
    rows = []
    ops = []
    Xg = d[feature_cols]
    Xc = pd.concat([d[["zebra_score"]], Xg], axis=1)
    for r, (tr, te) in enumerate(splitter.split(np.zeros(len(y)), y)):
        yt, ye = y[tr], y[te]
        ztr, zte = d.zebra_score.to_numpy(float)[tr], d.zebra_score.to_numpy(float)[te]
        pg = _fit_predict_matrix(Xg.iloc[tr], yt, Xg.iloc[te], seed + r, nonlinear=True)
        pc = _fit_predict_matrix(Xc.iloc[tr], yt, Xc.iloc[te], seed + 10_000 + r, nonlinear=True)
        mz_auc = safe_auc(ye, zte)
        mg = probs_metrics(ye, pg)
        mc = probs_metrics(ye, pc)
        rows.append({
            "split": r,
            "n_train": len(tr),
            "n_test": len(te),
            "cases_test": int(ye.sum()),
            "zebra_auc": mz_auc,
            "genomic_auc": mg["auc"],
            "combined_auc": mc["auc"],
            "combined_minus_zebra_auc": mc["auc"] - mz_auc,
            "genomic_minus_zebra_auc": mg["auc"] - mz_auc,
            "genomic_average_precision": mg["average_precision"],
            "combined_average_precision": mc["average_precision"],
            "genomic_brier": mg["brier"],
            "combined_brier": mc["brier"],
            "genomic_log_loss": mg["log_loss"],
            "combined_log_loss": mc["log_loss"],
        })
        for name, test_score in [("ZeBRA", zte), ("Genomic", pg), ("Combined", pc)]:
            if name == "ZeBRA":
                train_score = ztr
            elif name == "Genomic":
                train_score = _fit_predict_matrix(Xg.iloc[tr], yt, Xg.iloc[tr], seed + r, nonlinear=True)
            else:
                train_score = _fit_predict_matrix(Xc.iloc[tr], yt, Xc.iloc[tr], seed + 10_000 + r, nonlinear=True)
            for fpr in fprs:
                op = _operating_point(yt, train_score, ye, test_score, float(fpr))
                ops.append({"split": r, "model": name, "requested_fpr": float(fpr), **op})
    rdf = pd.DataFrame(rows)
    rdf.to_csv(out / "auc_by_split.csv", index=False)
    summary = []
    for c in ["zebra_auc", "genomic_auc", "combined_auc", "combined_minus_zebra_auc", "genomic_minus_zebra_auc"]:
        summary.append({"metric": c, **metric_summary(rdf[c])})
    pd.DataFrame(summary).to_csv(out / "auc_distribution_summary.csv", index=False)
    opdf = pd.DataFrame(ops)
    opdf.to_csv(out / "operating_points_by_split.csv", index=False)
    op_summary = opdf.groupby(["model", "requested_fpr"], as_index=False).agg(
        sensitivity_mean=("test_sensitivity", "mean"),
        fpr_mean=("test_fpr", "mean"),
        lr_plus_mean=("lr_plus", "mean"),
        sensitivity_median=("test_sensitivity", "median"),
        fpr_median=("test_fpr", "median"),
    )
    op_summary.to_csv(out / "operating_points_summary.csv", index=False)


def run_incremental_logistic(d: pd.DataFrame, feature_cols: list[str], outroot: Path, params: dict[str, Any]) -> None:
    out = ensure_dir(outroot / "02_INCREMENTAL_LOGISTIC")
    nrep = int(params["incremental_repeats"])
    train_size = float(params["train_size"])
    seed = int(params["random_state"])
    y = d.target.to_numpy(int)
    Xz = d[["zebra_score"]]
    Xzg = pd.concat([Xz, d[feature_cols]], axis=1)
    splitter = StratifiedShuffleSplit(n_splits=nrep, train_size=train_size, random_state=seed + 101)
    rows = []
    for r, (tr, te) in enumerate(splitter.split(np.zeros(len(y)), y)):
        yt, ye = y[tr], y[te]
        pz = _fit_predict_matrix(Xz.iloc[tr], yt, Xz.iloc[te], seed + r, nonlinear=False)
        pc = _fit_predict_matrix(Xzg.iloc[tr], yt, Xzg.iloc[te], seed + 5000 + r, nonlinear=False)
        a = probs_metrics(ye, pz)
        b = probs_metrics(ye, pc)
        rows.append({
            "split": r,
            **{f"zebra_{k}": v for k, v in a.items()},
            **{f"combined_{k}": v for k, v in b.items()},
            "delta_auc": b["auc"] - a["auc"],
            "delta_average_precision": b["average_precision"] - a["average_precision"],
            "delta_brier": b["brier"] - a["brier"],
            "delta_log_loss": b["log_loss"] - a["log_loss"],
        })
    rdf = pd.DataFrame(rows)
    rdf.to_csv(out / "paired_incremental_by_split.csv", index=False)
    sm = []
    for c in ["delta_auc", "delta_average_precision", "delta_brier", "delta_log_loss"]:
        sm.append({"metric": c, **metric_summary(rdf[c])})
    pd.DataFrame(sm).to_csv(out / "paired_incremental_summary.csv", index=False)
