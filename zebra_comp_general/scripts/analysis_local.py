#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
from scipy.stats import chi2

from analysis_global import ensure_dir, safe_auc
from analysis_gene import _called_gene


def sigmoid(x):
    x = np.clip(x, -35, 35)
    return 1 / (1 + np.exp(-x))


def logistic_fit(X, y, ridge=1e-8, max_iter=200, tol=1e-9):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    p = X.shape[1]
    beta = np.zeros(p)
    penalty = np.eye(p)
    penalty[0, 0] = 0
    for _ in range(max_iter):
        prob = sigmoid(X @ beta)
        w = np.clip(prob * (1 - prob), 1e-8, None)
        grad = X.T @ (y - prob) - ridge * penalty @ beta
        h = X.T @ (X * w[:, None]) + ridge * penalty
        try:
            step = np.linalg.solve(h, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(h) @ grad
        new = beta + step
        if np.max(np.abs(new - beta)) < tol:
            beta = new
            break
        beta = new
    prob = sigmoid(X @ beta)
    w = np.clip(prob * (1 - prob), 1e-8, None)
    fisher = X.T @ (X * w[:, None]) + ridge * penalty
    cov = np.linalg.pinv(fisher)
    pp = np.clip(prob, 1e-12, 1 - 1e-12)
    ll = float(np.sum(y * np.log(pp) + (1 - y) * np.log(1 - pp)))
    return {"beta": beta, "cov": cov, "prob": prob, "ll": ll}


def lrt(ll_large, ll_small, df=1):
    stat = max(2 * (ll_large - ll_small), 0.0)
    return float(stat), float(chi2.sf(stat, df))


def coef_or(fit, idx):
    b = float(fit["beta"][idx])
    se = float(np.sqrt(max(fit["cov"][idx, idx], 0)))
    lo, hi = b - 1.95996398454 * se, b + 1.95996398454 * se
    return float(np.exp(b)), float(np.exp(lo)), float(np.exp(hi))


def _analyze_local_subset(w: pd.DataFrame) -> dict[str, Any] | None:
    if len(w) < 100 or w.target.sum() < 8 or (w.target == 0).sum() < 20 or w.gene_signal.nunique() < 2:
        return None
    y = w.target.to_numpy(int)
    zr = w.zebra_score.to_numpy(float)
    g = w.gene_signal.to_numpy(int)
    sd = float(np.std(zr))
    z = (zr - np.mean(zr)) / sd if sd > 1e-12 else np.zeros_like(zr)
    i = np.ones(len(y))
    fz = logistic_fit(np.column_stack([i, z]), y)
    fg = logistic_fit(np.column_stack([i, g]), y)
    fb = logistic_fit(np.column_stack([i, z, g]), y)
    lz, pz = lrt(fb["ll"], fg["ll"])
    lg, pg = lrt(fb["ll"], fz["ll"])
    oz, ozl, ozh = coef_or(fb, 1)
    og, ogl, ogh = coef_or(fb, 2)
    return {
        "N": len(w),
        "cases": int(y.sum()),
        "controls": int((y == 0).sum()),
        "case_prevalence": float(y.mean()),
        "LRT_Z_given_G": lz,
        "LRT_Z_given_G_p": pz,
        "LRT_G_given_Z": lg,
        "LRT_G_given_Z_p": pg,
        "G_minus_Z_partial_LRT": lg - lz,
        "OR_Z_per_local_SD": oz,
        "OR_Z_ci_low": ozl,
        "OR_Z_ci_high": ozh,
        "OR_gene_positive": og,
        "OR_gene_ci_low": ogl,
        "OR_gene_ci_high": ogh,
        "AUC_Z": safe_auc(y, zr),
        "AUC_gene": safe_auc(y, g),
        "AUC_both_in_sample": safe_auc(y, fb["prob"]),
        "gene_positive_prevalence": float(g.mean()),
        "case_prevalence_gene_positive": float(y[g == 1].mean()) if (g == 1).any() else np.nan,
        "case_prevalence_gene_negative": float(y[g == 0].mean()) if (g == 0).any() else np.nan,
    }


def run_local_information(d: pd.DataFrame, outroot: Path, cfg: dict[str, Any]) -> None:
    out = ensure_dir(outroot / "04_LOCAL_INFORMATION")
    x = _called_gene(d)
    if x.empty or x.gene_signal.nunique() < 2:
        (out / "SKIPPED.txt").write_text("No usable binary gene signal.\n")
        return
    x["score_percentile"] = x.zebra_score.rank(method="average", pct=True) * 100
    aset = cfg.get("analysis", {})
    half = float(aset.get("local_window_half_width_percentile", 5.0))
    centers = aset.get("local_window_centers", list(np.arange(25.0, 97.6, 2.5)))
    rows = []
    for c in centers:
        w = x.loc[x.score_percentile.between(float(c) - half, float(c) + half)]
        r = _analyze_local_subset(w)
        if r:
            rows.append({
                "window_center_percentile": float(c),
                "window_lower_percentile": float(c) - half,
                "window_upper_percentile": float(c) + half,
                **r,
            })
    pd.DataFrame(rows).to_csv(out / "local_score_windows.csv", index=False)

    bands = aset.get("control_fpr_bands", [[0.02, 0.05], [0.05, 0.10], [0.10, 0.20], [0.20, 0.40]])
    controls = x.loc[x.target == 0]
    brows = []
    for lower, upper in bands:
        high_thr = float(controls.zebra_score.quantile(1 - float(lower)))
        low_thr = float(controls.zebra_score.quantile(1 - float(upper)))
        w = x.loc[x.zebra_score.between(low_thr, high_thr, inclusive="both")]
        r = _analyze_local_subset(w)
        if r:
            brows.append({
                "lower_fpr": float(lower),
                "upper_fpr": float(upper),
                "lower_score": low_thr,
                "upper_score": high_thr,
                **r,
            })
    pd.DataFrame(brows).to_csv(out / "control_fpr_bands.csv", index=False)
