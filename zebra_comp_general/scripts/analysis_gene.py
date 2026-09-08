#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler

from analysis_global import ensure_dir, safe_auc


def _called_gene(d: pd.DataFrame) -> pd.DataFrame:
    x = d.loc[d.gene_signal.notna()].copy()
    x["gene_signal"] = x.gene_signal.astype(int)
    return x


def run_gene_association(d: pd.DataFrame, outroot: Path, params: dict[str, Any]) -> None:
    out = ensure_dir(outroot / "03_GENE_SCORE_ASSOCIATION")
    x = _called_gene(d)
    if x.empty or x.gene_signal.nunique() < 2:
        (out / "SKIPPED.txt").write_text("Gene-specific analyses skipped: no usable binary gene signal.\n")
        return
    folds = int(params["gene_cv_folds"])
    seed = int(params["random_state"])
    rows = []
    for label, w in [("pooled", x), ("target_0", x.loc[x.target == 0]), ("target_1", x.loc[x.target == 1])]:
        if len(w) < 10 or w.gene_signal.nunique() < 2:
            continue
        wz = w.zebra_score.to_numpy(float)
        wg = w.gene_signal.to_numpy(int)
        auc = safe_auc(wg, wz)
        cv_auc = np.nan
        min_class = int(pd.Series(wg).value_counts().min())
        k = min(folds, min_class)
        if k >= 2:
            cv = StratifiedKFold(k, shuffle=True, random_state=seed)
            model = make_pipeline(
                SplineTransformer(n_knots=5, degree=3, include_bias=False),
                StandardScaler(),
                LogisticRegression(C=1.0, max_iter=5000),
            )
            prob = cross_val_predict(model, wz.reshape(-1, 1), wg, cv=cv, method="predict_proba")[:, 1]
            cv_auc = safe_auc(wg, prob)
        rows.append({
            "stratum": label,
            "N": len(w),
            "gene_positive_N": int(wg.sum()),
            "gene_positive_prevalence": float(wg.mean()),
            "raw_score_to_gene_auc": auc,
            "direction_independent_auc": max(auc, 1 - auc),
            "spline_cv_auc": cv_auc,
        })
    pd.DataFrame(rows).to_csv(out / "gene_score_association_summary.csv", index=False)

    erows = []
    for label, w in [("pooled", x), ("target_0", x.loc[x.target == 0]), ("target_1", x.loc[x.target == 1])]:
        if len(w) < 10:
            continue
        for pct in [90, 95, 99]:
            q = float(w.zebra_score.quantile(pct / 100))
            top = w.loc[w.zebra_score >= q]
            base = float(w.gene_signal.mean())
            tail = float(top.gene_signal.mean()) if len(top) else np.nan
            erows.append({
                "stratum": label,
                "score_percentile": pct,
                "threshold": q,
                "N_tail": len(top),
                "base_gene_prevalence": base,
                "tail_gene_prevalence": tail,
                "enrichment": tail / base if base > 0 else np.nan,
            })
    pd.DataFrame(erows).to_csv(out / "gene_tail_enrichment.csv", index=False)
