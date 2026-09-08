#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit

from analysis_global import ensure_dir
from analysis_gene import _called_gene


def _quantile_threshold(neg_scores: np.ndarray, fpr: float) -> float:
    return float(np.quantile(np.asarray(neg_scores, float), 1 - fpr, method="higher"))


def run_rescue(d: pd.DataFrame, outroot: Path, params: dict[str, Any], cfg: dict[str, Any]) -> None:
    out = ensure_dir(outroot / "06_BOUNDARY_RESCUE")
    x = _called_gene(d)
    if x.empty or x.gene_signal.nunique() < 2:
        (out / "SKIPPED.txt").write_text("No usable binary gene signal.\n")
        return
    aset = cfg.get("analysis", {})
    pairs = aset.get("rescue_fpr_pairs", [[0.005, 0.02], [0.005, 0.05], [0.005, 0.10], [0.01, 0.02], [0.01, 0.05], [0.01, 0.10]])
    reps = int(params["rescue_repeats"])
    train_size = float(params["train_size"])
    seed = int(params["random_state"])
    y = x.target.to_numpy(int)
    z = x.zebra_score.to_numpy(float)
    g = x.gene_signal.to_numpy(int)
    splitter = StratifiedShuffleSplit(n_splits=reps, train_size=train_size, random_state=seed + 707)
    rows = []
    for r, (tr, te) in enumerate(splitter.split(np.zeros(len(y)), 2 * y + g)):
        yt, ye = y[tr], y[te]
        zt, ze = z[tr], z[te]
        gt, ge = g[tr], g[te]
        neg = zt[yt == 0]
        for upper_fpr, lower_fpr in pairs:
            hi = _quantile_threshold(neg, float(upper_fpr))
            lo = _quantile_threshold(neg, float(lower_fpr))
            rescue_train = (zt >= hi) | ((zt >= lo) & (zt < hi) & (gt == 1))
            train_rescue_fpr = float(rescue_train[yt == 0].mean())
            matched_thr = _quantile_threshold(neg, train_rescue_fpr) if train_rescue_fpr > 0 else float(np.inf)
            rescue_test = (ze >= hi) | ((ze >= lo) & (ze < hi) & (ge == 1))
            matched_test = ze >= matched_thr
            cases = ye == 1
            controls = ye == 0
            rs = float(rescue_test[cases].mean()) if cases.any() else np.nan
            ms = float(matched_test[cases].mean()) if cases.any() else np.nan
            rf = float(rescue_test[controls].mean()) if controls.any() else np.nan
            mf = float(matched_test[controls].mean()) if controls.any() else np.nan
            rows.append({
                "split": r,
                "upper_target_fpr": float(upper_fpr),
                "lower_target_fpr": float(lower_fpr),
                "upper_threshold": hi,
                "lower_threshold": lo,
                "matched_threshold": matched_thr,
                "rescue_sensitivity": rs,
                "matched_sensitivity": ms,
                "rescue_minus_matched_sensitivity": rs - ms,
                "rescue_FPR": rf,
                "matched_FPR": mf,
                "rescue_minus_matched_FPR": rf - mf,
            })
    rdf = pd.DataFrame(rows)
    rdf.to_csv(out / "rescue_by_split.csv", index=False)
    summary = rdf.groupby(["upper_target_fpr", "lower_target_fpr"], as_index=False).agg(
        n_repeats=("split", "count"),
        rescue_sensitivity_mean=("rescue_sensitivity", "mean"),
        matched_sensitivity_mean=("matched_sensitivity", "mean"),
        rescue_minus_matched_sensitivity_mean=("rescue_minus_matched_sensitivity", "mean"),
        rescue_FPR_mean=("rescue_FPR", "mean"),
        matched_FPR_mean=("matched_FPR", "mean"),
        rescue_minus_matched_FPR_mean=("rescue_minus_matched_FPR", "mean"),
    )
    summary.to_csv(out / "rescue_summary.csv", index=False)
