#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler

from analysis_global import ensure_dir
from analysis_gene import _called_gene

# --- Direct score-level likelihood-ratio analysis ------------------------------
def _spline_model(knots: int, cval: float, seed: int):
    return Pipeline([
        ("spline", SplineTransformer(n_knots=knots, degree=3, knots="quantile", include_bias=False, extrapolation="constant")),
        ("scale", StandardScaler()),
        ("logistic", LogisticRegression(C=cval, solver="lbfgs", max_iter=5000, random_state=seed)),
    ])


def _fit_prob(xtr, ytr, xte, knots, cval, seed):
    xtr = np.asarray(xtr, float).reshape(-1, 1); xte = np.asarray(xte, float).reshape(-1, 1); ytr = np.asarray(ytr, int)
    if np.unique(ytr).size < 2:
        p = (ytr.sum() + .5) / (len(ytr) + 1.0); return np.full(len(xte), p)
    m = _spline_model(knots, cval, seed)
    with warnings.catch_warnings(): warnings.simplefilter("ignore"); m.fit(xtr, ytr)
    return m.predict_proba(xte)[:, 1]


def _lr_crossfit(x: pd.DataFrame, params: dict[str, Any]):
    y = x.target.to_numpy(int); g = x.gene_signal.to_numpy(int); zraw = x.zebra_score.to_numpy(float)
    # Logit only if bounded inside [0,1]; otherwise standardized rank-preserving raw score.
    if np.nanmin(zraw) >= 0 and np.nanmax(zraw) <= 1:
        s = np.clip(zraw, 1e-8, 1 - 1e-8); z = np.log(s / (1 - s))
    else:
        z = zraw.copy()
    strata = 2 * y + g
    folds = int(params["lr_folds"]); reps = int(params["lr_repeats"]); seed = int(params["random_state"])
    mincell = int(pd.Series(strata).value_counts().min()); folds = min(folds, mincell)
    if folds < 2:
        raise ValueError("A target x gene cell is too small for LR cross-fitting")
    knots = int(params["lr_spline_knots"]); cval = float(params["lr_logistic_c"])
    sy = np.zeros(len(x)); s1 = np.zeros(len(x)); s0 = np.zeros(len(x)); cnt = np.zeros(len(x), int)
    for r in range(reps):
        skf = StratifiedKFold(folds, shuffle=True, random_state=seed + r)
        for tr, te in skf.split(z, strata):
            yt, gt, zt = y[tr], g[tr], z[tr]
            sy[te] += _fit_prob(zt, yt, z[te], knots, cval, seed + r)
            s1[te] += _fit_prob(zt[yt == 1], gt[yt == 1], z[te], knots, cval, seed + 100 + r)
            s0[te] += _fit_prob(zt[yt == 0], gt[yt == 0], z[te], knots, cval, seed + 200 + r)
            cnt[te] += 1
    return sy / cnt, s1 / cnt, s0 / cnt


def _construct_lr(x: pd.DataFrame, py, q1, q0):
    eps = 1e-10
    out = x[["target", "zebra_score", "gene_signal"]].copy().reset_index(drop=True)
    out["score_percentile"] = out.zebra_score.rank(method="average", pct=True) * 100
    py, q1, q0 = [np.clip(np.asarray(a), eps, 1 - eps) for a in (py, q1, q0)]
    pi = float(out.target.mean()); prior = pi / (1 - pi)
    out["p_Y1_given_Z"] = py; out["p_G1_given_Y1_Z"] = q1; out["p_G1_given_Y0_Z"] = q0
    out["LR_Z_raw"] = (py / (1 - py)) / prior
    norm = float(out.loc[out.target == 0, "LR_Z_raw"].mean()); out["LR_Z"] = out.LR_Z_raw / norm
    out["LR_G1_given_Z"] = q1 / q0; out["LR_G0_given_Z"] = (1 - q1) / (1 - q0)
    out["LR_G_given_Z_observed"] = np.where(out.gene_signal == 1, out.LR_G1_given_Z, out.LR_G0_given_Z)
    out["LR_ZG"] = out.LR_Z * out.LR_G_given_Z_observed
    for c in ["LR_Z", "LR_G1_given_Z", "LR_G0_given_Z", "LR_G_given_Z_observed", "LR_ZG"]:
        out["log_" + c] = np.log(np.clip(out[c], 1e-300, None))
    out["max_abs_log_LR_G_statewise"] = np.maximum(np.abs(out.log_LR_G1_given_Z), np.abs(out.log_LR_G0_given_Z))
    return out


def run_direct_lr(d: pd.DataFrame, outroot: Path, params: dict[str, Any], cfg: dict[str, Any]) -> None:
    out = ensure_dir(outroot / "05_DIRECT_LR")
    x = _called_gene(d)
    if x.empty or x.gene_signal.nunique() < 2:
        (out / "SKIPPED.txt").write_text("No usable binary gene signal.\n"); return
    py, q1, q0 = _lr_crossfit(x, params); s = _construct_lr(x, py, q1, q0)
    try:
        s.to_parquet(out / "crossfitted_lr_rows.parquet", index=False)
    except Exception:
        s.to_csv(out / "crossfitted_lr_rows.csv", index=False)
    deltas = cfg.get("analysis", {}).get("lr_deltas", [0.10, 0.05, 0.025, 0.01, 0.005])
    controls = s.loc[s.target == 0]; vals = controls.max_abs_log_LR_G_statewise.to_numpy(float)
    st = []
    for delta in deltas:
        b = float(np.quantile(vals, 1 - float(delta)))
        st.append({"delta": float(delta), "b_delta_statewise_abs_logLR": b, "exp_b_delta": float(np.exp(b)), "empirical_tail_fraction": float(np.mean(vals > b))})
    pd.DataFrame(st).to_csv(out / "stochastic_bounds.csv", index=False)
    tails = []
    for pct in cfg.get("analysis", {}).get("tail_percentiles", [95, 97.5, 99]):
        w = s.loc[s.score_percentile >= float(pct)]
        tails.append({"score_percentile_cutoff": float(pct), "N": len(w), "cases": int(w.target.sum()),
                      "median_log_LR_Z": float(w.log_LR_Z.median()),
                      "q95_max_abs_log_LR_G_statewise": float(w.max_abs_log_LR_G_statewise.quantile(.95)),
                      "fraction_log_LR_Z_exceeds_statewise_genomic_scale": float((w.log_LR_Z > w.max_abs_log_LR_G_statewise).mean())})
    pd.DataFrame(tails).to_csv(out / "clinical_tail_summary.csv", index=False)

    # Decision flip localization at fixed clinical FPRs.
    flips = []
    for fpr in cfg.get("analysis", {}).get("operating_fprs", [0.005, 0.01, 0.02, 0.05, 0.10]):
        c = s.loc[s.target == 0]
        threshold = float(np.quantile(c.LR_Z.to_numpy(float), 1 - float(fpr), method="higher"))
        clinical = s.LR_Z >= threshold; joint = s.LR_ZG >= threshold
        changed = clinical != joint
        dist = np.abs(s.log_LR_Z - np.log(threshold))
        fd = dist.loc[changed]
        flips.append({"requested_clinical_FPR": float(fpr), "LR_threshold": threshold, "n_flips": int(changed.sum()),
                      "n_control_flips": int((changed & (s.target == 0)).sum()),
                      "median_boundary_distance_flips": float(fd.median()) if len(fd) else np.nan,
                      "q95_boundary_distance_flips": float(fd.quantile(.95)) if len(fd) else np.nan})
    pd.DataFrame(flips).to_csv(out / "decision_flip_localization.csv", index=False)

    # Prespecified estimator-sensitivity grid, mirroring the current ILD analysis.
    # This is a robustness analysis, not a model-selection exercise.
    aset = cfg.get("analysis", {})
    if aset.get("lr_robustness_grid", True):
        knot_grid = aset.get("lr_robustness_knots", [4, 5, 6])
        c_grid = aset.get("lr_robustness_c", [0.3, 1.0, 3.0])
        grid_rows = []
        tail_rows = []
        flip_rows = []
        for knots in knot_grid:
            for cval in c_grid:
                p2 = dict(params)
                p2["lr_spline_knots"] = int(knots)
                p2["lr_logistic_c"] = float(cval)
                py2, q12, q02 = _lr_crossfit(x, p2)
                s2 = _construct_lr(x, py2, q12, q02)
                c2 = s2.loc[s2.target == 0]
                vals2 = c2.max_abs_log_LR_G_statewise.to_numpy(float)
                row = {
                    "n_knots": int(knots),
                    "logistic_C": float(cval),
                    "N": len(s2),
                    "cases": int(s2.target.sum()),
                    "LR_Z_P0_mean": float(c2.LR_Z.mean()),
                    "LR_ZG_P0_mean": float(c2.LR_ZG.mean()),
                    "max_statewise_abs_logLR_G_controls": float(vals2.max()),
                }
                for delta in deltas:
                    tag = str(float(delta)).replace(".", "p")
                    b = float(np.quantile(vals2, 1 - float(delta)))
                    row[f"b_delta_{tag}"] = b
                    row[f"exp_b_delta_{tag}"] = float(np.exp(b))
                grid_rows.append(row)

                for pct in aset.get("tail_percentiles", [95, 97.5, 99]):
                    w = s2.loc[s2.score_percentile >= float(pct)]
                    tail_rows.append({
                        "n_knots": int(knots),
                        "logistic_C": float(cval),
                        "score_percentile_cutoff": float(pct),
                        "N": len(w),
                        "cases": int(w.target.sum()),
                        "median_log_LR_Z": float(w.log_LR_Z.median()),
                        "q95_max_abs_log_LR_G_statewise": float(w.max_abs_log_LR_G_statewise.quantile(.95)),
                        "fraction_log_LR_Z_exceeds_statewise_genomic_scale": float((w.log_LR_Z > w.max_abs_log_LR_G_statewise).mean()),
                    })

                for fpr in aset.get("operating_fprs", [0.005, 0.01, 0.02, 0.05, 0.10]):
                    threshold = float(np.quantile(c2.LR_Z.to_numpy(float), 1 - float(fpr), method="higher"))
                    changed = (s2.LR_Z >= threshold) != (s2.LR_ZG >= threshold)
                    dist = np.abs(s2.log_LR_Z - np.log(threshold)).loc[changed]
                    flip_rows.append({
                        "n_knots": int(knots),
                        "logistic_C": float(cval),
                        "requested_clinical_FPR": float(fpr),
                        "n_flips": int(changed.sum()),
                        "n_control_flips": int((changed & (s2.target == 0)).sum()),
                        "median_boundary_distance_flips": float(dist.median()) if len(dist) else np.nan,
                        "q95_boundary_distance_flips": float(dist.quantile(.95)) if len(dist) else np.nan,
                    })

        grid = pd.DataFrame(grid_rows)
        tails_grid = pd.DataFrame(tail_rows)
        flips_grid = pd.DataFrame(flip_rows)
        grid.to_csv(out / "lr_estimator_sensitivity_grid.csv", index=False)
        tails_grid.to_csv(out / "clinical_tail_robustness.csv", index=False)
        flips_grid.to_csv(out / "decision_flip_robustness.csv", index=False)

        ranges = []
        for delta in deltas:
            tag = str(float(delta)).replace(".", "p")
            col = f"b_delta_{tag}"
            v = grid[col].to_numpy(float)
            ranges.append({"metric": f"b_delta_{float(delta):g}", "min": float(v.min()), "median": float(np.median(v)), "max": float(v.max())})
        for pct in aset.get("tail_percentiles", [95, 97.5, 99]):
            q = tails_grid.loc[tails_grid.score_percentile_cutoff == float(pct)]
            for metric in ["median_log_LR_Z", "q95_max_abs_log_LR_G_statewise", "fraction_log_LR_Z_exceeds_statewise_genomic_scale"]:
                v = q[metric].to_numpy(float)
                ranges.append({"metric": f"tail_{float(pct):g}_{metric}", "min": float(v.min()), "median": float(np.median(v)), "max": float(v.max())})
        pd.DataFrame(ranges).to_csv(out / "robustness_range_summary.csv", index=False)
