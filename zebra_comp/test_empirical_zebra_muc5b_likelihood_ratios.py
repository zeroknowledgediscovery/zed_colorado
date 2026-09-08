#!/usr/bin/env python3
"""Cross-fitted empirical LR analysis for ZeBRA and MUC5B.

Run from zebra_comp/. Uses the same three processed inputs as the other
manuscript scripts and writes aggregate/OOF results to
RESULTS/EMPIRICAL_ZEBRA_MUC5B_LIKELIHOOD_RATIOS/.

This estimates score-level quantities for Z=ZeBRA and G=MUC5B T-carrier:
  Lambda_Z(z)      = p(z|Y=1)/p(z|Y=0)
  Lambda_G|Z(g;z)  = p(g|Y=1,z)/p(g|Y=0,z)
so Lambda_Z,G = Lambda_Z * Lambda_G|Z exactly at the observable-score level.
The outputs assess empirical scale/tail behavior and boundary localization;
a finite observed range must NOT be interpreted as proof of a population
uniform bound.
"""
from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler

BASE_DATA_FILE = "ILD_TOP_DRIVERS_DATA.csv"
TARGET_FILE = "REPHENOTYPES FOR IC.csv"
PRED_FILE = "PREDICTIONS_104W_PRED_WINDOW.parquet"
TARGET_NAME = "FILD or FILA ADJUDICATED"
GG_COLUMN, GT_COLUMN, TT_COLUMN = "rs35705950.1_G_2", "rs35705950.1_G_1", "rs35705950.1_G_0"
OUTDIR = Path("./RESULTS/EMPIRICAL_ZEBRA_MUC5B_LIKELIHOOD_RATIOS"); OUTDIR.mkdir(parents=True, exist_ok=True)
SEED, N_SPLITS, N_REPEATS = 1321, 5, 5
N_KNOTS, DEGREE, LOGISTIC_C, MAX_ITER = 5, 3, 1.0, 5000
EPS = 1e-10
DELTAS = [0.10, 0.05, 0.025, 0.01, 0.005]
BOOTSTRAPS = 1000
FPRS = [0.005, 0.01, 0.02, 0.05, 0.10]
WINDOW_CENTERS = np.arange(25.0, 97.6, 2.5)
WINDOW_HALF_WIDTH = 5.0


def numeric_target(s):
    return pd.to_numeric(s.replace({"Y":1,"y":1,"N":0,"n":0,True:1,False:0}), errors="coerce").fillna(0).astype(int)


def check_unique(df, name):
    if "patient_id" not in df: raise ValueError(f"{name} lacks patient_id")
    ndup = int(df.patient_id.duplicated().sum())
    if ndup: raise ValueError(f"{name} has {ndup} duplicated patient_id values")


def decode_muc5b(d):
    req = [GG_COLUMN, GT_COLUMN, TT_COLUMN]
    miss = [c for c in req if c not in d]
    if miss: raise ValueError("Missing MUC5B columns: " + ", ".join(miss))
    d = d.copy()
    for c in req: d[c] = pd.to_numeric(d[c], errors="coerce")
    sm = d[req].sum(axis=1)
    if (sm > 1).any(): raise ValueError("Some MUC5B rows have >1 active genotype state")
    d = d.loc[sm == 1].copy()  # all-zero is missing genotype in historical one-hot encoding
    d["MUC5B_T_carrier"] = ((d[GT_COLUMN] == 1) | (d[TT_COLUMN] == 1)).astype(int)
    return d


def collapse_prediction_rows(pred, relevant_patient_ids):
    """Return one ZeBRA prediction per relevant patient.

    The historical prediction parquet can contain duplicate patient_id rows.
    Existing manuscript scripts merge that file directly, but repeated rows are
    inappropriate for cross-fitting because the same patient could otherwise
    appear in both training and test folds. We therefore collapse duplicates
    only after restricting to patients in the genomic analysis cohort.

    Collapsing is allowed only when all non-missing predicted_risk values for a
    patient are identical (to numerical tolerance). If a relevant patient has
    genuinely different risk values, the script stops rather than choosing an
    arbitrary aggregation rule.
    """
    if not {"patient_id", "predicted_risk"}.issubset(pred.columns):
        raise ValueError("Prediction parquet needs patient_id,predicted_risk")

    p = pred[["patient_id", "predicted_risk"]].copy()
    p = p.loc[p.patient_id.isin(set(relevant_patient_ids))].copy()
    p["predicted_risk"] = pd.to_numeric(p.predicted_risk, errors="coerce")

    n_rows = len(p)
    n_ids = p.patient_id.nunique(dropna=False)
    n_extra = n_rows - n_ids
    if n_extra:
        print(
            f"Prediction parquet has {n_extra:,} extra duplicate rows among "
            f"{n_ids:,} genomic-cohort patient IDs; checking whether risks agree..."
        )

    conflict_ids = []
    for pid, s in p.groupby("patient_id", sort=False)["predicted_risk"]:
        vals = s.dropna().to_numpy(float)
        if vals.size > 1 and not np.allclose(vals, vals[0], rtol=0.0, atol=1e-12):
            conflict_ids.append(pid)
            if len(conflict_ids) >= 20:
                break

    if conflict_ids:
        examples = ", ".join(map(str, conflict_ids[:10]))
        raise ValueError(
            "Relevant patients have multiple different predicted_risk values in "
            f"{PRED_FILE}; refusing arbitrary aggregation. Example patient_id(s): {examples}. "
            "Inspect the prediction parquet and define the intended 104-week row selection."
        )

    # If duplicates are exact (or differ only by missing copies), retain the
    # single non-missing value. This preserves one independent row per patient.
    def first_nonmissing_or_nan(s):
        x = s.dropna()
        return float(x.iloc[0]) if len(x) else np.nan

    p = (
        p.groupby("patient_id", as_index=False, sort=False)["predicted_risk"]
         .agg(first_nonmissing_or_nan)
    )

    check_unique(p, "collapsed prediction table")
    return p


def load_data():
    base = pd.read_csv(BASE_DATA_FILE)
    if "target" in base: base = base.drop(columns="target")
    check_unique(base, BASE_DATA_FILE)
    tar = pd.read_csv(TARGET_FILE)
    if "arb_person_id" in tar and "patient_id" not in tar: tar = tar.rename(columns={"arb_person_id":"patient_id"})
    check_unique(tar, TARGET_FILE)
    if TARGET_NAME not in tar: raise ValueError(f"Missing target: {TARGET_NAME}")

    pred_raw = pd.read_parquet(PRED_FILE)
    pred = collapse_prediction_rows(pred_raw, base.patient_id)

    tar = tar[["patient_id",TARGET_NAME]].copy(); tar["target"] = numeric_target(tar[TARGET_NAME])
    d = base.merge(tar[["patient_id","target"]], on="patient_id", how="left").merge(pred, on="patient_id", how="left")
    d["target"] = pd.to_numeric(d.target, errors="coerce").fillna(0).astype(int)
    d["predicted_risk"] = pd.to_numeric(d.predicted_risk, errors="coerce")
    d = d.loc[d.predicted_risk.notna()].copy(); d = decode_muc5b(d)
    check_unique(d, "final LR analysis cohort")
    if not d.predicted_risk.between(0,1).all(): raise ValueError("predicted_risk must lie in [0,1]")
    s = np.clip(d.predicted_risk.to_numpy(float), 1e-8, 1-1e-8)
    d["z_model"] = np.log(s/(1-s))  # monotone expansion only; not a calibration claim
    d["ZeBRA_percentile"] = d.predicted_risk.rank(method="average", pct=True)*100
    return d.reset_index(drop=True)


def model():
    return Pipeline([
        ("spline", SplineTransformer(n_knots=N_KNOTS, degree=DEGREE, knots="quantile", include_bias=False, extrapolation="constant")),
        ("scale", StandardScaler()),
        ("logistic", LogisticRegression(C=LOGISTIC_C, solver="lbfgs", max_iter=MAX_ITER, random_state=SEED)),
    ])


def fit_prob(xtr, ytr, xte):
    xtr, xte, ytr = np.asarray(xtr,float).reshape(-1,1), np.asarray(xte,float).reshape(-1,1), np.asarray(ytr,int)
    if np.unique(ytr).size < 2:
        p = (ytr.sum()+0.5)/(len(ytr)+1.0); return np.full(len(xte), p)
    m = model()
    with warnings.catch_warnings(): warnings.simplefilter("ignore"); m.fit(xtr,ytr)
    return m.predict_proba(xte)[:,1]


def crossfit(d):
    y, g, z = d.target.to_numpy(int), d.MUC5B_T_carrier.to_numpy(int), d.z_model.to_numpy(float)
    strata = 2*y+g
    if pd.Series(strata).value_counts().min() < N_SPLITS: raise ValueError("A target x genotype cell is too small for 5-fold cross-fitting")
    sy = np.zeros(len(d)); s1 = np.zeros(len(d)); s0 = np.zeros(len(d)); cnt = np.zeros(len(d),int); splitrows=[]
    for r in range(N_REPEATS):
        skf = StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED+r)
        for f,(tr,te) in enumerate(skf.split(z,strata)):
            yt,gt,zt = y[tr],g[tr],z[tr]
            sy[te] += fit_prob(zt,yt,z[te])
            s1[te] += fit_prob(zt[yt==1],gt[yt==1],z[te])
            s0[te] += fit_prob(zt[yt==0],gt[yt==0],z[te]); cnt[te] += 1
            splitrows.append({"repeat":r,"fold":f,"train_n":len(tr),"test_n":len(te),"train_cases":int(yt.sum()),"train_carriers":int(gt.sum())})
    if not np.all(cnt==N_REPEATS): raise RuntimeError("Incomplete cross-fitting")
    return sy/cnt, s1/cnt, s0/cnt, pd.DataFrame(splitrows)


def construct_lr(d, py, q1, q0):
    out = d[["target","predicted_risk","z_model","ZeBRA_percentile","MUC5B_T_carrier"]].copy()
    out.insert(0,"analysis_row",np.arange(len(out)))
    raw_py,raw_q1,raw_q0 = np.asarray(py),np.asarray(q1),np.asarray(q0)
    py,q1,q0 = [np.clip(x,EPS,1-EPS) for x in (raw_py,raw_q1,raw_q0)]
    out["p_Y1_given_Z"],out["p_G1_given_Y1_Z"],out["p_G1_given_Y0_Z"] = py,q1,q0
    pi = out.target.mean(); prior_odds = pi/(1-pi)
    out["LR_Z_raw"] = (py/(1-py))/prior_odds
    norm = out.loc[out.target==0,"LR_Z_raw"].mean(); out["LR_Z"] = out.LR_Z_raw/norm
    out["LR_G1_given_Z"] = q1/q0
    out["LR_G0_given_Z"] = (1-q1)/(1-q0)
    out["LR_G_given_Z_observed"] = np.where(out.MUC5B_T_carrier==1,out.LR_G1_given_Z,out.LR_G0_given_Z)
    out["LR_ZG"] = out.LR_Z*out.LR_G_given_Z_observed
    for c in ["LR_Z_raw","LR_Z","LR_G1_given_Z","LR_G0_given_Z","LR_G_given_Z_observed","LR_ZG"]: out["log_"+c] = np.log(np.clip(out[c],1e-300,None))
    out["max_abs_log_LR_G_statewise"] = np.maximum(np.abs(out.log_LR_G1_given_Z),np.abs(out.log_LR_G0_given_Z))
    diag = {
        "n":len(out),"cases":int(out.target.sum()),"controls":int((out.target==0).sum()),"carriers":int(out.MUC5B_T_carrier.sum()),
        "prevalence":float(pi),"LR_Z_raw_P0_mean_before_normalization":float(norm),"LR_Z_P0_mean_after_normalization":float(out.loc[out.target==0,"LR_Z"].mean()),
        "LR_ZG_P0_mean":float(out.loc[out.target==0,"LR_ZG"].mean()),"probability_clip_epsilon":EPS,
        "n_pY_clipped":int(((raw_py<=EPS)|(raw_py>=1-EPS)).sum()),"n_q1_clipped":int(((raw_q1<=EPS)|(raw_q1>=1-EPS)).sum()),"n_q0_clipped":int(((raw_q0<=EPS)|(raw_q0>=1-EPS)).sum()),
        "n_splits":N_SPLITS,"n_repeats":N_REPEATS,"n_knots":N_KNOTS,"degree":DEGREE,"logistic_C":LOGISTIC_C,"seed":SEED,
    }
    return out,diag


def bootstrap_qci(x,q,rng):
    x=np.asarray(x,float); n=len(x); v=np.empty(BOOTSTRAPS)
    for b in range(BOOTSTRAPS): v[b]=np.quantile(x[rng.integers(0,n,n)],q)
    return np.quantile(v,[.025,.975])


def stochastic_summary(s):
    c=s.loc[s.target==0]; x=c.max_abs_log_LR_G_statewise.to_numpy(float); allstates=np.r_[c.log_LR_G0_given_Z,c.log_LR_G1_given_Z]; rng=np.random.default_rng(SEED+991); rows=[]
    for delta in DELTAS:
        b=float(np.quantile(x,1-delta)); lo,hi=bootstrap_qci(x,1-delta,rng); qlo,qhi=np.quantile(allstates,[delta/2,1-delta/2])
        rows.append({"delta":delta,"b_delta_statewise_abs_logLR":b,"b_delta_bootstrap_ci_low":lo,"b_delta_bootstrap_ci_high":hi,"exp_b_delta":np.exp(b),"empirical_tail_fraction":float(np.mean(x>b)),"central_state_logLR_lower":qlo,"central_state_logLR_upper":qhi,"central_state_m_delta":np.exp(qlo),"central_state_M_delta":np.exp(qhi)})
    return pd.DataFrame(rows)


def marginal_muc5b(s):
    rows=[]
    for state in [0,1]:
        vals=[]
        for y in [1,0]:
            g=s.loc[s.target==y,"MUC5B_T_carrier"].to_numpy(int); vals.append(((g==state).sum()+.5)/(len(g)+1.0))
        rows.append({"MUC5B_T_carrier_state":state,"P_G_given_Y1":vals[0],"P_G_given_Y0":vals[1],"marginal_LR_G":vals[0]/vals[1],"log_marginal_LR_G":np.log(vals[0]/vals[1])})
    return pd.DataFrame(rows)


def moving_windows(s):
    rows=[]
    for center in WINDOW_CENTERS:
        w=s.loc[s.ZeBRA_percentile.between(center-WINDOW_HALF_WIDTH,center+WINDOW_HALF_WIDTH)]
        if w.empty: continue
        rows.append({"window_center_percentile":center,"window_lower_percentile":center-WINDOW_HALF_WIDTH,"window_upper_percentile":center+WINDOW_HALF_WIDTH,"N":len(w),"cases":int(w.target.sum()),
                     "median_log_LR_Z":w.log_LR_Z.median(),"q10_log_LR_Z":w.log_LR_Z.quantile(.1),"q90_log_LR_Z":w.log_LR_Z.quantile(.9),
                     "median_log_LR_G1_given_Z":w.log_LR_G1_given_Z.median(),"median_log_LR_G0_given_Z":w.log_LR_G0_given_Z.median(),
                     "median_max_abs_log_LR_G_statewise":w.max_abs_log_LR_G_statewise.median(),"q95_max_abs_log_LR_G_statewise":w.max_abs_log_LR_G_statewise.quantile(.95),
                     "median_log_LR_G_observed":w.log_LR_G_given_Z_observed.median(),"median_log_LR_ZG":w.log_LR_ZG.median(),
                     "fraction_log_LR_Z_exceeds_statewise_genomic_scale":float(np.mean(w.log_LR_Z>w.max_abs_log_LR_G_statewise))})
    return pd.DataFrame(rows)


def clinical_tail(s):
    rows=[]
    for cutoff in [80,90,95,97.5,99]:
        w=s.loc[s.ZeBRA_percentile>=cutoff]
        rows.append({"ZeBRA_percentile_cutoff":cutoff,"N":len(w),"cases":int(w.target.sum()),"controls":int((w.target==0).sum()),"median_log_LR_Z":w.log_LR_Z.median(),"q90_log_LR_Z":w.log_LR_Z.quantile(.9),"max_log_LR_Z":w.log_LR_Z.max(),"median_max_abs_log_LR_G_statewise":w.max_abs_log_LR_G_statewise.median(),"q95_max_abs_log_LR_G_statewise":w.max_abs_log_LR_G_statewise.quantile(.95),"max_max_abs_log_LR_G_statewise":w.max_abs_log_LR_G_statewise.max(),"fraction_log_LR_Z_exceeds_statewise_genomic_scale":float(np.mean(w.log_LR_Z>w.max_abs_log_LR_G_statewise))})
    return pd.DataFrame(rows)


def qhigher(x,q):
    try: return float(np.quantile(x,q,method="higher"))
    except TypeError: return float(np.quantile(x,q,interpolation="higher"))


def decision_localization(s,stoch):
    y=s.target.to_numpy(int); lz=s.log_LR_Z.to_numpy(float); lg=s.log_LR_G_given_Z_observed.to_numpy(float); lj=s.log_LR_ZG.to_numpy(float); controls=y==0; bmap=dict(zip(stoch.delta,stoch.b_delta_statewise_abs_logLR)); rows=[]; events=[]
    for fpr in FPRS:
        lam=qhigher(s.loc[controls,"LR_Z"],1-fpr); ll=np.log(lam); cp=lz>=ll; jp=lj>=ll; flip=cp!=jp; up=(~cp)&jp; down=cp&(~jp); dist=np.abs(lz-ll); fd=dist[flip]
        row={"requested_clinical_FPR":fpr,"lambda_LR_Z":lam,"log_lambda_LR_Z":ll,"clinical_FPR":float(cp[controls].mean()),"joint_same_lambda_FPR":float(jp[controls].mean()),"clinical_sensitivity":float(cp[y==1].mean()),"joint_same_lambda_sensitivity":float(jp[y==1].mean()),"n_flips":int(flip.sum()),"n_upward_flips":int(up.sum()),"n_downward_flips":int(down.sum()),"overall_flip_fraction":float(flip.mean()),"control_flip_fraction":float(flip[controls].mean()),"case_flip_fraction":float(flip[y==1].mean())}
        for name,q in [("median",.5),("q90",.9),("q95",.95)]: row[f"{name}_boundary_distance_flips"] = float(np.quantile(fd,q)) if len(fd) else np.nan
        row["max_boundary_distance_flips"] = float(fd.max()) if len(fd) else np.nan
        row["max_flip_distance_minus_abs_observed_logLR_G"] = float(np.max(dist[flip]-np.abs(lg[flip]))) if len(fd) else np.nan
        for delta,b in bmap.items(): row[f"flip_fraction_outside_b_delta_{delta:g}"] = float(np.mean(fd>b)) if len(fd) else np.nan
        rows.append(row)
        for i in np.flatnonzero(flip): events.append({"analysis_row":int(s.iloc[i].analysis_row),"requested_clinical_FPR":fpr,"target":int(y[i]),"MUC5B_T_carrier":int(s.iloc[i].MUC5B_T_carrier),"ZeBRA_percentile":float(s.iloc[i].ZeBRA_percentile),"log_LR_Z":lz[i],"log_LR_G_given_Z_observed":lg[i],"log_LR_ZG":lj[i],"log_lambda":ll,"boundary_distance":dist[i],"flip_direction":"upward" if up[i] else "downward"})
    return pd.DataFrame(rows),pd.DataFrame(events)


def plots(s,w,stoch,dec):
    fig,ax=plt.subplots(figsize=(9.5,5.8)); ax.plot(w.window_center_percentile,w.median_log_LR_Z,marker="o",label="log Lambda_Z"); ax.plot(w.window_center_percentile,w.median_log_LR_G1_given_Z,marker="o",label="carrier log Lambda_G|Z"); ax.plot(w.window_center_percentile,w.median_log_LR_G0_given_Z,marker="o",label="non-carrier log Lambda_G|Z"); ax.axhline(0,ls="--",lw=1); ax.set(xlabel="ZeBRA risk percentile",ylabel="Estimated log likelihood ratio",title="Cross-fitted clinical and residual MUC5B likelihood ratios"); ax.legend(); fig.tight_layout(); fig.savefig(OUTDIR/"EMPIRICAL_LOG_LR_CHANNELS.png",dpi=220); plt.close(fig)
    c=s.loc[s.target==0]; x=np.sort(c.max_abs_log_LR_G_statewise); surv=np.maximum(1-np.arange(1,len(x)+1)/len(x),1/len(x)); fig,ax=plt.subplots(figsize=(8.4,5.5)); ax.step(x,surv,where="post"); ax.set_yscale("log"); ax.set(xlabel="b: statewise max |log Lambda_G|Z|",ylabel="Empirical P0(tail > b)",title="Residual-genomic empirical tail under operational controls"); fig.tight_layout(); fig.savefig(OUTDIR/"EMPIRICAL_RESIDUAL_GENOMIC_TAIL.png",dpi=220); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.4,5.5)); x=100*dec.requested_clinical_FPR; ax.plot(x,dec.median_boundary_distance_flips,marker="o",label="median"); ax.plot(x,dec.q95_boundary_distance_flips,marker="o",label="95th percentile"); ax.set_xscale("log"); ax.set_xticks([.5,1,2,5,10]); ax.set_xticklabels(["0.5","1","2","5","10"]); ax.set(xlabel="Requested clinical FPR (%)",ylabel="|log Lambda_Z - log lambda| among flips",title="Localization of MUC5B-induced fixed-LR decision changes"); ax.legend(); fig.tight_layout(); fig.savefig(OUTDIR/"EMPIRICAL_DECISION_FLIP_LOCALIZATION.png",dpi=220); plt.close(fig)


def main():
    d=load_data(); print(f"N={len(d):,}, cases={int(d.target.sum()):,}, MUC5B carriers={int(d.MUC5B_T_carrier.sum()):,}")
    py,q1,q0,splits=crossfit(d); s,diag=construct_lr(d,py,q1,q0); stoch=stochastic_summary(s); marg=marginal_muc5b(s); win=moving_windows(s); tail=clinical_tail(s); dec,events=decision_localization(s,stoch)
    s.to_csv(OUTDIR/"EMPIRICAL_LR_SUBJECT_LEVEL.csv",index=False); splits.to_csv(OUTDIR/"CROSSFIT_SPLIT_SUMMARY.csv",index=False); stoch.to_csv(OUTDIR/"STOCHASTIC_BOUND_SUMMARY.csv",index=False); marg.to_csv(OUTDIR/"MARGINAL_MUC5B_LR.csv",index=False); win.to_csv(OUTDIR/"EMPIRICAL_LR_MOVING_WINDOWS.csv",index=False); tail.to_csv(OUTDIR/"CLINICAL_TAIL_SUMMARY.csv",index=False); dec.to_csv(OUTDIR/"DECISION_FLIP_LOCALIZATION_SUMMARY.csv",index=False); events.to_csv(OUTDIR/"DECISION_FLIP_EVENTS.csv",index=False)
    (OUTDIR/"MODEL_DIAGNOSTICS.json").write_text(json.dumps(diag,indent=2,sort_keys=True)+"\n")
    plots(s,win,stoch,dec)
    lines=["EMPIRICAL ZEBRA-MUC5B LIKELIHOOD-RATIO ANALYSIS","Finite observed ranges are descriptive, not proof of a uniform population bound.",f"N={diag['n']}; cases={diag['cases']}; controls={diag['controls']}; carriers={diag['carriers']}",f"P0 mean LR_Z after normalization={diag['LR_Z_P0_mean_after_normalization']:.6f}; P0 mean LR_ZG={diag['LR_ZG_P0_mean']:.6f}","","Empirical stochastic residual-genomic scale under P0:"]
    for _,r in stoch.iterrows(): lines.append(f"  delta={r.delta:g}: b_delta={r.b_delta_statewise_abs_logLR:.4f}, exp(b)={r.exp_b_delta:.3f}, bootstrap95=[{r.b_delta_bootstrap_ci_low:.4f},{r.b_delta_bootstrap_ci_high:.4f}]")
    lines.append(""); lines.append("Fixed-LR decision localization:")
    for _,r in dec.iterrows(): lines.append(f"  FPR={100*r.requested_clinical_FPR:.1f}%: flips={int(r.n_flips)}, median distance={r.median_boundary_distance_flips:.4f}, q95={r.q95_boundary_distance_flips:.4f}")
    text="\n".join(lines)+"\n"; (OUTDIR/"SUMMARY.txt").write_text(text); print(text); print("Wrote",OUTDIR.resolve())

if __name__ == "__main__": main()
