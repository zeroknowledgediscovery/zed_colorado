#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GEN="$ROOT/RESULTS/ild_muc5b"
REF="$ROOT/../zebra_comp/RESULTS"
CONFIG="$ROOT/configs/ild_muc5b.json"

usage() {
  cat <<'EOF'
Usage:
  ./ild_reference.sh
      Compare an existing generalized ILD/MUC5B run against frozen zebra_comp results.

  ./ild_reference.sh --run full
      Run the full generalized ILD/MUC5B pipeline, then compare.

  ./ild_reference.sh --run fast
      Run the fast smoke-test profile, then compare approximately.
      FAST mode is useful for code/path validation but is not a formal parity test.

Exit status:
  0  no parity-critical FAILs (PASS/WARN/INFO are allowed)
  1  at least one parity-critical FAIL
  2  usage/input error
EOF
}

if [[ $# -gt 0 ]]; then
  if [[ "$1" != "--run" ]]; then
    usage >&2
    exit 2
  fi
  MODE="${2:-full}"
  case "$MODE" in
    full|fast) ;;
    *) echo "Mode must be full or fast" >&2; exit 2 ;;
  esac
  if [[ $# -gt 2 ]]; then
    usage >&2
    exit 2
  fi
  echo "Running generalized ILD/MUC5B pipeline in '$MODE' mode..."
  "$ROOT/run_all.sh" "$CONFIG" "$MODE"
fi

if [[ ! -d "$GEN" ]]; then
  echo "Missing generalized results: $GEN" >&2
  echo "Run: ./ild_reference.sh --run full" >&2
  exit 2
fi
if [[ ! -d "$REF" ]]; then
  echo "Missing frozen reference results: $REF" >&2
  exit 2
fi

python - "$GEN" "$REF" <<'PY'
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GEN = Path(sys.argv[1])
REF = Path(sys.argv[2])
OUT = GEN / "ILD_REFERENCE_COMPARISON.tsv"

rows = []
critical_failures = 0


def add(section, metric, ref, new, pass_tol=None, warn_tol=None, critical=True, note=""):
    global critical_failures
    try:
        r = float(ref)
        n = float(new)
        delta = n - r
        ad = abs(delta)
        if pass_tol is None:
            status = "INFO"
        elif ad <= pass_tol:
            status = "PASS"
        elif warn_tol is not None and ad <= warn_tol:
            status = "WARN"
        else:
            status = "FAIL" if critical else "WARN"
        if status == "FAIL" and critical:
            critical_failures += 1
        rows.append({
            "section": section,
            "metric": metric,
            "reference": r,
            "generalized": n,
            "delta": delta,
            "status": status,
            "critical": critical,
            "note": note,
        })
    except Exception as e:
        if critical:
            critical_failures += 1
        rows.append({
            "section": section,
            "metric": metric,
            "reference": ref,
            "generalized": new,
            "delta": np.nan,
            "status": "FAIL" if critical else "WARN",
            "critical": critical,
            "note": f"comparison error: {e}; {note}",
        })


def add_exact(section, metric, ref, new, critical=True, note=""):
    global critical_failures
    ok = ref == new
    status = "PASS" if ok else ("FAIL" if critical else "WARN")
    if status == "FAIL" and critical:
        critical_failures += 1
    rows.append({
        "section": section,
        "metric": metric,
        "reference": ref,
        "generalized": new,
        "delta": (float(new) - float(ref)) if isinstance(ref, (int,float,np.number)) and isinstance(new, (int,float,np.number)) else np.nan,
        "status": status,
        "critical": critical,
        "note": note,
    })


def need(path: Path):
    global critical_failures
    if not path.exists():
        critical_failures += 1
        rows.append({
            "section": "FILES",
            "metric": str(path.relative_to(GEN if str(path).startswith(str(GEN)) else REF)),
            "reference": "present",
            "generalized": "missing",
            "delta": np.nan,
            "status": "FAIL",
            "critical": True,
            "note": "required comparison file missing",
        })
        return False
    return True


# ------------------------------------------------------------------
# Run mode / manifest
# ------------------------------------------------------------------
manifest_path = GEN / "RUN_MANIFEST.json"
if not need(manifest_path):
    mode = "unknown"
    manifest = {}
else:
    manifest = json.loads(manifest_path.read_text())
    mode = manifest.get("mode", "unknown")

print(f"Generalized run mode: {mode}")
if mode == "fast":
    print("NOTE: fast mode is a smoke test; strict numerical parity should be judged from a full run.\n")

# ------------------------------------------------------------------
# 1. Global ZeBRA discrimination.
# Genomic/combined models are informational because generalized model
# fitting intentionally does not reproduce notebook-32 search exactly.
# ------------------------------------------------------------------
glob_ref_p = REF / "032_REPEATED_SPLITS_03_TARGET_LOGIC" / "auc_distribution_summary.csv"
glob_new_p = GEN / "01_GLOBAL_COMPARISON" / "auc_distribution_summary.csv"
if need(glob_ref_p) and need(glob_new_p):
    rr = pd.read_csv(glob_ref_p)
    nn = pd.read_csv(glob_new_p)
    rr = rr.loc[rr.target_name == "FILD or FILA ADJUDICATED"].copy()
    ref_auc = {r.model: float(r["mean"]) for _, r in rr.iterrows()}
    new_auc = {r.metric: float(r["mean"]) for _, r in nn.iterrows()}
    add("GLOBAL", "ZeBRA mean AUC", ref_auc["zebra"], new_auc["zebra_auc"], .003 if mode == "full" else .010, .008 if mode == "full" else .025, True,
        "same score/cohort target; should be close")
    add("GLOBAL", "Genomic mean AUC", ref_auc["genomic"], new_auc["genomic_auc"], None, None, False,
        "INFO only: generalized nonlinear model is not notebook-32-identical")
    add("GLOBAL", "Combined mean AUC", ref_auc["combined"], new_auc["combined_auc"], None, None, False,
        "INFO only: generalized nonlinear model is not notebook-32-identical")

# ------------------------------------------------------------------
# 2. Paired incremental regularized analysis: parity-critical.
# ------------------------------------------------------------------
inc_ref_p = REF / "033_INCREMENTAL_LOGISTIC_32_COHORT" / "paired_incremental_summary.csv"
inc_new_p = GEN / "02_INCREMENTAL_LOGISTIC" / "paired_incremental_summary.csv"
if need(inc_ref_p) and need(inc_new_p):
    rr = pd.read_csv(inc_ref_p)
    nn = pd.read_csv(inc_new_p)
    rr = rr.loc[rr.target_name == "FILD or FILA ADJUDICATED"]
    ref_map = {r.metric: float(r.mean_delta) for _, r in rr.iterrows()}
    new_map = {str(r.metric).replace("delta_", ""): float(r["mean"]) for _, r in nn.iterrows()}
    add("INCREMENTAL", "mean delta AUC", ref_map["auc"], new_map["auc"], .003 if mode == "full" else .010, .008 if mode == "full" else .025, True)
    add("INCREMENTAL", "mean delta average precision", ref_map["average_precision"], new_map["average_precision"], .015 if mode == "full" else .035, .035 if mode == "full" else .070, True)
    add("INCREMENTAL", "mean delta Brier", ref_map["brier"], new_map["brier"], .0010 if mode == "full" else .0030, .0030 if mode == "full" else .0080, True)
    add("INCREMENTAL", "mean delta log loss", ref_map["log_loss"], new_map["log_loss"], .006 if mode == "full" else .015, .015 if mode == "full" else .035, True)

# ------------------------------------------------------------------
# 3. Local-information bands: compare the directly corresponding
# FPR-defined bands and their principal quantities.
# ------------------------------------------------------------------
loc_ref_p = REF / "LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES" / "FPR_DEFINED_LOCAL_INFORMATION_BANDS.csv"
loc_new_p = GEN / "04_LOCAL_INFORMATION" / "control_fpr_bands.csv"
if need(loc_ref_p) and need(loc_new_p):
    rr = pd.read_csv(loc_ref_p)
    nn = pd.read_csv(loc_new_p)
    # frozen labels correspond to (0.02,0.05), (0.05,0.10), (0.10,0.20), (0.20,0.40)
    mapping = {
        "2_to_5pct_FPR_band": (0.02, 0.05),
        "5_to_10pct_FPR_band": (0.05, 0.10),
        "10_to_20pct_FPR_band": (0.10, 0.20),
        "20_to_40pct_FPR_band": (0.20, 0.40),
    }
    for label, (lo, hi) in mapping.items():
        r = rr.loc[rr.band == label]
        n = nn.loc[np.isclose(nn.lower_fpr, lo) & np.isclose(nn.upper_fpr, hi)]
        if r.empty or n.empty:
            add_exact("LOCAL", f"{label} row present", True, False, True)
            continue
        r = r.iloc[0]; n = n.iloc[0]
        add("LOCAL", f"{label}: LRT G|Z", r.LRT_G_given_Z, n.LRT_G_given_Z, 2.0 if mode == "full" else 5.0, 5.0 if mode == "full" else 10.0, True)
        add("LOCAL", f"{label}: AUC ZeBRA", r.AUC_Z, n.AUC_Z, .035 if mode == "full" else .070, .080 if mode == "full" else .15, True)
        add("LOCAL", f"{label}: AUC gene", r.AUC_MUC5B, n.AUC_gene, .035 if mode == "full" else .070, .080 if mode == "full" else .15, True)

# ------------------------------------------------------------------
# 4. Direct-LR primary specification and robustness grid.
# This is the strongest parity check.
# ------------------------------------------------------------------
lr_ref_p = REF / "EMPIRICAL_ZEBRA_MUC5B_LR_ROBUSTNESS" / "LR_ESTIMATOR_SENSITIVITY_GRID.csv"
lr_new_p = GEN / "05_DIRECT_LR" / "lr_estimator_sensitivity_grid.csv"
if need(lr_ref_p) and need(lr_new_p):
    rr = pd.read_csv(lr_ref_p)
    nn = pd.read_csv(lr_new_p)
    r0 = rr.loc[(rr.n_knots == 5) & np.isclose(rr.logistic_C, 1.0)].iloc[0]
    n0 = nn.loc[(nn.n_knots == 5) & np.isclose(nn.logistic_C, 1.0)].iloc[0]
    add_exact("DIRECT_LR", "called-genotype N", int(r0.N), int(n0.N), True)
    add_exact("DIRECT_LR", "called-genotype cases", int(r0.cases), int(n0.cases), True)
    add("DIRECT_LR", "P0 mean LR_Z", r0.LR_Z_P0_mean, n0.LR_Z_P0_mean, 1e-8, 1e-5, True)
    add("DIRECT_LR", "P0 mean LR_ZG", r0.LR_ZG_P0_mean, n0.LR_ZG_P0_mean, .02 if mode == "full" else .05, .05 if mode == "full" else .10, True)
    for col, label in [
        ("b_delta_0p05", "b_hat_0.05"),
        ("b_delta_0p01", "b_hat_0.01"),
        ("b_delta_0p005", "b_hat_0.005"),
    ]:
        add("DIRECT_LR", label, r0[col], n0[col], .03 if mode == "full" else .08, .08 if mode == "full" else .15, True)

    # Across-grid range agreement: max absolute deviation for b_0.05 and b_0.01.
    merged = rr.merge(nn, on=["n_knots", "logistic_C"], suffixes=("_ref", "_new"))
    for col, label in [("b_delta_0p05", "grid max |delta b_0.05|"), ("b_delta_0p01", "grid max |delta b_0.01|")]:
        md = float(np.max(np.abs(merged[f"{col}_new"] - merged[f"{col}_ref"])))
        # Compare against zero directly.
        add("DIRECT_LR", label, 0.0, md, .04 if mode == "full" else .10, .10 if mode == "full" else .20, True)

# ------------------------------------------------------------------
# 5. Boundary rescue: compare the six corresponding threshold pairs.
# ------------------------------------------------------------------
res_ref_p = REF / "MUC5B_ZEBRA_RESCUE_RULE" / "REPEATED_SPLIT_MUC5B_RESCUE_SUMMARY.csv"
res_new_p = GEN / "06_BOUNDARY_RESCUE" / "rescue_summary.csv"
if need(res_ref_p) and need(res_new_p):
    rr = pd.read_csv(res_ref_p)
    nn = pd.read_csv(res_new_p)
    for _, n in nn.iterrows():
        main = float(n.upper_target_fpr)
        lower = float(n.lower_target_fpr)
        r = rr.loc[np.isclose(rr.main_target_fpr, main) & np.isclose(rr.lower_target_fpr, lower)]
        if r.empty:
            add_exact("RESCUE", f"pair {main:g}/{lower:g} present", True, False, True)
            continue
        r = r.iloc[0]
        tag = f"{100*main:g}%->{100*lower:g}%"
        add("RESCUE", f"{tag}: delta sensitivity", r.rescue_minus_matched_sensitivity_mean, n.rescue_minus_matched_sensitivity_mean,
            .020 if mode == "full" else .050, .050 if mode == "full" else .10, True)
        add("RESCUE", f"{tag}: delta FPR", r.rescue_minus_matched_FPR_mean, n.rescue_minus_matched_FPR_mean,
            .003 if mode == "full" else .008, .008 if mode == "full" else .015, True)

# ------------------------------------------------------------------
# Output
# ------------------------------------------------------------------
result = pd.DataFrame(rows)
result.to_csv(OUT, sep="\t", index=False)

# Pretty terminal table without requiring tabulate.
def fmt(v):
    if isinstance(v, (float, np.floating)):
        if not math.isfinite(float(v)):
            return str(v)
        return f"{float(v):.6g}"
    return str(v)

print("=" * 110)
print("ILD GENERALIZED PIPELINE vs FROZEN MANUSCRIPT REFERENCE")
print("=" * 110)
print(f"{'STATUS':<6} {'SECTION':<12} {'METRIC':<46} {'REFERENCE':>12} {'GENERAL':>12} {'DELTA':>12}")
print("-" * 110)
for _, r in result.iterrows():
    print(f"{r.status:<6} {str(r.section):<12} {str(r.metric)[:46]:<46} {fmt(r.reference):>12} {fmt(r.generalized):>12} {fmt(r.delta):>12}")
print("-" * 110)
counts = result.status.value_counts().to_dict()
print("Summary:", ", ".join(f"{k}={counts.get(k,0)}" for k in ["PASS", "WARN", "FAIL", "INFO"]))
print(f"Detailed TSV: {OUT}")

if mode == "fast":
    print("\nFAST MODE: use './ild_reference.sh --run full' before declaring numerical parity.")

if critical_failures:
    print(f"\nPARITY RESULT: FAIL ({critical_failures} parity-critical comparison(s) failed)")
    sys.exit(1)
else:
    print("\nPARITY RESULT: no parity-critical failures")
    sys.exit(0)
PY
