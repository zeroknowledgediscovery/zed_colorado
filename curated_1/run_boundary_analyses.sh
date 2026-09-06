#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

rm -rf \
  RESULTS/LOCAL_ZEBRA_MUC5B_INFORMATION_CURVES \
  RESULTS/MUC5B_ZEBRA_RESCUE_RULE \
  RESULTS/ZEBRA_HYBRID_ROC_CONVEX_HULL

python test_local_zebra_genomics_predictive_curves.py
python test_muc5b_zebra_rescue_matched_fpr_lr.py
python test_zebra_hybrid_roc_convex_hull_zedstat.py

echo "Regenerated boundary-analysis results under curated_1/RESULTS/."
