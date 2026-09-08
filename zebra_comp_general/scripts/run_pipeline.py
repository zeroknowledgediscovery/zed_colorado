#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn

from analysis_global import run_global_comparison, run_incremental_logistic
from analysis_gene import run_gene_association
from analysis_local import run_local_information
from analysis_lr import run_direct_lr
from analysis_rescue import run_rescue
from common import build_analysis_frame, load_config, resolve_path

ROOT = Path(__file__).resolve().parents[1]
PARAM_FILE = ROOT / "run_parameters.json"


def main() -> None:
    ap = argparse.ArgumentParser(description="Run disease-agnostic ZeBRA + genomics comparison pipeline")
    ap.add_argument("--config", required=True, help="Disease/gene JSON config")
    ap.add_argument("--mode", choices=["full", "fast"], default="full")
    ap.add_argument("--keep-existing", action="store_true", help="Do not delete an existing analysis output directory")
    args = ap.parse_args()

    cfg = load_config(args.config)
    params_all = json.loads(PARAM_FILE.read_text())
    params = params_all[args.mode]

    analysis_name = cfg.get("analysis_name", Path(args.config).stem)
    result_base = resolve_path(cfg, cfg.get("results_root", str(ROOT / "RESULTS")))
    outroot = result_base / analysis_name
    if outroot.exists() and not args.keep_existing:
        shutil.rmtree(outroot)
    outroot.mkdir(parents=True, exist_ok=True)

    print(f"Analysis: {analysis_name}")
    print(f"Mode: {args.mode}")
    print(f"Output: {outroot}")
    print("Building standardized cohort...")
    d, feature_cols, meta = build_analysis_frame(cfg)
    if d.target.nunique() < 2:
        raise SystemExit("Primary phenotype contains only one class after filtering")
    if len(feature_cols) == 0:
        raise SystemExit("No genomic features available after selection/encoding")

    print(json.dumps(meta, indent=2))
    try:
        d.to_parquet(outroot / "analysis_frame.parquet", index=False)
    except Exception:
        d.to_csv(outroot / "analysis_frame.csv", index=False)
    pd.DataFrame({"genomic_feature": feature_cols}).to_csv(outroot / "selected_genomic_features.csv", index=False)

    manifest = {
        "analysis_name": analysis_name,
        "mode": args.mode,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": cfg["_config_path"],
        "config": {k: v for k, v in cfg.items() if not k.startswith("_")},
        "run_parameters": params,
        "cohort": meta,
        "versions": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
        },
    }
    (outroot / "RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2))

    run_global_comparison(d, feature_cols, outroot, params, cfg)
    run_incremental_logistic(d, feature_cols, outroot, params)
    run_gene_association(d, outroot, params)
    run_local_information(d, outroot, cfg)
    run_direct_lr(d, outroot, params, cfg)
    run_rescue(d, outroot, params, cfg)
    print(f"Completed: {outroot}")


if __name__ == "__main__":
    main()
