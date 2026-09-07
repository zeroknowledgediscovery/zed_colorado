#!/usr/bin/env python3
"""Execute the curated manuscript analyses using run_parameters.json.

The source notebooks/scripts are never modified. Parameter values are patched
into temporary execution copies, which are removed at the end of the run.
All analyses execute with curated_2 as the working directory and therefore
write their normal RESULTS/... outputs in place.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PARAM_FILE = ROOT / "run_parameters.json"

NOTEBOOK32 = ROOT / "32_COMBINED_CLASSIFIERS.ipynb"
NOTEBOOK33 = ROOT / "33_INCREMENTAL_LOGISTIC_ANALYSIS.ipynb"
APPARENT_SCRIPT = ROOT / "test_zebra_muc5b_apparent_association.py"
STRATIFIED_SCRIPT = ROOT / "testinverse_muc5b_stratified_by_fild.py"
LOCAL_SCRIPT = ROOT / "test_local_zebra_genomics_predictive_curves.py"
RESCUE_SCRIPT = ROOT / "test_muc5b_zebra_rescue_matched_fpr_lr.py"
HULL_SCRIPT = ROOT / "test_zebra_hybrid_roc_convex_hull_zedstat.py"

REQUIRED_INPUTS = [
    ROOT / "ILD_TOP_DRIVERS_DATA.csv",
    ROOT / "REPHENOTYPES FOR IC.csv",
    ROOT / "PREDICTIONS_104W_PRED_WINDOW.parquet",
]
MUC5B_COLUMNS = [
    "rs35705950.1_G_0",
    "rs35705950.1_G_1",
    "rs35705950.1_G_2",
]


def replace_assignment(text: str, name: str, value) -> tuple[str, int]:
    value_text = repr(value)
    pattern = rf"(?m)(\b{re.escape(name)}\s*=\s*)([-+0-9.eE]+)"
    return re.subn(pattern, rf"\g<1>{value_text}", text)


def patch_text(text: str, assignments: dict[str, object], label: str) -> str:
    for name, value in assignments.items():
        text, n = replace_assignment(text, name, value)
        if n == 0:
            raise RuntimeError(f"Could not find {name}=... in {label}")
    return text


def patch_notebook(src: Path, dst: Path, assignments: dict[str, object]) -> None:
    nb = json.loads(src.read_text())
    found = {k: 0 for k in assignments}
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", "")
        was_list = isinstance(source, list)
        text = "".join(source) if was_list else source
        for name, value in assignments.items():
            new_text, n = replace_assignment(text, name, value)
            if n:
                text = new_text
                found[name] += n
        cell["source"] = text.splitlines(keepends=True) if was_list else text
    missing = [name for name, n in found.items() if n == 0]
    if missing:
        raise RuntimeError(
            f"Could not patch {', '.join(missing)} in notebook {src.name}"
        )
    dst.write_text(json.dumps(nb, indent=1))


def patch_script(src: Path, dst: Path, assignments: dict[str, object]) -> None:
    dst.write_text(patch_text(src.read_text(), assignments, src.name))


def run(cmd: list[str], cwd: Path) -> None:
    print("\n+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def validate_inputs() -> None:
    missing = [p.name for p in REQUIRED_INPUTS if not p.exists()]
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))

    import pandas as pd

    columns = list(pd.read_csv(REQUIRED_INPUTS[0], nrows=0).columns)
    if "patient_id" not in columns:
        raise SystemExit(
            "ILD_TOP_DRIVERS_DATA.csv must contain the join key 'patient_id'."
        )
    missing_muc5b = [c for c in MUC5B_COLUMNS if c not in columns]
    if missing_muc5b:
        raise SystemExit(
            "ILD_TOP_DRIVERS_DATA.csv is missing required MUC5B one-hot columns: "
            + ", ".join(missing_muc5b)
        )
    if len(columns) <= 4:
        raise SystemExit(
            "ILD_TOP_DRIVERS_DATA.csv does not appear to contain a genomic feature matrix."
        )

    pred_cols = set(pd.read_parquet(REQUIRED_INPUTS[2]).columns)
    if not {"patient_id", "predicted_risk"}.issubset(pred_cols):
        raise SystemExit(
            "PREDICTIONS_104W_PRED_WINDOW.parquet must contain patient_id and predicted_risk."
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["full", "fast"], default="full")
    args = parser.parse_args()

    params_all = json.loads(PARAM_FILE.read_text())
    if args.mode not in params_all:
        raise SystemExit(f"Profile {args.mode!r} is not defined in {PARAM_FILE}")
    p = params_all[args.mode]

    validate_inputs()
    print("Run profile:", args.mode)
    print(json.dumps(p, indent=2))

    with tempfile.TemporaryDirectory(prefix=".curated2_run_", dir=ROOT) as td:
        tmp = Path(td)

        n32 = tmp / NOTEBOOK32.name
        patch_notebook(
            NOTEBOOK32,
            n32,
            {
                "N_REPEATS": p["notebook32_repeats"],
                "TRAIN_SIZE": p["train_size"],
                "INNER_CV_FOLDS": p["inner_cv_folds"],
                "N_PARAM_TRIALS": p["model_search_trials"],
            },
        )

        n33 = tmp / NOTEBOOK33.name
        patch_notebook(
            NOTEBOOK33,
            n33,
            {
                "N_REPEATS": p["notebook33_repeats"],
                "TRAIN_SIZE": p["train_size"],
                "INNER_CV_FOLDS": p["inner_cv_folds"],
            },
        )

        apparent = tmp / APPARENT_SCRIPT.name
        patch_script(
            APPARENT_SCRIPT,
            apparent,
            {
                "CV_FOLDS": p["apparent_association_cv_folds"],
                "CV_SEED": p["random_state"],
            },
        )

        stratified = tmp / STRATIFIED_SCRIPT.name
        patch_script(
            STRATIFIED_SCRIPT,
            stratified,
            {
                "N_PERMUTATIONS": p["stratified_permutations"],
                "PERMUTATION_SEED": p["random_state"],
                "CV_FOLDS": p["stratified_cv_folds"],
                "CV_SEED": p["random_state"],
            },
        )

        # This analysis is deterministic and has no random-state parameter.
        local = tmp / LOCAL_SCRIPT.name
        patch_script(LOCAL_SCRIPT, local, {})

        # The rescue script currently draws independent split seeds internally;
        # only its exposed repeat count and train fraction are configurable.
        rescue = tmp / RESCUE_SCRIPT.name
        patch_script(
            RESCUE_SCRIPT,
            rescue,
            {
                "N_REPEATS": p["rescue_repeats"],
                "TRAIN_SIZE": p["train_size"],
            },
        )

        # The hull analysis uses K-fold OOF predictions, not a train fraction.
        hull = tmp / HULL_SCRIPT.name
        patch_script(
            HULL_SCRIPT,
            hull,
            {
                "N_FOLDS": p["hull_folds"],
                "RANDOM_STATE": p["random_state"],
            },
        )

        staged = []
        try:
            sources = [n32, n33, apparent, stratified, local, rescue, hull]
            for src in sources:
                dst = ROOT / (".run_" + src.name)
                shutil.copy2(src, dst)
                staged.append(dst)

            for nb in staged[:2]:
                run(
                    [
                        "jupyter", "nbconvert", "--to", "notebook", "--execute",
                        "--inplace", "--ExecutePreprocessor.timeout=-1", nb.name,
                    ],
                    ROOT,
                )

            for script in staged[2:]:
                run([sys.executable, script.name], ROOT)
        finally:
            for path in staged:
                path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
