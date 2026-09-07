#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-full}"
case "$MODE" in
  full|fast) ;;
  *)
    echo "Usage: $0 [full|fast]" >&2
    exit 2
    ;;
esac

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

for f in \
  ILD_TOP_DRIVERS_DATA.csv \
  "REPHENOTYPES FOR IC.csv" \
  PREDICTIONS_104W_PRED_WINDOW.parquet \
  run_parameters.json
do
  if [[ ! -f "$f" ]]; then
    echo "Missing required input: $f" >&2
    exit 1
  fi
done

echo "Removing previous RESULTS/ ..."
rm -rf RESULTS
mkdir -p RESULTS

echo "Running curated_2 in '$MODE' mode ..."
python scripts/run_with_config.py --mode "$MODE"

echo
echo "Completed. Results are under:"
find RESULTS -maxdepth 1 -mindepth 1 -type d -print | sort
