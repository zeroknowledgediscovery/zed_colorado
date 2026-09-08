#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 CONFIG.json [full|fast]" >&2
  exit 2
fi

CONFIG="$1"
MODE="${2:-full}"
case "$MODE" in full|fast) ;; *) echo "Mode must be full or fast" >&2; exit 2;; esac

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
python scripts/run_pipeline.py --config "$CONFIG" --mode "$MODE"
