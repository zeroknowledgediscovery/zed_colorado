#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 genomicdataheader.csv [output.csv]" >&2
  exit 2
fi
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${2:-HFREF_genomic_header_matches.csv}"
python "$ROOT/scripts/find_phenotype_genomic_headers.py" "$1" --panel HFREF --output "$OUT"
