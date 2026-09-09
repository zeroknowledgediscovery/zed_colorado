#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 /path/to/genomicdataheader.csv /path/to/genomicdata.csv [ID_COLUMN]" >&2
  exit 2
fi

HEADER="$1"
GENOMIC="$2"
IDCOL="${3:-FID}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python "$HERE/01_find_apoe_columns.py" "$HEADER" --output APOE_candidates.csv
python "$HERE/02_inspect_apoe_columns.py" "$GENOMIC" APOE_candidates.csv --id-column "$IDCOL"

echo
echo "Done with discovery/inspection."
echo "Now inspect APOE_candidates.csv and APOE_column_value_counts.csv."
echo "Then run 03_make_apoe4_carrier.py with the exact two column names and dosage alleles."
