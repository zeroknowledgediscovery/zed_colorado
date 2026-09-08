#!/usr/bin/env python3
"""Generate the manuscript's exhaustive retained genomic-panel inventory.

The inventory is derived from the header of ILD_TOP_DRIVERS_DATA.csv, which is
the processed participant-by-feature matrix used by the broad genomic-only and
ZeBRA+genomics analyses. This script does not attempt to reconstruct the
historical biological_drivers.csv versus MORE_LOCI.csv membership because the
former source list is not retained in the current repository.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from pathlib import Path

COLUMN_RE = re.compile(r"^(rs\d+(?:\.\d+)?)_([ACGT])_([012])$")

KNOWN_MANUSCRIPT_LOCI = {
    "rs35705950": "MUC5B",
    "rs2736100": "TERT",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        type=Path,
        default=root / "ILD_TOP_DRIVERS_DATA.csv",
        help="Processed genomic matrix (header is sufficient).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root.parent / "tex" / "supplementary_genomic_panel_inventory.csv",
        help="Output CSV path.",
    )
    return parser.parse_args()


def read_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle))


def build_inventory(columns: list[str]) -> list[dict[str, object]]:
    grouped: OrderedDict[tuple[str, str, str], set[str]] = OrderedDict()

    for column in columns:
        if column in {"patient_id", "target"}:
            continue
        match = COLUMN_RE.match(column)
        if not match:
            raise ValueError(f"Unexpected genomic feature column: {column}")
        locus_token, allele_tag, state = match.groups()
        rsid = locus_token.split(".", 1)[0]
        key = (rsid, locus_token, allele_tag)
        grouped.setdefault(key, set()).add(state)

    rows: list[dict[str, object]] = []
    for index, ((rsid, locus_token, allele_tag), states) in enumerate(grouped.items(), start=1):
        if states != {"0", "1", "2"}:
            raise ValueError(
                f"{locus_token}_{allele_tag} does not have exactly states 0/1/2: {sorted(states)}"
            )
        role = (
            "broad panel + focal mechanistic T-carrier channel"
            if rsid == "rs35705950"
            else "broad array-derived genomic panel"
        )
        rows.append(
            {
                "variant_index": index,
                "rsid": rsid,
                "processed_locus_token": locus_token,
                "allele_tag": allele_tag,
                "known_gene_or_locus": KNOWN_MANUSCRIPT_LOCI.get(rsid, ""),
                "manuscript_role": role,
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    rows = build_inventory(read_header(args.matrix))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    n_loci = len(rows)
    n_one_hot = 3 * n_loci
    print(f"Wrote {n_loci} retained variant loci ({n_one_hot} one-hot columns) to {args.output}")


if __name__ == "__main__":
    main()
