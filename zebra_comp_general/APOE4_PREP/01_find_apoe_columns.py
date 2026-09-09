#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

TARGETS = {
    "rs429358": {"chrom": "19", "pos_grch37": 45411941, "ref_alt": "T/C"},
    "rs7412": {"chrom": "19", "pos_grch37": 45412079, "ref_alt": "C/T"},
}

RS_RE = re.compile(r"(?<![A-Za-z0-9])(rs\d+)(?!\d)", re.I)
JHU_RE = re.compile(r"JHU_(?:chr)?([0-9XYMT]+)[.:](\d+)", re.I)
COORD_RE = re.compile(r"(?:^|[^0-9A-Za-z])(?:chr)?([0-9XYMT]+):(\d+)", re.I)


def read_headers(path: Path) -> list[str]:
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    raw = pd.read_csv(path, sep=sep, header=None, dtype=str)
    if raw.shape[1] == 1:
        vals = raw.iloc[:, 0].dropna().astype(str).str.strip().tolist()
        if vals and vals[0].lower() in {"header_name", "column", "column_name", "genomic_header"}:
            vals = vals[1:]
        return vals
    # Also support a one-row CSV containing the full header.
    return [str(v).strip() for v in raw.iloc[0].dropna().tolist()]


def parse_header(header: str):
    rsids = [m.group(1).lower() for m in RS_RE.finditer(header)]
    m = JHU_RE.search(header) or COORD_RE.search(header)
    if m:
        return rsids, m.group(1).upper().replace("CHR", ""), int(m.group(2))
    return rsids, None, None


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Find exact rs429358/rs7412 columns or their GRCh37 coordinate representations."
    )
    ap.add_argument("header_file", help="Full genomic header CSV/TSV")
    ap.add_argument("--output", default="APOE_candidates.csv")
    args = ap.parse_args()

    headers = read_headers(Path(args.header_file).resolve())
    rows = []

    for h in headers:
        rsids, chrom, pos = parse_header(h)
        for rsid, meta in TARGETS.items():
            match_type = None
            if rsid.lower() in rsids:
                match_type = "exact_rsid"
            elif chrom == meta["chrom"] and pos == meta["pos_grch37"]:
                match_type = "grch37_coordinate"
            if match_type:
                rows.append(
                    {
                        "target_rsid": rsid,
                        "header_name": h,
                        "match_type": match_type,
                        "grch37_coordinate": f"chr{meta['chrom']}:{meta['pos_grch37']}",
                        "expected_ref_alt": meta["ref_alt"],
                    }
                )

    out = pd.DataFrame(
        rows,
        columns=[
            "target_rsid",
            "header_name",
            "match_type",
            "grch37_coordinate",
            "expected_ref_alt",
        ],
    )
    out.to_csv(args.output, index=False)

    print(out.to_string(index=False) if len(out) else "NO APOE DEFINING COLUMNS FOUND")
    print(f"\nWrote {args.output}")
    for rsid in TARGETS:
        n = int((out["target_rsid"] == rsid).sum()) if len(out) else 0
        print(f"{rsid}: {n} candidate column(s)")

    if len(out) == 0 or any((out["target_rsid"] == r).sum() == 0 for r in TARGETS):
        raise SystemExit(
            "At least one APOE-defining SNP was not found. Do not substitute another APOE-region variant."
        )


if __name__ == "__main__":
    main()
