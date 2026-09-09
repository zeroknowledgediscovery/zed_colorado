#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description="Inspect candidate APOE rs429358/rs7412 genotype columns.")
    ap.add_argument("genomic_csv", help="Full genomic matrix CSV")
    ap.add_argument("candidates_csv", help="Output from 01_find_apoe_columns.py")
    ap.add_argument("--id-column", default="FID")
    ap.add_argument("--counts-output", default="APOE_column_value_counts.csv")
    ap.add_argument("--sample-output", default="APOE_first_rows.csv")
    ap.add_argument("--sample-n", type=int, default=20)
    args = ap.parse_args()

    cand = pd.read_csv(args.candidates_csv, dtype=str)
    if cand.empty:
        raise SystemExit("Candidate file is empty.")
    if "header_name" not in cand.columns:
        raise SystemExit("Candidate file lacks header_name.")

    cols = [args.id_column] + list(dict.fromkeys(cand["header_name"].dropna().tolist()))
    print("Loading only these columns:")
    for c in cols:
        print("  ", c)

    g = pd.read_csv(args.genomic_csv, usecols=cols)

    rows = []
    for c in cols[1:]:
        vc = g[c].value_counts(dropna=False).sort_index()
        print(f"\n=== {c} ===")
        print(vc.to_string())
        for value, n in vc.items():
            rows.append({"header_name": c, "value": value, "n": int(n)})

    pd.DataFrame(rows).to_csv(args.counts_output, index=False)
    g.head(args.sample_n).to_csv(args.sample_output, index=False)

    print(f"\nWrote {args.counts_output}")
    print(f"Wrote {args.sample_output}")
    print("\nSTOP HERE and verify which allele each 0/1/2 column counts before running step 03.")


if __name__ == "__main__":
    main()
