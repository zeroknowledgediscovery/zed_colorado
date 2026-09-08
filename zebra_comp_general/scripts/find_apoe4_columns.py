#!/usr/bin/env python3
"""Locate columns containing the APOE epsilon-defining SNPs rs429358 and rs7412.

This utility intentionally does NOT infer APOE genotype or APOE4 carrier status.
Allele/genotype encoding varies across genomic exports, so the matched columns
must be inspected before constructing the focal APOE4_carrier variable.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

RSIDS = ("rs429358", "rs7412")


def read_header_names(path: Path) -> list[str]:
    # Supports the project's usual one-column genomicdataheader CSV and also a
    # single-row CSV containing many header names.
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    raw = pd.read_csv(path, sep=sep, header=None, dtype=str)
    if raw.shape[1] == 1:
        vals = raw.iloc[:, 0].dropna().astype(str).str.strip().tolist()
        if vals and vals[0].lower() in {"header_name", "column", "column_name", "genomic_header"}:
            vals = vals[1:]
        return vals
    return [str(v).strip() for v in raw.iloc[0].dropna().tolist()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Locate rs429358/rs7412 genomic header columns")
    ap.add_argument("header_file", help="One-column or single-row genomic header CSV/TSV")
    ap.add_argument("--output", default="APOE4_defining_columns.csv")
    args = ap.parse_args()

    headers = read_header_names(Path(args.header_file).resolve())
    rows = []
    for rsid in RSIDS:
        rgx = re.compile(re.escape(rsid), re.I)
        found = [h for h in headers if rgx.search(h)]
        if not found:
            rows.append({"rsid": rsid, "header_name": "", "status": "NOT_FOUND"})
        else:
            rows.extend({"rsid": rsid, "header_name": h, "status": "FOUND"} for h in found)

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)
    print(out.to_string(index=False))
    print(f"\nWrote {args.output}")
    print("Do not infer APOE4 carrier status until the allele/genotype encoding of these columns is verified.")


if __name__ == "__main__":
    main()
