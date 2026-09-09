#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate APOE4 carrier counts, optionally stratified by ADRD phenotype.")
    ap.add_argument("apoe_csv")
    ap.add_argument("--id-column", default="FID")
    ap.add_argument("--phenotype-file")
    ap.add_argument("--phenotype-id-column", default="FID")
    ap.add_argument("--phenotype-column")
    ap.add_argument("--positive-values", default="1,Y,y")
    args = ap.parse_args()

    x = pd.read_csv(args.apoe_csv)
    if args.id_column not in x.columns or "APOE4_carrier" not in x.columns:
        raise SystemExit(f"Expected columns {args.id_column!r} and 'APOE4_carrier'.")

    called = x["APOE4_carrier"].notna()
    carriers = x["APOE4_carrier"] == 1

    print(f"Total patients: {len(x):,}")
    print(f"APOE4-called patients: {int(called.sum()):,}")
    print(f"APOE4 carriers: {int(carriers.sum()):,}")
    if called.any():
        print(f"APOE4 carrier prevalence: {x.loc[called, 'APOE4_carrier'].mean():.4f}")

    if not args.phenotype_file:
        return
    if not args.phenotype_column:
        raise SystemExit("--phenotype-column is required when --phenotype-file is provided.")

    p = pd.read_csv(args.phenotype_file)
    needed = [args.phenotype_id_column, args.phenotype_column]
    missing = [c for c in needed if c not in p.columns]
    if missing:
        raise SystemExit("Phenotype file missing columns: " + ", ".join(missing))

    p = p[needed].copy()
    if args.phenotype_id_column != args.id_column:
        p = p.rename(columns={args.phenotype_id_column: args.id_column})

    d = x.merge(p, on=args.id_column, how="inner")
    pos = {v.strip().lower() for v in args.positive_values.split(",")}
    ytxt = d[args.phenotype_column].astype(str).str.strip().str.lower()
    d["_case"] = ytxt.isin(pos)
    dcalled = d[d["APOE4_carrier"].notna()].copy()

    print(f"\nJoined phenotype patients: {len(d):,}")
    print(f"ADRD cases among APOE-called patients: {int(dcalled['_case'].sum()):,}")

    for label, m in [("cases", dcalled["_case"]), ("controls", ~dcalled["_case"])]:
        z = dcalled.loc[m, "APOE4_carrier"]
        print(f"{label}: n={len(z):,}, APOE4 carriers={int((z == 1).sum()):,}, prevalence={z.mean():.4f}" if len(z) else f"{label}: n=0")


if __name__ == "__main__":
    main()
