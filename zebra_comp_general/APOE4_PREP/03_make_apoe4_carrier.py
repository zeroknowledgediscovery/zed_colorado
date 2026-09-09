#!/usr/bin/env python3
from __future__ import annotations

import argparse
import numpy as np
import pandas as pd


def to_c_dosage(x: pd.Series, dosage_allele: str) -> pd.Series:
    y = pd.to_numeric(x, errors="coerce")
    bad = y.dropna()[~y.dropna().isin([0, 1, 2])]
    if len(bad):
        raise ValueError(f"Observed non-0/1/2 genotype dosage values: {sorted(bad.unique())[:20]}")
    allele = dosage_allele.upper()
    if allele == "C":
        return y
    if allele == "T":
        return 2 - y
    raise ValueError("dosage allele must be C or T")


def main() -> None:
    ap = argparse.ArgumentParser(description="Construct conservative APOE4 carrier status from rs429358 and rs7412 dosages.")
    ap.add_argument("genomic_csv")
    ap.add_argument("--id-column", default="FID")
    ap.add_argument("--rs429358-column", required=True)
    ap.add_argument("--rs429358-dosage-allele", required=True, choices=["C", "T", "c", "t"])
    ap.add_argument("--rs7412-column", required=True)
    ap.add_argument("--rs7412-dosage-allele", required=True, choices=["C", "T", "c", "t"])
    ap.add_argument("--double-het-e2e4", action="store_true", help="Treat C-dosage pair (1,1) as epsilon-2/epsilon-4 carrier instead of ambiguous.")
    ap.add_argument("--output", default="APOE4_carrier.csv")
    args = ap.parse_args()

    usecols = [args.id_column, args.rs429358_column, args.rs7412_column]
    g = pd.read_csv(args.genomic_csv, usecols=usecols)

    d429 = to_c_dosage(g[args.rs429358_column], args.rs429358_dosage_allele)
    d741 = to_c_dosage(g[args.rs7412_column], args.rs7412_dosage_allele)

    out = pd.DataFrame({
        args.id_column: g[args.id_column],
        "rs429358_C_dosage": d429,
        "rs7412_C_dosage": d741,
    })
    out["APOE_genotype"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out["APOE4_carrier"] = np.nan

    called = d429.notna() & d741.notna()

    # Canonical two-SNP APOE states using C dosages at rs429358 and rs7412.
    states = {
        (0, 0): ("e2/e2", 0.0),
        (0, 1): ("e2/e3", 0.0),
        (0, 2): ("e3/e3", 0.0),
        (1, 2): ("e3/e4", 1.0),
        (2, 2): ("e4/e4", 1.0),
    }
    if args.double_het_e2e4:
        states[(1, 1)] = ("e2/e4", 1.0)

    for (a, b), (gt, carrier) in states.items():
        m = called & (d429 == a) & (d741 == b)
        out.loc[m, "APOE_genotype"] = gt
        out.loc[m, "APOE4_carrier"] = carrier

    # Double heterozygotes remain explicitly ambiguous unless user opts in.
    if not args.double_het_e2e4:
        m = called & (d429 == 1) & (d741 == 1)
        out.loc[m, "APOE_genotype"] = "ambiguous_double_het"

    # Any other called dosage pair is noncanonical and remains missing for carrier status.
    m_noncanonical = called & out["APOE_genotype"].isna()
    out.loc[m_noncanonical, "APOE_genotype"] = "noncanonical"

    out.to_csv(args.output, index=False)

    print(f"N total: {len(out):,}")
    print(f"Two-SNP called: {int(called.sum()):,}")
    print(f"APOE4 carrier called: {int(out['APOE4_carrier'].notna().sum()):,}")
    print(f"APOE4 carriers: {int((out['APOE4_carrier'] == 1).sum()):,}")
    if out["APOE4_carrier"].notna().any():
        print(f"APOE4 prevalence among called: {out['APOE4_carrier'].mean():.4f}")
    print("\nGenotype counts:")
    print(out["APOE_genotype"].value_counts(dropna=False).to_string())
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
