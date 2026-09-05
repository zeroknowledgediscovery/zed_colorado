#!/usr/bin/env python3

"""
Map phenotype-associated genes to column/header names in
genomicdataheader.csv.

The genomicdataheader.csv file contains ONE variant/header name per row.

Examples:
    rs2076295_T
    rs35705950.1_G
    JHU_2.179634520_C
    2:179634520-CT_C

Outputs:
    <PANEL>_genomic_header_matches.csv

with exactly:

    gene_name,header_name

Uses Ensembl GRCh37:
    gene -> genomic interval
    gene interval -> known rsIDs

It additionally directly maps coordinate-form header names to genes.
"""

import argparse
import csv
import json
import re
import time
from pathlib import Path
from collections import defaultdict

import pandas as pd
import requests


# ================================================================
# PHENOTYPE PANELS
# ================================================================

PANELS = {

    # ------------------------------------------------------------
    # HFrEF / dilated cardiomyopathy
    # ------------------------------------------------------------
    "HFREF": [
        "TTN",
        "BAG3",
        "LMNA",
        "RBM20",
        "FLNC",
        "DSP",
        "PLN",
        "MYH7",
        "TNNT2",
        "TNNC1",
        "SCN5A",
        "DES",
        "ACTN2",
        "HSPB7",
        "NKX2-5",
    ],

    # ------------------------------------------------------------
    # IPF / fibrotic ILD
    #
    # Combines:
    #   common IPF susceptibility loci
    #   telomere biology genes
    #   familial pulmonary fibrosis genes
    #   surfactant genes
    # ------------------------------------------------------------
    "IPF": [
        # strongest common IPF locus
        "MUC5B",

        # telomere biology
        "TERT",
        "TERC",
        "RTEL1",
        "PARN",
        "STN1",       # formerly OBFC1
        "DKC1",
        "TINF2",
        "NAF1",
        "ZCCHC8",

        # surfactant / familial pulmonary fibrosis
        "SFTPA1",
        "SFTPA2",
        "SFTPC",
        "ABCA3",

        # replicated common IPF susceptibility loci
        "DSP",
        "FAM13A",
        "DPP9",
        "TOLLIP",
        "ATP11A",

        # additional repeatedly implicated loci
        "AKAP13",
        "DEPTOR",
        "KIF15",
        "MAD1L1",
        "IVD",
        "ZKSCAN1",
    ],

    # ------------------------------------------------------------
    # Alzheimer's disease + related dementias
    #
    # Includes:
    #   Mendelian AD
    #   late-onset AD risk loci
    #   microglial/lipid AD genes
    #   FTD / tauopathy
    #   Lewy-body dementia
    #   selected hereditary vascular dementia genes
    # ------------------------------------------------------------
    "ADRD": [
        # core Alzheimer's
        "APOE",
        "APP",
        "PSEN1",
        "PSEN2",

        # strong AD rare/common risk genes
        "TREM2",
        "SORL1",
        "ABCA7",
        "BIN1",
        "CR1",
        "CLU",
        "PICALM",
        "CD33",

        # additional strong GWAS / functional loci
        "INPP5D",
        "PLCG2",
        "ABI3",
        "CD2AP",
        "EPHA1",
        "FERMT2",
        "PTK2B",
        "ADAM10",
        "ACE",
        "SPI1",
        "MS4A6A",
        "ABCA1",

        # FTD / tau-related dementia
        "MAPT",
        "GRN",
        "C9ORF72",
        "TMEM106B",
        "TBK1",

        # Lewy-body / synuclein-related dementia
        "SNCA",
        "GBA1",

        # hereditary vascular dementia
        "NOTCH3",
        "HTRA1",
    ],
}


ENSEMBL = "https://grch37.rest.ensembl.org"

HTTP_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}


# ================================================================
# Ensembl HTTP
# ================================================================

def ensembl_get(endpoint, params=None, retries=6):

    url = ENSEMBL + endpoint

    for attempt in range(retries):

        r = requests.get(
            url,
            params=params,
            headers=HTTP_HEADERS,
            timeout=120,
        )

        if r.status_code == 200:
            return r.json()

        if r.status_code == 429:
            wait = float(r.headers.get("Retry-After", 2))
            time.sleep(wait)
            continue

        if r.status_code >= 500:
            time.sleep(2 ** attempt)
            continue

        raise RuntimeError(
            f"\nEnsembl returned {r.status_code}"
            f"\nURL: {r.url}"
            f"\n{r.text[:1000]}"
        )

    raise RuntimeError(
        f"Ensembl failed after {retries} attempts: {url}"
    )


# ================================================================
# Gene coordinates
# ================================================================

def get_gene_region(gene):

    data = ensembl_get(
        f"/lookup/symbol/homo_sapiens/{gene}"
    )

    return {
        "gene": gene,
        "ensembl_id": data["id"],
        "chromosome":
            str(data["seq_region_name"]).replace("chr", ""),
        "start": int(data["start"]),
        "end": int(data["end"]),
    }


# ================================================================
# Known variants overlapping gene
# ================================================================

def get_gene_variants(region):

    chrom = region["chromosome"]
    start = region["start"]
    end = region["end"]

    # Ensembl overlap API has a maximum interval size,
    # so split very large genes if necessary.
    chunk_size = 4_000_000

    ids = set()

    chunk_start = start

    while chunk_start <= end:

        chunk_end = min(
            chunk_start + chunk_size - 1,
            end
        )

        interval = (
            f"{chrom}:{chunk_start}-{chunk_end}"
        )

        data = ensembl_get(
            f"/overlap/region/homo_sapiens/{interval}",
            params={"feature": "variation"},
        )

        for variant in data:

            vid = variant.get("id")

            if vid:
                ids.add(str(vid).lower())

        chunk_start = chunk_end + 1

    return ids


# ================================================================
# Parse genomicdataheader names
# ================================================================

RS_RE = re.compile(
    r"(rs\d+)",
    re.IGNORECASE
)

# examples:
# JHU_2.179634520_C
# JHU_6.7563232_G

JHU_RE = re.compile(
    r"JHU_(?:chr)?([0-9XYMT]+)[.:](\d+)",
    re.IGNORECASE
)

# examples:
# 2:179634520-CT_C
# 6:7568027-CT_C
# chr2:179634520-A_G

COORD_RE = re.compile(
    r"(?:chr)?([0-9XYMT]+):(\d+)",
    re.IGNORECASE
)


def parse_header(header):

    header = header.strip()

    rsid = None
    chrom = None
    pos = None

    m = RS_RE.search(header)

    if m:
        rsid = m.group(1).lower()

    m = JHU_RE.search(header)

    if m:

        chrom = m.group(1).upper()
        pos = int(m.group(2))

        return rsid, chrom, pos

    m = COORD_RE.search(header)

    if m:

        chrom = m.group(1).upper()
        pos = int(m.group(2))

    return rsid, chrom, pos


# ================================================================
# Construct gene / variant lookup
# ================================================================

def build_gene_reference(
    genes,
    cache_file
):

    if cache_file.exists():

        print(
            f"Loading cache: {cache_file}"
        )

        with open(cache_file) as f:
            cache = json.load(f)

    else:
        cache = {}

    regions = {}

    rsid_to_genes = defaultdict(set)

    changed = False

    unavailable_genes = []

    for gene in genes:

        print(f"\n[{gene}]")

        try:

            if gene in cache:

                print("  cached")

                region = cache[gene]["region"]

                variants = set(
                    cache[gene]["variants"]
                )

            else:

                print(
                    "  finding GRCh37 interval..."
                )

                region = get_gene_region(gene)

                print(
                    f"  chr{region['chromosome']}:"
                    f"{region['start']:,}-"
                    f"{region['end']:,}"
                )

                print(
                    "  retrieving overlapping variants..."
                )

                variants = get_gene_variants(
                    region
                )

                print(
                    f"  {len(variants):,} variants"
                )

                cache[gene] = {
                    "region": region,
                    "variants": sorted(variants),
                }

                changed = True

                time.sleep(0.15)

            regions[gene] = region

            for variant in variants:

                if variant.startswith("rs"):

                    rsid_to_genes[
                        variant
                    ].add(gene)

        except Exception as e:

            print(
                f"  WARNING: could not resolve "
                f"{gene}: {e}"
            )

            unavailable_genes.append(gene)

    if changed:

        with open(cache_file, "w") as f:

            json.dump(
                cache,
                f
            )

        print(
            f"\nSaved cache: {cache_file}"
        )

    if unavailable_genes:

        print(
            "\nUnresolved genes:"
        )

        print(
            ", ".join(unavailable_genes)
        )

    return regions, rsid_to_genes


# ================================================================
# Coordinate matching
# ================================================================

def genes_at_coordinate(
    chrom,
    pos,
    regions,
    padding=0
):

    if chrom is None or pos is None:
        return []

    chrom = (
        chrom.replace("chr", "").upper()
    )

    hits = []

    for gene, region in regions.items():

        rchrom = (
            str(region["chromosome"])
            .replace("chr", "")
            .upper()
        )

        if chrom != rchrom:
            continue

        if (
            region["start"] - padding
            <= pos
            <= region["end"] + padding
        ):

            hits.append(gene)

    return hits


# ================================================================
# Scan one-column genomicdataheader.csv
# ================================================================

def scan_headers(
    header_file,
    regions,
    rsid_to_genes,
    padding=0
):

    results = set()

    n = 0

    with open(
        header_file,
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as f:

        reader = csv.reader(f)

        for row in reader:

            if not row:
                continue

            header_name = row[0].strip()

            if not header_name:
                continue

            n += 1

            if n % 100000 == 0:

                print(
                    f"Scanned {n:,} headers; "
                    f"{len(results):,} matches"
                )

            rsid, chrom, pos = (
                parse_header(header_name)
            )

            matched_genes = set()

            # rsID-based matching
            if rsid is not None:

                matched_genes.update(
                    rsid_to_genes.get(
                        rsid,
                        set()
                    )
                )

            # coordinate-based matching
            if (
                chrom is not None
                and pos is not None
            ):

                matched_genes.update(
                    genes_at_coordinate(
                        chrom,
                        pos,
                        regions,
                        padding,
                    )
                )

            for gene in matched_genes:

                results.add(
                    (
                        gene,
                        header_name,
                    )
                )

    print(
        f"\nScanned {n:,} headers"
    )

    return results


# ================================================================
# Run one panel
# ================================================================

def run_panel(
    panel_name,
    header_file,
    padding
):

    genes = PANELS[panel_name]

    print("\n")
    print("=" * 70)
    print(panel_name)
    print("=" * 70)

    print(
        f"\n{len(genes)} genes:"
    )

    print(
        ", ".join(genes)
    )

    cache_file = Path(
        f"{panel_name}_"
        f"ensembl_grch37_cache.json"
    )

    regions, rsid_to_genes = (
        build_gene_reference(
            genes,
            cache_file
        )
    )

    print(
        f"\nKnown rsIDs indexed: "
        f"{len(rsid_to_genes):,}"
    )

    matches = scan_headers(
        header_file,
        regions,
        rsid_to_genes,
        padding,
    )

    df = pd.DataFrame(
        sorted(matches),
        columns=[
            "gene_name",
            "header_name",
        ],
    )

    output_file = (
        f"{panel_name}_"
        f"genomic_header_matches.csv"
    )

    df.to_csv(
        output_file,
        index=False
    )

    print("\nMATCHES BY GENE")
    print("-" * 50)

    counts = (
        df.groupby("gene_name")
        .size()
        .reindex(
            genes,
            fill_value=0
        )
        .rename("n_headers")
    )

    print(
        counts.to_string()
    )

    print(
        f"\nTOTAL MATCHES: "
        f"{len(df):,}"
    )

    print(
        f"Saved: {output_file}"
    )

    return df


# ================================================================
# Main
# ================================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "header_file",
        help=(
            "genomicdataheader.csv "
            "(one header per row)"
        ),
    )

    parser.add_argument(
        "--panel",
        choices=[
            "HFREF",
            "IPF",
            "ADRD",
            "ALL",
        ],
        default="IPF",
    )

    parser.add_argument(
        "--padding",
        type=int,
        default=0,
        help=(
            "Include this many bp upstream/"
            "downstream of each gene. "
            "Default = 0."
        ),
    )

    args = parser.parse_args()

    if args.panel == "ALL":

        panels = [
            "HFREF",
            "IPF",
            "ADRD",
        ]

    else:

        panels = [
            args.panel
        ]

    for panel in panels:

        run_panel(
            panel,
            args.header_file,
            args.padding,
        )


if __name__ == "__main__":
    main()
