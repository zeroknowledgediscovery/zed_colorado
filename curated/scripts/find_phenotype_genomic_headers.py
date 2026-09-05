#!/usr/bin/env python3
"""
Map phenotype-associated genes to variant/header names in a one-column
``genomicdataheader.csv`` file.

The Colorado genomic header contains one variant identifier per row, for example::

    rs2076295_T
    rs35705950.1_G
    JHU_2.179634520_C
    2:179634520-CT_C

The script uses Ensembl GRCh37 to map each curated gene to its genomic interval
and to retrieve known rsIDs overlapping that interval. It then scans the full
header and writes exactly two columns::

    gene_name,header_name

Named panels are provided for HFrEF/dilated cardiomyopathy, IPF/fibrotic ILD,
and Alzheimer's disease and related dementias (ADRD).
"""

import argparse
import csv
import json
import re
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests


PANELS = {
    "HFREF": [
        "TTN", "BAG3", "LMNA", "RBM20", "FLNC", "DSP", "PLN", "MYH7",
        "TNNT2", "TNNC1", "SCN5A", "DES", "ACTN2", "HSPB7", "NKX2-5",
    ],
    "IPF": [
        "MUC5B",
        "TERT", "TERC", "RTEL1", "PARN", "STN1", "DKC1", "TINF2", "NAF1",
        "ZCCHC8",
        "SFTPA1", "SFTPA2", "SFTPC", "ABCA3",
        "DSP", "FAM13A", "DPP9", "TOLLIP", "ATP11A",
        "AKAP13", "DEPTOR", "KIF15", "MAD1L1", "IVD", "ZKSCAN1",
    ],
    "ADRD": [
        "APOE", "APP", "PSEN1", "PSEN2",
        "TREM2", "SORL1", "ABCA7", "BIN1", "CR1", "CLU", "PICALM", "CD33",
        "INPP5D", "PLCG2", "ABI3", "CD2AP", "EPHA1", "FERMT2", "PTK2B",
        "ADAM10", "ACE", "SPI1", "MS4A6A", "ABCA1",
        "MAPT", "GRN", "C9ORF72", "TMEM106B", "TBK1",
        "SNCA", "GBA1",
        "NOTCH3", "HTRA1",
    ],
}

ENSEMBL = "https://grch37.rest.ensembl.org"
HTTP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

RS_RE = re.compile(r"(rs\d+)", re.IGNORECASE)
JHU_RE = re.compile(r"JHU_(?:chr)?([0-9XYMT]+)[.:](\d+)", re.IGNORECASE)
COORD_RE = re.compile(r"(?:chr)?([0-9XYMT]+):(\d+)", re.IGNORECASE)


def ensembl_get(endpoint, params=None, retries=6):
    """GET JSON from the Ensembl GRCh37 REST API with basic retry handling."""
    url = ENSEMBL + endpoint
    for attempt in range(retries):
        response = requests.get(url, params=params, headers=HTTP_HEADERS, timeout=120)
        if response.status_code == 200:
            return response.json()
        if response.status_code == 429:
            time.sleep(float(response.headers.get("Retry-After", 2)))
            continue
        if response.status_code >= 500:
            time.sleep(2 ** attempt)
            continue
        raise RuntimeError(
            f"Ensembl returned {response.status_code}\nURL: {response.url}\n"
            f"{response.text[:1000]}"
        )
    raise RuntimeError(f"Ensembl failed after {retries} attempts: {url}")


def get_gene_region(gene):
    data = ensembl_get(f"/lookup/symbol/homo_sapiens/{gene}")
    return {
        "gene": gene,
        "ensembl_id": data["id"],
        "chromosome": str(data["seq_region_name"]).replace("chr", ""),
        "start": int(data["start"]),
        "end": int(data["end"]),
    }


def get_gene_variants(region):
    """Return known variant IDs overlapping a GRCh37 gene interval."""
    chrom = region["chromosome"]
    start = region["start"]
    end = region["end"]

    chunk_size = 4_000_000
    ids = set()
    chunk_start = start
    while chunk_start <= end:
        chunk_end = min(chunk_start + chunk_size - 1, end)
        interval = f"{chrom}:{chunk_start}-{chunk_end}"
        data = ensembl_get(
            f"/overlap/region/homo_sapiens/{interval}",
            params={"feature": "variation"},
        )
        for variant in data:
            variant_id = variant.get("id")
            if variant_id:
                ids.add(str(variant_id).lower())
        chunk_start = chunk_end + 1
    return ids


def parse_header(header):
    """Extract an rsID and/or chromosome/position from a Colorado header."""
    header = header.strip()
    rsid = None
    chrom = None
    pos = None

    match = RS_RE.search(header)
    if match:
        rsid = match.group(1).lower()

    match = JHU_RE.search(header)
    if match:
        return rsid, match.group(1).upper(), int(match.group(2))

    match = COORD_RE.search(header)
    if match:
        chrom = match.group(1).upper()
        pos = int(match.group(2))

    return rsid, chrom, pos


def build_gene_reference(genes, cache_file):
    """Build gene intervals and an rsID -> gene lookup, caching Ensembl calls."""
    if cache_file.exists():
        print(f"Loading cache: {cache_file}")
        with cache_file.open() as handle:
            cache = json.load(handle)
    else:
        cache = {}

    regions = {}
    rsid_to_genes = defaultdict(set)
    changed = False
    unresolved = []

    for gene in genes:
        print(f"\n[{gene}]")
        try:
            if gene in cache:
                print("  cached")
                region = cache[gene]["region"]
                variants = set(cache[gene]["variants"])
            else:
                print("  finding GRCh37 interval...")
                region = get_gene_region(gene)
                print(
                    f"  chr{region['chromosome']}:"
                    f"{region['start']:,}-{region['end']:,}"
                )
                print("  retrieving overlapping variants...")
                variants = get_gene_variants(region)
                print(f"  {len(variants):,} variants")
                cache[gene] = {"region": region, "variants": sorted(variants)}
                changed = True
                time.sleep(0.15)

            regions[gene] = region
            for variant in variants:
                if variant.startswith("rs"):
                    rsid_to_genes[variant].add(gene)
        except Exception as exc:
            print(f"  WARNING: could not resolve {gene}: {exc}")
            unresolved.append(gene)

    if changed:
        with cache_file.open("w") as handle:
            json.dump(cache, handle)
        print(f"\nSaved cache: {cache_file}")

    if unresolved:
        print("\nUnresolved genes: " + ", ".join(unresolved))

    return regions, rsid_to_genes


def genes_at_coordinate(chrom, pos, regions, padding=0):
    if chrom is None or pos is None:
        return []
    chrom = chrom.replace("chr", "").upper()
    hits = []
    for gene, region in regions.items():
        region_chrom = str(region["chromosome"]).replace("chr", "").upper()
        if chrom != region_chrom:
            continue
        if region["start"] - padding <= pos <= region["end"] + padding:
            hits.append(gene)
    return hits


def scan_headers(header_file, regions, rsid_to_genes, padding=0):
    """Scan every row of the one-column genomic header file."""
    results = set()
    n = 0

    with open(header_file, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            header_name = row[0].strip()
            if not header_name:
                continue

            n += 1
            if n % 100_000 == 0:
                print(f"Scanned {n:,} headers; {len(results):,} matches")

            rsid, chrom, pos = parse_header(header_name)
            matched_genes = set()

            if rsid is not None:
                matched_genes.update(rsid_to_genes.get(rsid, set()))
            if chrom is not None and pos is not None:
                matched_genes.update(
                    genes_at_coordinate(chrom, pos, regions, padding=padding)
                )

            for gene in matched_genes:
                results.add((gene, header_name))

    print(f"\nScanned {n:,} headers")
    return results


def run_panel(panel_name, header_file, padding, output_dir):
    genes = PANELS[panel_name]
    print("\n" + "=" * 70)
    print(panel_name)
    print("=" * 70)
    print(f"\n{len(genes)} genes:")
    print(", ".join(genes))

    output_dir.mkdir(parents=True, exist_ok=True)
    cache_file = output_dir / f"{panel_name}_ensembl_grch37_cache.json"
    regions, rsid_to_genes = build_gene_reference(genes, cache_file)
    print(f"\nKnown rsIDs indexed: {len(rsid_to_genes):,}")

    matches = scan_headers(header_file, regions, rsid_to_genes, padding=padding)
    df = pd.DataFrame(sorted(matches), columns=["gene_name", "header_name"])
    output_file = output_dir / f"{panel_name}_genomic_header_matches.csv"
    df.to_csv(output_file, index=False)

    print("\nMATCHES BY GENE")
    print("-" * 50)
    if df.empty:
        counts = pd.Series(0, index=genes, name="n_headers")
    else:
        counts = (
            df.groupby("gene_name")
            .size()
            .reindex(genes, fill_value=0)
            .rename("n_headers")
        )
    print(counts.to_string())
    print(f"\nTOTAL MATCHES: {len(df):,}")
    print(f"Saved: {output_file}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Map curated phenotype genes to Colorado genomic header names."
    )
    parser.add_argument(
        "header_file",
        help="genomicdataheader.csv containing one header/variant name per row",
    )
    parser.add_argument(
        "--panel",
        choices=["HFREF", "IPF", "ADRD", "ALL"],
        default="IPF",
        help="Curated gene panel to map (default: IPF)",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=0,
        help="bp upstream/downstream of each gene to include (default: 0)",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for match CSVs and Ensembl caches (default: current directory)",
    )
    args = parser.parse_args()

    panels = ["HFREF", "IPF", "ADRD"] if args.panel == "ALL" else [args.panel]
    for panel in panels:
        run_panel(
            panel,
            args.header_file,
            padding=args.padding,
            output_dir=Path(args.output_dir),
        )


if __name__ == "__main__":
    main()
