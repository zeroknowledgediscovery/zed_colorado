#!/usr/bin/env python3
"""Map disease-associated genes to raw genomic header names using Ensembl GRCh37.

Input is a one-column CSV/TSV containing raw genomic header names. Headers may
contain rsIDs or GRCh37 coordinates (e.g. JHU_2.179634520_C or 2:179634520-CT_C).
The gene panels are stored in ../gene_panels.json so they can be edited without
changing code.
"""
from __future__ import annotations

import argparse
import json
import re
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
PANEL_FILE = ROOT / "gene_panels.json"
ENSEMBL = "https://grch37.rest.ensembl.org"
HTTP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
RS_RE = re.compile(r"(rs\d+)", re.I)
JHU_RE = re.compile(r"JHU_(?:chr)?([0-9XYMT]+)[.:](\d+)", re.I)
COORD_RE = re.compile(r"(?:chr)?([0-9XYMT]+):(\d+)", re.I)


def load_panels() -> dict[str, list[str]]:
    return json.loads(PANEL_FILE.read_text())


def ensembl_get(endpoint, params=None, retries=6):
    url = ENSEMBL + endpoint
    for attempt in range(retries):
        r = requests.get(url, params=params, headers=HTTP_HEADERS, timeout=120)
        if r.status_code == 200:
            return r.json()
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 2)))
            continue
        if r.status_code >= 500:
            time.sleep(2 ** attempt)
            continue
        raise RuntimeError(f"Ensembl {r.status_code}: {r.url}\n{r.text[:500]}")
    raise RuntimeError(f"Ensembl failed after {retries} attempts: {url}")


def get_gene_region(gene):
    d = ensembl_get(f"/lookup/symbol/homo_sapiens/{gene}")
    return {
        "gene": gene,
        "ensembl_id": d["id"],
        "chromosome": str(d["seq_region_name"]).replace("chr", ""),
        "start": int(d["start"]),
        "end": int(d["end"]),
    }


def get_gene_variants(region):
    ids, chunk_size = set(), 4_000_000
    s, end, chrom = region["start"], region["end"], region["chromosome"]
    while s <= end:
        e = min(s + chunk_size - 1, end)
        data = ensembl_get(
            f"/overlap/region/homo_sapiens/{chrom}:{s}-{e}",
            params={"feature": "variation"},
        )
        ids.update(str(v["id"]).lower() for v in data if v.get("id"))
        s = e + 1
    return ids


def parse_header(header):
    header = header.strip()
    mrs = RS_RE.search(header)
    rsid = mrs.group(1).lower() if mrs else None
    m = JHU_RE.search(header) or COORD_RE.search(header)
    if m:
        return rsid, m.group(1).upper(), int(m.group(2))
    return rsid, None, None


def build_reference(panel: str, genes: list[str], cache_dir: Path):
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{panel}_ensembl_grch37_cache.json"
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    regions, rsid_to_genes, changed = {}, defaultdict(set), False
    for gene in genes:
        if gene in cache:
            region = cache[gene]["region"]
            variants = set(cache[gene]["variants"])
        else:
            print(f"[{gene}] resolving GRCh37 interval and variants", flush=True)
            region = get_gene_region(gene)
            variants = get_gene_variants(region)
            cache[gene] = {"region": region, "variants": sorted(variants)}
            changed = True
            time.sleep(.15)
        regions[gene] = region
        for v in variants:
            if v.startswith("rs"):
                rsid_to_genes[v].add(gene)
    if changed:
        cache_file.write_text(json.dumps(cache))
    return regions, rsid_to_genes


def genes_at_coordinate(chrom, pos, regions, padding):
    if chrom is None or pos is None:
        return []
    chrom = chrom.replace("chr", "").upper()
    return [
        g
        for g, r in regions.items()
        if str(r["chromosome"]).replace("chr", "").upper() == chrom
        and r["start"] - padding <= pos <= r["end"] + padding
    ]


def read_headers(path: Path) -> list[str]:
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    raw = pd.read_csv(path, sep=sep, header=None)
    vals = [str(v).strip() for v in raw.iloc[:, 0].dropna().tolist()]
    if vals and vals[0].lower() in {"header_name", "column", "column_name", "genomic_header"}:
        vals = vals[1:]
    return vals


def run_panel(panel, header_file, padding, output, cache_dir):
    panels = load_panels()
    genes = panels[panel]
    regions, rsmap = build_reference(panel, genes, cache_dir)
    matches = set()
    for name in read_headers(header_file):
        rsid, chrom, pos = parse_header(name)
        gs = set(rsmap.get(rsid, set())) if rsid else set()
        gs.update(genes_at_coordinate(chrom, pos, regions, padding))
        matches.update((g, name) for g in gs)
    out = pd.DataFrame(sorted(matches), columns=["gene_name", "header_name"])
    out.to_csv(output, index=False)
    print(f"{panel}: {len(out):,} gene/header pairs -> {output}")


def main():
    panels = load_panels()
    ap = argparse.ArgumentParser()
    ap.add_argument("header_file")
    ap.add_argument("--panel", choices=sorted(list(panels) + ["ALL"]), default="IPF")
    ap.add_argument("--padding", type=int, default=0)
    ap.add_argument("--output", help="Output CSV; with --panel ALL, interpreted as output directory")
    ap.add_argument("--cache-dir", default=str(ROOT / ".ensembl_cache"))
    a = ap.parse_args()
    header = Path(a.header_file).resolve()
    cache = Path(a.cache_dir).resolve()
    chosen = list(panels) if a.panel == "ALL" else [a.panel]
    for panel in chosen:
        if a.panel == "ALL":
            od = Path(a.output).resolve() if a.output else Path.cwd()
            od.mkdir(parents=True, exist_ok=True)
            output = od / f"{panel}_genomic_header_matches.csv"
        else:
            output = Path(a.output).resolve() if a.output else Path.cwd() / f"{panel}_genomic_header_matches.csv"
        run_panel(panel, header, a.padding, output, cache)


if __name__ == "__main__":
    main()
