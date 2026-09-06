#!/usr/bin/env python3
"""Map phenotype-associated genes to rows in a one-column genomic header CSV.

The header file is expected to contain one variant/header name per row, e.g.
rs2076295_T, JHU_2.179634520_C, or 2:179634520-CT_C. Gene intervals and
overlapping known variants are obtained from the Ensembl GRCh37 REST API.
Output contains exactly: gene_name,header_name.
"""

from __future__ import annotations
import argparse, csv, json, re, time
from pathlib import Path
from collections import defaultdict
import pandas as pd
import requests

PANELS = {
    "HFREF": ["TTN","BAG3","LMNA","RBM20","FLNC","DSP","PLN","MYH7","TNNT2","TNNC1","SCN5A","DES","ACTN2","HSPB7","NKX2-5"],
    "IPF": ["MUC5B","TERT","TERC","RTEL1","PARN","STN1","DKC1","TINF2","NAF1","ZCCHC8","SFTPA1","SFTPA2","SFTPC","ABCA3","DSP","FAM13A","DPP9","TOLLIP","ATP11A","AKAP13","DEPTOR","KIF15","MAD1L1","IVD","ZKSCAN1"],
    "ADRD": ["APOE","APP","PSEN1","PSEN2","TREM2","SORL1","ABCA7","BIN1","CR1","CLU","PICALM","CD33","INPP5D","PLCG2","ABI3","CD2AP","EPHA1","FERMT2","PTK2B","ADAM10","ACE","SPI1","MS4A6A","ABCA1","MAPT","GRN","C9ORF72","TMEM106B","TBK1","SNCA","GBA1","NOTCH3","HTRA1"],
}

ENSEMBL = "https://grch37.rest.ensembl.org"
HTTP_HEADERS = {"Content-Type":"application/json", "Accept":"application/json"}
RS_RE = re.compile(r"(rs\d+)", re.I)
JHU_RE = re.compile(r"JHU_(?:chr)?([0-9XYMT]+)[.:](\d+)", re.I)
COORD_RE = re.compile(r"(?:chr)?([0-9XYMT]+):(\d+)", re.I)

def ensembl_get(endpoint, params=None, retries=6):
    url = ENSEMBL + endpoint
    for attempt in range(retries):
        r = requests.get(url, params=params, headers=HTTP_HEADERS, timeout=120)
        if r.status_code == 200:
            return r.json()
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 2))); continue
        if r.status_code >= 500:
            time.sleep(2 ** attempt); continue
        raise RuntimeError(f"Ensembl {r.status_code}: {r.url}\n{r.text[:500]}")
    raise RuntimeError(f"Ensembl failed after {retries} attempts: {url}")

def get_gene_region(gene):
    d = ensembl_get(f"/lookup/symbol/homo_sapiens/{gene}")
    return {"gene":gene, "ensembl_id":d["id"], "chromosome":str(d["seq_region_name"]).replace("chr",""), "start":int(d["start"]), "end":int(d["end"])}

def get_gene_variants(region):
    ids, chunk_size = set(), 4_000_000
    s, end, chrom = region["start"], region["end"], region["chromosome"]
    while s <= end:
        e = min(s + chunk_size - 1, end)
        data = ensembl_get(f"/overlap/region/homo_sapiens/{chrom}:{s}-{e}", params={"feature":"variation"})
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

def build_reference(genes, cache_file):
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    regions, rsid_to_genes, changed = {}, defaultdict(set), False
    for gene in genes:
        if gene in cache:
            region = cache[gene]["region"]
            variants = set(cache[gene]["variants"])
        else:
            print(f"[{gene}] resolving GRCh37 interval and variants")
            region = get_gene_region(gene)
            variants = get_gene_variants(region)
            cache[gene] = {"region":region, "variants":sorted(variants)}
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
    chrom = chrom.replace("chr","").upper()
    return [g for g,r in regions.items() if str(r["chromosome"]).replace("chr","").upper() == chrom and r["start"]-padding <= pos <= r["end"]+padding]

def run_panel(panel, header_file, padding):
    genes = PANELS[panel]
    regions, rsmap = build_reference(genes, Path(f"{panel}_ensembl_grch37_cache.json"))
    matches = set()
    with open(header_file, newline="", encoding="utf-8-sig") as f:
        for row in csv.reader(f):
            if not row or not row[0].strip():
                continue
            name = row[0].strip()
            rsid, chrom, pos = parse_header(name)
            gs = set(rsmap.get(rsid, set())) if rsid else set()
            gs.update(genes_at_coordinate(chrom, pos, regions, padding))
            matches.update((g, name) for g in gs)
    out = pd.DataFrame(sorted(matches), columns=["gene_name","header_name"])
    outfile = f"{panel}_genomic_header_matches.csv"
    out.to_csv(outfile, index=False)
    print(f"{panel}: {len(out):,} gene/header pairs -> {outfile}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("header_file")
    ap.add_argument("--panel", choices=["HFREF","IPF","ADRD","ALL"], default="IPF")
    ap.add_argument("--padding", type=int, default=0)
    a = ap.parse_args()
    panels = list(PANELS) if a.panel == "ALL" else [a.panel]
    for panel in panels:
        run_panel(panel, a.header_file, a.padding)

if __name__ == "__main__":
    main()
