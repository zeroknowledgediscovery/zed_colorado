#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


def load_config(path: str | Path) -> dict[str, Any]:
    p = Path(path).resolve()
    cfg = json.loads(p.read_text())
    cfg["_config_path"] = str(p)
    cfg["_config_dir"] = str(p.parent)
    return cfg


def resolve_path(cfg: dict[str, Any], value: str | Path) -> Path:
    p = Path(value)
    if p.is_absolute():
        return p
    return (Path(cfg["_config_dir"]) / p).resolve()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    if suffix in {".csv", ".txt", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        return pd.read_csv(path, sep=sep)
    raise ValueError(f"Unsupported table type: {path}")


def normalize_id(df: pd.DataFrame, source_id: str, canonical: str = "patient_id") -> pd.DataFrame:
    if source_id not in df.columns:
        raise ValueError(f"Missing ID column {source_id!r}")
    out = df.copy()
    if source_id != canonical:
        if canonical in out.columns:
            raise ValueError(f"Both {source_id!r} and canonical ID {canonical!r} are present")
        out = out.rename(columns={source_id: canonical})
    return out


def check_unique_id(df: pd.DataFrame, name: str, id_col: str = "patient_id") -> None:
    if id_col not in df.columns:
        raise ValueError(f"{name} lacks {id_col}")
    ndup = int(df[id_col].duplicated().sum())
    if ndup:
        raise ValueError(f"{name} has {ndup} duplicated {id_col} values")


def collapse_score_rows(
    df: pd.DataFrame,
    id_col: str,
    score_col: str,
    policy: str = "identical",
) -> pd.DataFrame:
    x = df[[id_col, score_col]].copy()
    x[score_col] = pd.to_numeric(x[score_col], errors="coerce")
    if not x[id_col].duplicated().any():
        return x
    if policy == "error":
        raise ValueError("Score table has duplicated patient IDs")
    if policy == "mean":
        return x.groupby(id_col, as_index=False, sort=False)[score_col].mean()
    if policy == "first":
        return x.groupby(id_col, as_index=False, sort=False)[score_col].first()
    if policy != "identical":
        raise ValueError(f"Unknown score duplicate policy {policy!r}")

    conflicts: list[Any] = []
    for pid, s in x.groupby(id_col, sort=False)[score_col]:
        vals = s.dropna().to_numpy(float)
        if vals.size > 1 and not np.allclose(vals, vals[0], rtol=0.0, atol=1e-12):
            conflicts.append(pid)
            if len(conflicts) >= 10:
                break
    if conflicts:
        raise ValueError(
            "Score table contains non-identical duplicate values for patient IDs: "
            + ", ".join(map(str, conflicts))
        )

    def first_nonmissing(s: pd.Series) -> float:
        y = s.dropna()
        return float(y.iloc[0]) if len(y) else np.nan

    return x.groupby(id_col, as_index=False, sort=False)[score_col].agg(first_nonmissing)


def encode_binary_label(s: pd.Series, spec: dict[str, Any]) -> pd.Series:
    pos = set(str(v).strip().lower() for v in spec.get("positive_values", [1, "1", "Y", "y", True]))
    neg = set(str(v).strip().lower() for v in spec.get("negative_values", [0, "0", "N", "n", False]))
    missing = spec.get("missing", "drop")

    def enc(v: Any) -> float:
        if pd.isna(v):
            if missing == "negative":
                return 0.0
            if missing == "positive":
                return 1.0
            return np.nan
        t = str(v).strip().lower()
        if t in pos:
            return 1.0
        if t in neg:
            return 0.0
        try:
            q = float(v)
            if q == 1:
                return 1.0
            if q == 0:
                return 0.0
        except Exception:
            pass
        if spec.get("unknown", "error") == "drop":
            return np.nan
        raise ValueError(f"Cannot encode phenotype value {v!r} for {spec.get('name', spec.get('column'))}")

    return s.map(enc)


def apply_absent_row_policy(
    d: pd.DataFrame,
    phenos: list[dict[str, Any]],
    presence_column: str = "_phenotype_row_present",
) -> pd.DataFrame:
    """Apply phenotype semantics to genomic patients absent from the phenotype file.

    `missing` controls a null cell in an existing phenotype row. `absent_row`
    controls a patient who is not represented in the phenotype table at all.
    If `absent_row` is omitted it inherits `missing`, which preserves the common
    case-control convention that a sparse case/adjudication table plus
    `missing=negative` means patients absent from that table are controls.
    """
    out = d.copy()
    if presence_column not in out:
        raise ValueError(f"Internal phenotype-presence marker {presence_column!r} missing")
    absent = out[presence_column].isna()
    for spec in phenos:
        name = spec["name"]
        policy = spec.get("absent_row", spec.get("missing", "drop"))
        if policy == "negative":
            out.loc[absent, name] = 0.0
        elif policy == "positive":
            out.loc[absent, name] = 1.0
        elif policy == "drop":
            pass
        elif policy == "error":
            if absent.any():
                raise ValueError(
                    f"{int(absent.sum())} genomic/scored patients are absent from phenotype table "
                    f"for phenotype {name!r}; absent_row='error'"
                )
        else:
            raise ValueError(
                f"Unknown absent_row policy {policy!r} for phenotype {name!r}; "
                "use negative, positive, drop, or error"
            )
    return out


def _read_selector_file(cfg: dict[str, Any], selector: dict[str, Any]) -> list[str]:
    p = resolve_path(cfg, selector["path"])
    x = read_table(p)
    field = selector.get("column_field", "header_name")
    if field not in x.columns:
        if len(x.columns) == 1:
            field = x.columns[0]
        else:
            raise ValueError(f"Selector file {p} lacks field {field!r}")
    return [str(v) for v in x[field].dropna().tolist()]


def select_columns(columns: Iterable[str], cfg: dict[str, Any], spec: dict[str, Any]) -> list[str]:
    cols = list(columns)
    selectors = spec.get("include")
    if selectors is None:
        selectors = [spec]
    selected: list[str] = []
    for sel in selectors:
        typ = sel.get("type", "all")
        found: list[str]
        if typ == "all":
            found = cols.copy()
        elif typ == "columns":
            wanted = list(sel.get("columns", []))
            missing = [c for c in wanted if c not in cols]
            if missing and sel.get("allow_missing", False) is False:
                raise ValueError("Missing explicitly selected genomic columns: " + ", ".join(missing[:20]))
            found = [c for c in wanted if c in cols]
        elif typ == "file":
            wanted = _read_selector_file(cfg, sel)
            wanted_set = set(wanted)
            found = [c for c in cols if c in wanted_set]
        elif typ == "regex":
            rgx = re.compile(sel["pattern"])
            found = [c for c in cols if rgx.search(c)]
        elif typ == "prefix":
            pref = tuple(sel.get("prefixes", []))
            found = [c for c in cols if c.startswith(pref)]
        elif typ == "slice":
            start, end = sel["start_column"], sel["end_column"]
            if start not in cols or end not in cols:
                raise ValueError(f"Genomic slice endpoints {start!r}, {end!r} not both present")
            i, j = cols.index(start), cols.index(end)
            if i > j:
                i, j = j, i
            found = cols[i : j + 1]
        else:
            raise ValueError(f"Unknown genomic selector type {typ!r}")
        for c in found:
            if c not in selected:
                selected.append(c)

    excluded = set(spec.get("exclude_columns", []))
    if spec.get("exclude_regex"):
        rgx = re.compile(spec["exclude_regex"])
        excluded.update(c for c in selected if rgx.search(c))
    selected = [c for c in selected if c not in excluded]
    if not selected:
        raise ValueError("Genomic selector resolved to zero columns")
    return selected


def decode_gene_signal(df: pd.DataFrame, spec: dict[str, Any]) -> pd.Series:
    if not spec or not spec.get("enabled", True):
        return pd.Series(np.nan, index=df.index, name="gene_signal")
    typ = spec.get("type", "binary_column")
    if typ == "binary_column":
        c = spec["column"]
        if c not in df:
            raise ValueError(f"Missing gene signal column {c!r}")
        x = df[c]
        if "positive_values" in spec:
            pos = set(str(v).strip().lower() for v in spec["positive_values"])
            neg = set(str(v).strip().lower() for v in spec.get("negative_values", []))
            out = pd.Series(np.nan, index=x.index, dtype=float)
            txt = x.astype(str).str.strip().str.lower()
            out.loc[txt.isin(pos)] = 1.0
            if neg:
                out.loc[txt.isin(neg)] = 0.0
            else:
                out.loc[x.notna() & ~txt.isin(pos)] = 0.0
            return out.rename("gene_signal")
        return pd.to_numeric(x, errors="coerce").rename("gene_signal")
    if typ == "dosage_column":
        c = spec["column"]
        x = pd.to_numeric(df[c], errors="coerce")
        return (x >= float(spec.get("positive_threshold", 1))).where(x.notna()).astype(float).rename("gene_signal")
    if typ == "genotype_column":
        c = spec["column"]
        pos = set(str(v).upper() for v in spec["positive_values"])
        called = set(str(v).upper() for v in spec.get("called_values", []))
        x = df[c]
        txt = x.astype(str).str.upper()
        out = pd.Series(np.nan, index=x.index, dtype=float)
        out.loc[x.notna() & txt.isin(pos)] = 1.0
        if called:
            out.loc[x.notna() & txt.isin(called - pos)] = 0.0
        else:
            out.loc[x.notna() & ~txt.isin(pos)] = 0.0
        return out.rename("gene_signal")
    if typ == "onehot_any":
        called_cols = list(spec["called_columns"])
        positive_cols = list(spec["positive_columns"])
        miss = [c for c in called_cols if c not in df.columns]
        if miss:
            raise ValueError("Missing gene one-hot columns: " + ", ".join(miss))
        x = df[called_cols].apply(pd.to_numeric, errors="coerce")
        if spec.get("require_exactly_one", True):
            sm = x.sum(axis=1, min_count=1)
            if (sm > 1).any():
                raise ValueError("Gene-state one-hot encoding contains multi-hot rows")
            called = sm == 1
        else:
            called = x.notna().any(axis=1)
        pos = x[positive_cols].fillna(0).sum(axis=1) > 0
        out = pd.Series(np.nan, index=df.index, dtype=float)
        out.loc[called] = pos.loc[called].astype(float)
        return out.rename("gene_signal")
    raise ValueError(f"Unknown gene signal type {typ!r}")


def encode_genomic_matrix(raw: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    for c in columns:
        s = raw[c]
        num = pd.to_numeric(s, errors="coerce")
        nonmissing = int(s.notna().sum())
        convertible = int(num.notna().sum())
        if pd.api.types.is_numeric_dtype(s) or nonmissing == 0 or convertible / max(nonmissing, 1) >= 0.95:
            pieces.append(pd.DataFrame({c: num.astype(float)}, index=raw.index))
        else:
            d = pd.get_dummies(s.astype("string"), prefix=c, prefix_sep="__", dummy_na=False, dtype=float)
            pieces.append(d)
    out = pd.concat(pieces, axis=1)
    nunique = out.nunique(dropna=True)
    out = out.loc[:, nunique > 1]
    if out.shape[1] == 0:
        raise ValueError("All selected genomic features are invariant")
    return out


def build_analysis_frame(cfg: dict[str, Any]) -> tuple[pd.DataFrame, list[str], dict[str, Any]]:
    inputs = cfg["inputs"]
    id_name = cfg.get("canonical_id", "patient_id")

    gspec = inputs["genomics"]
    pspec = inputs["phenotypes"]
    zspec = inputs["zebra"]
    graw = normalize_id(read_table(resolve_path(cfg, gspec["path"])), gspec.get("id_column", id_name), id_name)
    praw = normalize_id(read_table(resolve_path(cfg, pspec["path"])), pspec.get("id_column", id_name), id_name)
    zraw = normalize_id(read_table(resolve_path(cfg, zspec["path"])), zspec.get("id_column", id_name), id_name)

    check_unique_id(graw, "genomics", id_name)
    check_unique_id(praw, "phenotypes", id_name)
    score_col = zspec["score_column"]
    if score_col not in zraw:
        raise ValueError(f"Score table lacks {score_col!r}")
    zraw = collapse_score_rows(zraw, id_name, score_col, zspec.get("duplicate_policy", "identical"))
    if score_col != "zebra_score":
        zraw = zraw.rename(columns={score_col: "zebra_score"})

    primary_name = cfg.get("primary_phenotype")
    phenos = cfg.get("phenotypes", [])
    if not phenos:
        raise ValueError("Config needs at least one phenotype specification")
    if primary_name is None:
        primary_name = phenos[0]["name"]
    pcols = [p["column"] for p in phenos]
    missing_p = [c for c in pcols if c not in praw]
    if missing_p:
        raise ValueError("Missing phenotype columns: " + ", ".join(missing_p))

    gene = decode_gene_signal(graw, cfg.get("gene_signal", {}))
    graw = graw.copy()
    graw["gene_signal"] = gene

    selector = cfg.get("genomic_features", {"type": "all"})
    raw_feature_cols = select_columns(
        [c for c in graw.columns if c not in {id_name, "gene_signal"}], cfg, selector
    )
    forbidden = {p["column"] for p in phenos} | {"target", "zebra_score", "predicted_risk"}
    raw_feature_cols = [c for c in raw_feature_cols if c not in forbidden]
    Xg = encode_genomic_matrix(graw, raw_feature_cols)
    Xg.insert(0, id_name, graw[id_name].values)
    Xg["gene_signal"] = gene.values

    p = praw[[id_name] + pcols].copy()
    p["_phenotype_row_present"] = 1
    for spec in phenos:
        p[spec["name"]] = encode_binary_label(p[spec["column"]], spec)
    p = p[[id_name, "_phenotype_row_present"] + [x["name"] for x in phenos]]

    d = Xg.merge(p, on=id_name, how="left").merge(zraw[[id_name, "zebra_score"]], on=id_name, how="left")
    d["zebra_score"] = pd.to_numeric(d["zebra_score"], errors="coerce")
    if cfg.get("drop_missing_score", True):
        d = d.loc[d.zebra_score.notna()].copy()
    if cfg.get("score_range"):
        lo, hi = cfg["score_range"]
        if not d.zebra_score.dropna().between(float(lo), float(hi)).all():
            raise ValueError(f"ZeBRA score outside configured range [{lo}, {hi}]")

    # Apply phenotype-row absence semantics *after* the left merge. This is the
    # critical distinction for sparse case/adjudication tables: an absent row
    # may mean control rather than unknown. The frozen ILD analysis uses this
    # convention explicitly via fillna(0) after merging the phenotype table.
    d = apply_absent_row_policy(d, phenos)
    phenotype_row_present_n = int(d["_phenotype_row_present"].notna().sum())
    phenotype_row_absent_n = int(d["_phenotype_row_present"].isna().sum())
    d = d.drop(columns="_phenotype_row_present")

    if primary_name not in d:
        raise ValueError(f"Primary phenotype {primary_name!r} not present")
    d = d.loc[d[primary_name].notna()].copy()
    d["target"] = d[primary_name].astype(int)
    d["gene_signal"] = pd.to_numeric(d["gene_signal"], errors="coerce")

    for filt in cfg.get("filters", []):
        c, op, v = filt["column"], filt.get("op", "eq"), filt.get("value")
        if c not in d:
            raise ValueError(f"Filter column {c!r} not available")
        if op == "eq": d = d.loc[d[c] == v]
        elif op == "ne": d = d.loc[d[c] != v]
        elif op == "notnull": d = d.loc[d[c].notna()]
        elif op == "in": d = d.loc[d[c].isin(v)]
        else: raise ValueError(f"Unsupported filter op {op!r}")

    feature_cols = [c for c in Xg.columns if c not in {id_name, "gene_signal"}]
    meta = {
        "analysis_name": cfg.get("analysis_name", "zebra_genomics"),
        "primary_phenotype": primary_name,
        "raw_genomic_columns": raw_feature_cols,
        "encoded_genomic_columns": feature_cols,
        "n": int(len(d)),
        "cases": int(d.target.sum()),
        "controls": int((d.target == 0).sum()),
        "phenotype_row_present_n_after_score_filter": phenotype_row_present_n,
        "phenotype_row_absent_n_after_score_filter": phenotype_row_absent_n,
        "gene_called_n": int(d.gene_signal.notna().sum()),
        "gene_positive_n": int(d.gene_signal.fillna(0).sum()) if d.gene_signal.notna().any() else 0,
    }
    return d.reset_index(drop=True), feature_cols, meta


def metric_summary(values: pd.Series) -> dict[str, float]:
    x = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    if len(x) == 0:
        return {"mean": np.nan, "sd": np.nan, "median": np.nan, "q025": np.nan, "q975": np.nan}
    return {
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "median": float(np.median(x)),
        "q025": float(np.quantile(x, .025)),
        "q975": float(np.quantile(x, .975)),
    }


def finite_or_nan(v: Any) -> float:
    try:
        x = float(v)
        return x if math.isfinite(x) else np.nan
    except Exception:
        return np.nan
