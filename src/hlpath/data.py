"""Data loading and feature engineering for hearing loss variant pathogenicity models."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

BASIC_FEATURES = [
    "variant_length",
    "has_del",
    "has_ins",
    "has_dup",
    "has_fs",
    "has_stop",
    "has_splice",
]

CONSERVATION_FEATURES = ["ensembl_conservation"]

GENE_CONSTRAINT_FEATURES = [
    "gene_constraint_oe_lof_upper",
    "gene_constraint_oe_mis_upper",
]


def load_table(path: str | Path) -> pd.DataFrame:
    """Load a CSV/TSV/parquet table based on file extension."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="	")
    return pd.read_csv(path)


def choose_hgvs_column(df: pd.DataFrame) -> Optional[str]:
    """Return the most likely HGVS/name column."""
    for col in ["hgvs", "HGVS", "Hgvs", "Name"]:
        if col in df.columns:
            return col
    return None


def make_basic_flags(text: object) -> tuple[int, int, int, int, int]:
    """Parse deletion/insertion/duplication/frameshift/stop flags from HGVS-like text."""
    if not isinstance(text, str):
        text = ""
    lower = text.lower()
    return (
        int("del" in lower),
        int("ins" in lower),
        int("dup" in lower),
        int(("fs" in lower) or ("frameshift" in lower)),
        int(("*" in text) or ("ter" in lower) or ("stop" in lower)),
    )


def has_splice_from_hgvs(text: object) -> int:
    """Detect canonical splice-site patterns from HGVS-like text.

    This follows the stricter definition used in the later notebooks:
    IVS, +1, +2, -1, or -2.
    """
    if not isinstance(text, str):
        return 0
    upper = text.upper()
    return int(
        ("IVS" in upper)
        or ("+1" in upper)
        or ("+2" in upper)
        or ("-1" in upper)
        or ("-2" in upper)
    )


def add_variant_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add variant-level engineered features if missing."""
    out = df.copy()

    if "label" in out.columns:
        out = out.dropna(subset=["label"]).copy()
        out["label"] = pd.to_numeric(out["label"], errors="coerce").astype(int)

    if "variant_length" not in out.columns:
        if {"Start", "Stop"}.issubset(out.columns):
            start = pd.to_numeric(out["Start"], errors="coerce")
            stop = pd.to_numeric(out["Stop"], errors="coerce")
            out["variant_length"] = (stop - start).abs().fillna(1).clip(lower=1).astype(int)
        else:
            out["variant_length"] = 1

    hgvs_col = choose_hgvs_column(out)
    if hgvs_col is not None:
        flags = out[hgvs_col].apply(make_basic_flags)
        out[["has_del", "has_ins", "has_dup", "has_fs", "has_stop"]] = pd.DataFrame(
            flags.tolist(), index=out.index
        )
        out["has_splice"] = out[hgvs_col].apply(has_splice_from_hgvs)
    else:
        for col in ["has_del", "has_ins", "has_dup", "has_fs", "has_stop", "has_splice"]:
            if col not in out.columns:
                out[col] = 0

    out["is_splice_variant"] = (out["has_splice"] == 1).astype(int)
    out["is_truncating"] = ((out["has_fs"] == 1) | (out["has_stop"] == 1)).astype(int)
    out["variant_class"] = np.where(out["is_truncating"] == 1, "truncating", "non_truncating")

    for col in BASIC_FEATURES:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    return out


def feature_columns(
    use_conservation: bool = True,
    use_gene_constraint: bool = False,
    df: Optional[pd.DataFrame] = None,
) -> list[str]:
    """Return feature list, optionally filtering by columns present in df."""
    cols = list(BASIC_FEATURES)
    if use_conservation:
        cols += CONSERVATION_FEATURES
    if use_gene_constraint:
        cols += GENE_CONSTRAINT_FEATURES

    if df is not None:
        cols = [c for c in cols if c in df.columns]
    return cols


def prepare_model_table(path: str | Path, require_gene: bool = True) -> pd.DataFrame:
    """Load table, add features, and run minimal sanity checks."""
    df = load_table(path)
    df = add_variant_features(df)
    if require_gene and "GeneSymbol" not in df.columns:
        raise ValueError("Missing required column: GeneSymbol")
    if "label" not in df.columns:
        raise ValueError("Missing required column: label")
    return df.dropna(subset=["label", "GeneSymbol"]).copy()
