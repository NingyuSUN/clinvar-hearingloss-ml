"""Load the per-variant modelling matrix."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import MODELING_MATRIX
from .protocol import ALL_FEATURES

REQUIRED = ["VariationID", "y", "p4_gene_group", "p4_review_star", "GeneSymbol",
            "ReviewStatus", "consequence_class", "p4_headline_cohort",
            "p4_in_modeling_cohort", "p4_sens_ge2star", "p4_expert_panel", *ALL_FEATURES]


def load_matrix(path: str | Path | None = None) -> pd.DataFrame:
    """Read ``data/modeling_matrix.csv`` (or a ``.gz`` beside it) and check columns."""
    path = Path(path) if path else MODELING_MATRIX
    if not path.exists() and path.with_suffix(path.suffix + ".gz").exists():
        path = path.with_suffix(path.suffix + ".gz")
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Ship data/modeling_matrix.csv(.gz) or regenerate it "
            "with hlpath.protocol.build_matrix (see docs/data.md)."
        )
    df = pd.read_csv(path, low_memory=False)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"modelling matrix is missing columns: {missing}")
    df["y"] = pd.to_numeric(df["y"], errors="coerce").astype("Int64")
    return df
