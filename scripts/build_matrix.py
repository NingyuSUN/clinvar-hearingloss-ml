#!/usr/bin/env python3
"""Rebuild data/modeling_matrix.csv(.gz) from a raw annotated feature table.

    python scripts/build_matrix.py --raw path/to/annotated_features.csv

The raw table needs one row per ClinVar variant with: VariationID, GeneSymbol,
ClinicalSignificance, ReviewStatus, p1_status, p2_af_status, p2_max_alt_af_pass,
p3_variant_is_coding, and the six feature groups (see docs/data.md). The upstream
VEP / gnomAD annotation steps are external to this repository.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import pandas as pd

from hlpath.config import MODELING_MATRIX
from hlpath.protocol import build_matrix


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=MODELING_MATRIX)
    ap.add_argument("--gzip", action="store_true", help="also write a .gz copy")
    args = ap.parse_args()

    raw = pd.read_csv(args.raw, low_memory=False)
    matrix = build_matrix(raw)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.out, index=False)
    print(f"wrote {args.out}  ({len(matrix)} rows, {len(matrix.columns)} columns)")
    for flag in ["p4_headline_cohort", "p4_in_modeling_cohort", "p4_sens_ge2star"]:
        print(f"  {flag}: {int((matrix[flag] == 1).sum())}")
    if args.gzip:
        with open(args.out, "rb") as f, gzip.open(str(args.out) + ".gz", "wb") as g:
            shutil.copyfileobj(f, g)
        print(f"wrote {args.out}.gz")


if __name__ == "__main__":
    main()
