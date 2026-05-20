#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hlpath.data import BASIC_FEATURES, prepare_model_table
from hlpath.modeling import compare_feature_sets


def main():
    parser = argparse.ArgumentParser(description="Evaluate baseline vs conservation under gene-based split.")
    parser.add_argument("--data", required=True, help="Processed CSV with label, GeneSymbol, and ensembl_conservation.")
    parser.add_argument("--out", default="results/tables/conservation_gene_split_results.csv")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    df = prepare_model_table(args.data)
    feature_sets = {
        "Baseline (no conservation)": BASIC_FEATURES,
        "+Conservation": BASIC_FEATURES + ["ensembl_conservation"],
    }
    results = compare_feature_sets(df, feature_sets, seed=args.seed)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print(results.to_string(index=False))
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
