#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hlpath.data import BASIC_FEATURES, prepare_model_table
from hlpath.modeling import result_row, run_gene_split_xgb, run_random_split_xgb


def main():
    parser = argparse.ArgumentParser(description="Compare random row split with strict gene-based split.")
    parser.add_argument("--data", required=True, help="CSV with label, GeneSymbol, and variant-level features.")
    parser.add_argument("--out", default="results/tables/random_vs_gene_split_results.csv")
    parser.add_argument("--use-conservation", action="store_true", help="Include ensembl_conservation if present.")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    df = prepare_model_table(args.data)
    features = list(BASIC_FEATURES)
    if args.use_conservation and "ensembl_conservation" in df.columns:
        features.append("ensembl_conservation")

    random_res = run_random_split_xgb(df, features, experiment="Random split", seed=args.seed)
    gene_res = run_gene_split_xgb(df, features, experiment="Gene-based split", seed=args.seed)

    results = pd.DataFrame([result_row(random_res), result_row(gene_res)])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print(results.to_string(index=False))
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
