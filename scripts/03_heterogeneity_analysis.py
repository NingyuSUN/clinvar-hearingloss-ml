#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from hlpath.data import BASIC_FEATURES, prepare_model_table
from hlpath.heterogeneity import conservation_bins, evaluate_subgroups


def main():
    parser = argparse.ArgumentParser(description="Run subgroup and conservation-bin analyses.")
    parser.add_argument("--data", required=True)
    parser.add_argument("--outdir", default="results/tables")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = prepare_model_table(args.data)
    feature_sets = {
        "Baseline": BASIC_FEATURES,
        "+Conservation": BASIC_FEATURES + ["ensembl_conservation"],
    }

    bins = conservation_bins(df)
    bins.to_csv(outdir / "conservation_bins.csv", index=False)

    for col in ["variant_class", "is_splice_variant"]:
        res = evaluate_subgroups(df, feature_sets, subgroup_col=col, seed=args.seed)
        res.to_csv(outdir / f"subgroup_{col}.csv", index=False)

    print(f"Saved subgroup analyses to {outdir}")


if __name__ == "__main__":
    main()
