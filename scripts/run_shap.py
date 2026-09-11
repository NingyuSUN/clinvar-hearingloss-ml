#!/usr/bin/env python3
"""TreeSHAP feature attribution on the headline cohort and the missense subset.

    python scripts/run_shap.py --out results/shap/

Requires shap and matplotlib (not core dependencies -- see docs/shap.md).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap

from hlpath.config import RESULTS_DIR
from hlpath.features import load_matrix
from hlpath.shap_analysis import (FEATS, build_shap_long, by_consequence, conservation_bands,
                                  group_importance, importance_table, per_variant)


def beeswarm(agg: pd.DataFrame, path: Path, title: str) -> None:
    sv = agg[[f"shap__{f}" for f in FEATS]].to_numpy()
    xv = agg[[f"val__{f}" for f in FEATS]].to_numpy()
    expl = shap.Explanation(values=sv, data=xv, feature_names=FEATS,
                            base_values=agg["shap_base"].to_numpy())
    plt.figure(figsize=(7.5, 6))
    shap.plots.beeswarm(expl, max_display=18, show=False)
    plt.title(title, fontsize=11)
    plt.tight_layout()
    plt.savefig(path, dpi=140)
    plt.close()


def dependence(agg: pd.DataFrame, feat: str, path: Path, title: str) -> None:
    x, s = agg[f"val__{feat}"].to_numpy(), agg[f"shap__{feat}"].to_numpy()
    plt.figure(figsize=(5.4, 4))
    plt.scatter(x, s, s=9, alpha=0.5, color="#0c7d6d")
    plt.axhline(0, color="#999", lw=0.8)
    plt.xlabel(feat)
    plt.ylabel(f"SHAP value for {feat}")
    plt.title(title, fontsize=10)
    plt.tight_layout()
    plt.savefig(path, dpi=140)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / "shap")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    fig_dir = args.out / "figures"
    fig_dir.mkdir(exist_ok=True)

    matrix = load_matrix()
    cohorts = {"headline": ("p4_headline_cohort", None),
              "missense": ("p4_headline_cohort", ["coding_nontruncating"])}
    aggs = {}
    for tag, (flag, consequence) in cohorts.items():
        long = build_shap_long(matrix, flag, consequence)
        agg = per_variant(long)
        aggs[tag] = agg
        importance_table(agg).to_csv(args.out / f"importance_{tag}.csv", index=False)
        group_importance(agg).to_csv(args.out / f"group_importance_{tag}.csv", index=False)
        beeswarm(agg, fig_dir / f"beeswarm_{tag}.png", f"TreeSHAP -- {tag}")
        print(f"{tag}: {len(agg)} variants")

    by_consequence(aggs["headline"]).to_csv(args.out / "by_consequence.csv", index=False)
    conservation_bands(aggs["missense"]).to_csv(args.out / "conservation_bands.csv", index=False)
    for feat in ["ensembl_conservation", "p4_freq_log10"]:
        dependence(aggs["missense"], feat, fig_dir / f"dependence_missense_{feat}.png",
                  f"missense: {feat}")
    print(f"wrote tables + figures to {args.out}")


if __name__ == "__main__":
    main()
