#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hlpath.data import BASIC_FEATURES, prepare_model_table
from hlpath.modeling import run_gene_split_xgb
from hlpath.shap_analysis import gray_zone, save_dependence_plot, save_summary_plot, tree_shap_values


def main():
    parser = argparse.ArgumentParser(description="Train +Conservation model and save SHAP plots.")
    parser.add_argument("--data", required=True)
    parser.add_argument("--figdir", default="results/figures")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    df = prepare_model_table(args.data)
    features = BASIC_FEATURES + ["ensembl_conservation"]
    res = run_gene_split_xgb(df, features, experiment="+Conservation", seed=args.seed)

    X_test = res.test_pred[features].copy()
    # Use median imputation for SHAP plotting consistency.
    X_test = X_test.fillna(X_test.median(numeric_only=True))

    explainer, shap_values = tree_shap_values(res.booster, X_test)
    save_summary_plot(shap_values, X_test, figdir / "shap_global_summary.png")
    save_dependence_plot(shap_values, X_test, "ensembl_conservation", figdir / "shap_conservation_dependence.png")

    X_gray = gray_zone(X_test)
    if len(X_gray) > 0:
        _, shap_gray = tree_shap_values(res.booster, X_gray)
        save_summary_plot(shap_gray, X_gray, figdir / "shap_gray_zone_summary.png")

    print(f"AUC={res.auc:.3f}, R90={res.r90}")
    print(f"Saved SHAP figures to {figdir}")


if __name__ == "__main__":
    main()
