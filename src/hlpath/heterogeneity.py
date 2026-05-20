"""Subgroup and boundary analyses."""
from __future__ import annotations

import pandas as pd

from .modeling import result_row, run_gene_split_xgb


def evaluate_subgroups(
    df: pd.DataFrame,
    feature_sets: dict[str, list[str]],
    subgroup_col: str,
    min_rows: int = 100,
    seed: int = 42,
) -> pd.DataFrame:
    """Run gene-split evaluation within each subgroup when enough rows exist."""
    rows = []
    for subgroup_value, sub_df in df.groupby(subgroup_col):
        if len(sub_df) < min_rows:
            continue
        for exp_name, features in feature_sets.items():
            try:
                res = run_gene_split_xgb(
                    sub_df,
                    features,
                    experiment=f"{exp_name} | {subgroup_col}={subgroup_value}",
                    seed=seed,
                )
                row = result_row(res)
                row[subgroup_col] = subgroup_value
                row["n_rows"] = len(sub_df)
                rows.append(row)
            except ValueError as exc:
                rows.append({
                    "Experiment": exp_name,
                    subgroup_col: subgroup_value,
                    "n_rows": len(sub_df),
                    "error": str(exc),
                })
    return pd.DataFrame(rows)


def conservation_bins(df: pd.DataFrame, col: str = "ensembl_conservation") -> pd.DataFrame:
    """Summarize label rate across conservation bins."""
    d = df.dropna(subset=[col, "label"]).copy()
    d["conservation_bin"] = pd.cut(
        d[col],
        bins=[-float("inf"), -2, -1, 0, 1, 2, float("inf")],
        labels=["<-2", "[-2,-1)", "[-1,0)", "[0,1)", "[1,2)", ">=2"],
    )
    return (
        d.groupby("conservation_bin", observed=True)
        .agg(n=("label", "size"), pathogenic_rate=("label", "mean"))
        .reset_index()
    )
