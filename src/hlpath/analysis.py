"""Turn raw per-fold output into the reported tables: the additive feature ladder
with paired tests, the pooled out-of-fold confusion with a bootstrap CI, the
consequence-stratified breakdown, and the gene-split / random-split gap.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from .metrics import bootstrap_ci, paired_delta

LADDER = ["base_conservation_only", "base_consequence_only", "base_freq_only",
          "L1_allele", "L2_base", "L3_consequence", "L4_constraint", "L5_frequency",
          "L6_full", "F_minus_consequence", "F_minus_allele", "F_minus_conservation",
          "F_minus_constraint", "F_minus_frequency", "F_minus_domain"]
LADDER_STEPS = [("L2_base", "L3_consequence"), ("L3_consequence", "L4_constraint"),
                ("L4_constraint", "L5_frequency"), ("L5_frequency", "L6_full"),
                ("L6_full", "F_minus_consequence"), ("L6_full", "F_minus_frequency"),
                ("L6_full", "F_minus_constraint"), ("L6_full", "F_minus_conservation"),
                ("L6_full", "F_minus_domain"), ("L6_full", "F_minus_allele")]


def _pooled(per_fold: pd.DataFrame, oof: pd.DataFrame) -> pd.DataFrame:
    """Attach each row's own frozen fold threshold and majority-vote it across seeds."""
    thr = per_fold[["seed", "outer_fold", "threshold"]]
    merged = oof.merge(thr, on=["seed", "outer_fold"], how="left")
    merged["hit"] = (merged["prob"] >= merged["threshold"]).astype(float)
    return merged.groupby("VariationID").agg(
        y=("y", "first"), prob=("prob", "mean"), hit=("hit", "mean"),
        consequence_class=("consequence_class", "first"),
    ).reset_index()


def ladder_table(per_fold: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name in LADDER:
        g = per_fold[per_fold["feature_set"] == name]
        if not len(g):
            continue
        rows.append({
            "feature_set": name,
            "n_features": int(g["n_features"].iloc[0]),
            "n_folds": len(g),
            "auc_per_fold_mean": round(g["auc"].mean(), 4),
            "auc_per_fold_sd": round(g["auc"].std(ddof=1), 4),
            "best_iteration_median": int(g["best_iteration"].median()),
            "R90_recall_mean": round(g.get("R90_recall", pd.Series(dtype=float)).mean(), 4),
            "R90_precision_mean": round(g.get("R90_precision", pd.Series(dtype=float)).mean(), 4),
        })
    return pd.DataFrame(rows)


def paired_ladder_table(per_fold: pd.DataFrame) -> pd.DataFrame:
    rows = [paired_delta(per_fold, a, b) for a, b in LADDER_STEPS]
    return pd.DataFrame([r for r in rows if r])


def pooled_oof_table(per_fold: pd.DataFrame, oof_by_set: dict, feature_set: str = "L6_full",
                     n_boot: int = 2000) -> dict:
    pooled = _pooled(per_fold[per_fold["feature_set"] == feature_set], oof_by_set[feature_set])
    pred = (pooled["hit"] >= 0.5).to_numpy().astype(int)
    y = pooled["y"].to_numpy().astype(int)
    ci = bootstrap_ci(y, pooled["prob"].to_numpy(), pred, n_boot=n_boot)
    return {
        "feature_set": feature_set,
        "n_variants": int(len(pooled)),
        "bootstrap_ci": ci,
        "confusion": {
            "TP": int(((pred == 1) & (y == 1)).sum()),
            "FP": int(((pred == 1) & (y == 0)).sum()),
            "TN": int(((pred == 0) & (y == 0)).sum()),
            "FN": int(((pred == 0) & (y == 1)).sum()),
        },
    }


def stratified_table(per_fold: pd.DataFrame, oof_by_set: dict,
                     feature_set: str = "L6_full") -> pd.DataFrame:
    pooled = _pooled(per_fold[per_fold["feature_set"] == feature_set], oof_by_set[feature_set])
    pooled["pred"] = (pooled["hit"] >= 0.5).astype(int)
    rows = []
    for cc, g in pooled.groupby("consequence_class"):
        auc = roc_auc_score(g["y"], g["prob"]) if g["y"].nunique() > 1 else np.nan
        positives = g[g["y"] == 1]
        called = g[g["pred"] == 1]
        rows.append({
            "consequence_class": cc,
            "n": len(g),
            "prevalence": round(g["y"].mean(), 3),
            "auc": round(auc, 4) if pd.notna(auc) else None,
            "R90_recall": round(positives["pred"].mean(), 3) if len(positives) else None,
            "R90_precision": round(called["y"].mean(), 3) if len(called) else None,
        })
    return pd.DataFrame(rows)


def split_gap(gene_per_fold, gene_oof, random_per_fold, random_oof,
              feature_set: str = "L6_full") -> dict:
    def pooled_auc(pf, oof):
        p = _pooled(pf[pf["feature_set"] == feature_set], oof[feature_set])
        return float(roc_auc_score(p["y"], p["prob"]))

    gene = pooled_auc(gene_per_fold, gene_oof)
    rand = pooled_auc(random_per_fold, random_oof)
    return {"feature_set": feature_set, "gene_held_out_auc": round(gene, 4),
            "random_split_auc": round(rand, 4), "gap": round(rand - gene, 4)}
