"""TreeSHAP feature attribution on the frozen full model (`L6_full`).

Every fold is retrained exactly as `pipeline.run_experiment` trains it (same
balanced split, same post-hoc log-loss iteration selection), so the explained
model is identical to the one the reported AUCs come from. SHAP values are taken
with XGBoost's native `pred_contribs=True` (exact TreeSHAP; identical to
`shap.TreeExplainer`) on each fold's **out-of-fold outer-test rows** — an honest,
not in-sample, attribution.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import xgboost as xgb

from .config import INNER_VALID_FRAC, K_OUTER, NUM_BOOST_ROUND, POSTHOC_MIN_ITER, SEEDS, xgb_params
from .pipeline import _impute, balanced_group_folds
from .protocol import FEATURE_SETS

FEATS = FEATURE_SETS["L6_full"]
GROUPS = {
    "p1_has_frameshift": "consequence", "p1_has_stop_gained": "consequence",
    "p1_has_canonical_splice": "consequence", "p1_has_splice_region": "consequence",
    "variant_length": "allele", "has_del": "allele", "has_ins": "allele", "has_dup": "allele",
    "ensembl_conservation": "conservation",
    "p3_gene_constraint_oe_lof_upper": "constraint", "p3_gene_constraint_oe_mis_upper": "constraint",
    "p4_freq_log10": "frequency", "p2_af_known_missing": "frequency",
    "p2_not_observed_in_gnomad_r4": "frequency", "p2_low_an_flag": "frequency",
    "p3_in_specific_domain": "domain", "p3_specific_domain_count": "domain",
    "p3_variant_is_coding": "domain",
}


def _fold_shap(d, outer, inner, out_fold, seed):
    te = d[d["_grp"].map(outer) == out_fold]
    tr_all = d[d["_grp"].map(outer) != out_fold]
    iv = tr_all[tr_all["_grp"].map(inner) == 1]
    it = tr_all[tr_all["_grp"].map(inner) == 0]
    if not len(iv) or not len(it):
        return None
    xit, xiv, xte = _impute(it[FEATS], iv[FEATS], te[FEATS])
    yit = it["y"].astype(int).to_numpy()
    yiv = iv["y"].astype(int).to_numpy()
    hist: dict = {}
    bst = xgb.train(xgb_params(seed), xgb.DMatrix(xit, label=yit, feature_names=FEATS),
                    num_boost_round=NUM_BOOST_ROUND,
                    evals=[(xgb.DMatrix(xiv, label=yiv, feature_names=FEATS), "iv")],
                    evals_result=hist, verbose_eval=False)
    logloss = np.asarray(hist["iv"]["logloss"], dtype=float)
    best_iter = (int(np.argmin(logloss[POSTHOC_MIN_ITER:]) + POSTHOC_MIN_ITER)
                 if len(logloss) > POSTHOC_MIN_ITER else int(np.argmin(logloss)))
    dte = xgb.DMatrix(xte, feature_names=FEATS)
    prob = bst.predict(dte, iteration_range=(0, best_iter + 1))
    contribs = bst.predict(dte, pred_contribs=True, iteration_range=(0, best_iter + 1))

    out = pd.DataFrame({"VariationID": te["VariationID"].to_numpy(), "seed": seed,
                        "outer_fold": out_fold, "y": te["y"].astype(int).to_numpy(),
                        "consequence_class": te["consequence_class"].to_numpy(),
                        "review_star": te["p4_review_star"].to_numpy(),
                        "prob": prob, "shap_base": contribs[:, -1]})
    for i, f in enumerate(FEATS):
        out[f"val__{f}"] = xte[f].to_numpy()
        out[f"shap__{f}"] = contribs[:, i]
    return out


def build_shap_long(matrix: pd.DataFrame, cohort_flag: str,
                    consequence: list[str] | None = None) -> pd.DataFrame:
    """SHAP values for every (variant, seed): one row per outer-test appearance."""
    d = matrix[matrix[cohort_flag] == 1].reset_index(drop=True).copy()
    if consequence is not None:
        d = d[d["consequence_class"].isin(consequence)].reset_index(drop=True)
    d["_grp"] = d["p4_gene_group"].astype(str)
    sizes = d.groupby("_grp").size().to_dict()
    prev = d.groupby("_grp")["y"].mean().to_dict()
    inner_k = max(2, round(1 / INNER_VALID_FRAC))
    chunks = []
    for seed in SEEDS:
        outer = balanced_group_folds(sizes, prev, K_OUTER, seed)
        for out_fold in range(K_OUTER):
            train_groups = [g for g, f in outer.items() if f != out_fold]
            tsz = {g: sizes[g] for g in train_groups}
            tpv = {g: prev[g] for g in train_groups}
            ia = balanced_group_folds(tsz, tpv, inner_k, seed * 97 + out_fold)
            inner = {g: (1 if ia[g] == 0 else 0) for g in train_groups}
            res = _fold_shap(d, outer, inner, out_fold, seed)
            if res is not None:
                chunks.append(res)
    return pd.concat(chunks, ignore_index=True)


def per_variant(long: pd.DataFrame) -> pd.DataFrame:
    """Average each variant's SHAP values and feature values over the 10 seeds."""
    shap_cols = [f"shap__{f}" for f in FEATS]
    val_cols = [f"val__{f}" for f in FEATS]
    return long.groupby("VariationID").agg(
        {**{c: "mean" for c in shap_cols + val_cols + ["prob", "shap_base", "y"]},
         "consequence_class": "first", "review_star": "first"}).reset_index()


def importance_table(agg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for f in FEATS:
        s = agg[f"shap__{f}"]
        corr = np.corrcoef(agg[f"val__{f}"], s)[0, 1]
        rows.append({"feature": f, "group": GROUPS[f], "mean_abs_shap": float(s.abs().mean()),
                     "mean_shap": float(s.mean()),
                     "direction": "higher value -> " + ("pathogenic" if corr > 0 else "benign")
                     if np.isfinite(corr) else "n/a (constant)"})
    t = pd.DataFrame(rows).sort_values("mean_abs_shap", ascending=False)
    t["share_pct"] = (100 * t["mean_abs_shap"] / t["mean_abs_shap"].sum()).round(1)
    return t


def group_importance(agg: pd.DataFrame) -> pd.DataFrame:
    t = importance_table(agg)
    g = t.groupby("group")["mean_abs_shap"].sum().sort_values(ascending=False)
    return (g / g.sum() * 100).round(1).rename("share_pct").reset_index()


def group_shap_share(sub: pd.DataFrame) -> dict:
    g = {gr: 0.0 for gr in set(GROUPS.values())}
    for f in FEATS:
        g[GROUPS[f]] += sub[f"shap__{f}"].abs().mean()
    tot = sum(g.values()) or 1.0
    return {gr: round(100 * v / tot, 1) for gr, v in g.items()}


def by_consequence(agg: pd.DataFrame) -> pd.DataFrame:
    rows = [{"consequence_class": cc, "n": len(sub), **group_shap_share(sub)}
            for cc, sub in agg.groupby("consequence_class")]
    return pd.DataFrame(rows)


def conservation_bands(agg: pd.DataFrame) -> pd.DataFrame:
    """Feature-group |SHAP| share in the low / mid / high conservation tercile —
    tests whether other feature groups compensate in an intermediate band."""
    c = agg["val__ensembl_conservation"]
    lo, hi = c.quantile(0.33), c.quantile(0.66)
    bands = {"low_conservation": c <= lo, "mid_conservation": (c > lo) & (c < hi),
             "high_conservation": c >= hi}
    rows = []
    for name, mask in bands.items():
        rows.append({"band": name, "n": int(mask.sum()),
                     "conservation_range": f"[{c[mask].min():.2f}, {c[mask].max():.2f}]",
                     **group_shap_share(agg[mask])})
    return pd.DataFrame(rows)
