"""Gene-held-out cross-validation.

  * outer   : grouped K-fold on the gene groups, balanced on row count and label
              prevalence
  * inner   : a balanced gene-group hold-out inside each training fold, for
              iteration selection and the R90 threshold (never rows, never the
              outer-test genes)
  * model   : XGBoost trained to a cap with no in-loop early stopping; the best
              iteration is the argmin of the inner-validation log-loss curve
  * threshold : R90 chosen on that inner-validation, frozen, applied to the outer
              test using the same fitted model
  * impute  : continuous features -> inner-train median; binary indicators -> 0,
              recomputed inside every split
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score

from .config import (INNER_VALID_FRAC, K_OUTER, NUM_BOOST_ROUND, POSTHOC_MIN_ITER,
                     SEEDS, xgb_params)
from .metrics import metrics_at, r90_threshold
from .protocol import CONTINUOUS, FEATURE_SETS


def balanced_group_folds(sizes: dict, prev: dict, k: int, seed: int) -> dict:
    """Assign gene groups to ``k`` folds, greedily balancing total rows and pooled
    label prevalence (a longest-processing-time pass plus one refinement pass).
    Deterministic given ``seed``."""
    rng = np.random.default_rng(seed)
    groups = list(sizes)
    rng.shuffle(groups)
    groups.sort(key=lambda g: sizes[g], reverse=True)
    total = max(1, sum(sizes.values()))
    target = total / k
    overall_prev = sum(prev[g] * sizes[g] for g in groups) / total
    fold_rows, fold_pos, fold_tot = [0] * k, [0.0] * k, [0] * k
    assign: dict[str, int] = {}
    for g in groups:
        s, p = sizes[g], prev[g]
        best, best_cost = 0, None
        for f in range(k):
            after_tot = fold_tot[f] + s
            prev_after = (fold_pos[f] + p * s) / after_tot
            cost = (fold_rows[f] + s) / target + 3.0 * abs(prev_after - overall_prev)
            if best_cost is None or cost < best_cost:
                best, best_cost = f, cost
        assign[g] = best
        fold_rows[best] += s
        fold_pos[best] += p * s
        fold_tot[best] += s
    for _ in range(50):
        hi, lo = int(np.argmax(fold_rows)), int(np.argmin(fold_rows))
        if fold_rows[hi] - fold_rows[lo] <= target * 0.06:
            break
        moved = False
        for g in sorted((g for g, f in assign.items() if f == hi), key=lambda g: sizes[g]):
            if sizes[g] <= fold_rows[hi] - fold_rows[lo]:
                assign[g] = lo
                for arr, val in ((fold_rows, sizes[g]), (fold_pos, prev[g] * sizes[g]),
                                 (fold_tot, sizes[g])):
                    arr[hi] -= val
                    arr[lo] += val
                moved = True
                break
        if not moved:
            break
    return assign


def _impute(train_x: pd.DataFrame, *others: pd.DataFrame):
    med = {c: float(train_x[c].median()) for c in train_x.columns if c in CONTINUOUS}

    def fix(df):
        df = df.copy()
        for c in df.columns:
            df[c] = df[c].fillna(med.get(c, 0))
        return df

    return (fix(train_x), *[fix(o) for o in others])


def _one_fold(data, feats, outer, inner, out_fold, seed):
    te = data[data["_grp"].map(outer) == out_fold]
    tr = data[data["_grp"].map(outer) != out_fold]
    iv = tr[tr["_grp"].map(inner) == 1]
    it = tr[tr["_grp"].map(inner) == 0]
    if not len(iv) or not len(it):
        return None
    xit, xiv, xte = _impute(it[feats], iv[feats], te[feats])
    yit, yiv, yte = (it["y"].astype(int).to_numpy(), iv["y"].astype(int).to_numpy(),
                     te["y"].astype(int).to_numpy())
    hist: dict = {}
    dit = xgb.DMatrix(xit, label=yit, feature_names=feats)
    div = xgb.DMatrix(xiv, label=yiv, feature_names=feats)
    dte = xgb.DMatrix(xte, label=yte, feature_names=feats)
    booster = xgb.train(xgb_params(seed), dit, num_boost_round=NUM_BOOST_ROUND,
                        evals=[(div, "iv")], evals_result=hist, verbose_eval=False)
    logloss = np.asarray(hist["iv"]["logloss"], dtype=float)
    if len(logloss) > POSTHOC_MIN_ITER:
        best_iter = int(np.argmin(logloss[POSTHOC_MIN_ITER:]) + POSTHOC_MIN_ITER)
    else:
        best_iter = int(np.argmin(logloss))
    piv = booster.predict(div, iteration_range=(0, best_iter + 1))
    pte = booster.predict(dte, iteration_range=(0, best_iter + 1))
    thr, ok, reason = r90_threshold(yiv, piv)

    row = {"seed": seed, "outer_fold": out_fold, "n_test": len(te), "best_iteration": best_iter,
           "threshold_ok": ok, "threshold_reason": reason,
           "auc": float(roc_auc_score(yte, pte)) if len(np.unique(yte)) > 1 else np.nan}
    if ok and thr is not None:
        row["threshold"] = float(thr)
        for key, val in metrics_at(yte, pte, thr).items():
            row[f"R90_{key}"] = val
    oof = pd.DataFrame({
        "VariationID": te["VariationID"].to_numpy(), "y": yte, "prob": pte,
        "seed": seed, "outer_fold": out_fold,
        "GeneSymbol": te["GeneSymbol"].to_numpy(),
        "gene_group": te["_grp"].to_numpy(),
        "consequence_class": te["consequence_class"].to_numpy(),
        "review_star": te["p4_review_star"].to_numpy(),
    })
    return row, oof


def run_experiment(matrix: pd.DataFrame, cohort_flag: str, feature_sets: list[str],
                   seeds: list[int] = SEEDS, random_split: bool = False):
    """Run every ``feature_set`` over ``seeds`` x ``K_OUTER`` folds on the cohort
    selected by ``cohort_flag``. Returns ``(per_fold_df, {feature_set: oof_df})``."""
    d = matrix[matrix[cohort_flag] == 1].reset_index(drop=True).copy()
    d["_grp"] = ([f"row{i}" for i in range(len(d))] if random_split
                 else d["p4_gene_group"].astype(str))
    sizes = d.groupby("_grp").size().to_dict()
    prev = d.groupby("_grp")["y"].mean().to_dict()

    rows: list = []
    oof: dict = {}
    for seed in seeds:
        outer = balanced_group_folds(sizes, prev, K_OUTER, seed)
        for out_fold in range(K_OUTER):
            train_groups = [g for g, f in outer.items() if f != out_fold]
            tsz = {g: sizes[g] for g in train_groups}
            tpv = {g: prev[g] for g in train_groups}
            inner_k = max(2, round(1 / INNER_VALID_FRAC))
            inner_folds = balanced_group_folds(tsz, tpv, inner_k, seed * 97 + out_fold)
            inner = {g: (1 if inner_folds[g] == 0 else 0) for g in train_groups}
            for name in feature_sets:
                res = _one_fold(d, FEATURE_SETS[name], outer, inner, out_fold, seed)
                if res is None:
                    continue
                fold_row, fold_oof = res
                fold_row |= {"feature_set": name, "n_features": len(FEATURE_SETS[name]),
                             "cohort": cohort_flag}
                rows.append(fold_row)
                oof.setdefault(name, []).append(fold_oof)
    per_fold = pd.DataFrame(rows)
    return per_fold, {name: pd.concat(chunks, ignore_index=True) for name, chunks in oof.items()}
