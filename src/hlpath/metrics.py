"""Evaluation metrics: the R90 operating point, confusion at a threshold, a
gene-group conditional bootstrap interval, and descriptive paired differences.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import TARGET_RECALL


def metrics_at(y_true, y_prob, threshold: float) -> dict:
    y = np.asarray(y_true).astype(int)
    pred = (np.asarray(y_prob, dtype=float) >= threshold).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    return {
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "recall": tp / (tp + fn) if (tp + fn) else float("nan"),
        "precision": tp / (tp + fp) if (tp + fp) else float("nan"),
        "specificity": tn / (tn + fp) if (tn + fp) else float("nan"),
        "accuracy": (tp + tn) / len(y),
    }


def r90_threshold(y_true, y_prob, target: float = TARGET_RECALL):
    """Smallest-FP threshold with recall >= ``target`` on the validation scores.

    Returns ``(threshold, ok, reason)``. ``ok`` is False, with no fabricated
    recall = 1, when the validation set is single-class or the scores are not finite.
    """
    y = np.asarray(y_true).astype(int)
    p = np.asarray(y_prob, dtype=float)
    if not np.isfinite(p).all():
        return None, False, "non_finite_scores"
    if len(np.unique(y)) < 2:
        return None, False, "one_class_validation"
    n_pos = int((y == 1).sum())
    best = None
    for thr in np.unique(p):
        pred = (p >= thr).astype(int)
        tp = int(((pred == 1) & (y == 1)).sum())
        fp = int(((pred == 1) & (y == 0)).sum())
        if tp / n_pos >= target:
            key = (fp, -thr)
            if best is None or key < best[0]:
                best = (key, float(thr))
    if best is None:
        return float(p.min()), True, "target_unreachable_use_min_threshold"
    return best[1], True, "ok"


def _rank_auc(y, s) -> float:
    """Mann-Whitney AUC, fast for repeated bootstrap draws."""
    y = np.asarray(y)
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = pd.Series(np.asarray(s, dtype=float)).rank().to_numpy()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def bootstrap_ci(y_true, prob, pred=None, n_boot: int = 2000, seed: int = 0,
                 *, groups) -> dict:
    """Gene-group percentile intervals conditional on fixed OOF predictions.

    Sample G whole biological groups with replacement. Scores and thresholded
    predictions stay fixed; seeds are not independent observations. This does
    not estimate retraining, selection, or all cross-fold prediction dependence.
    Metrics remain variant-weighted (not a mean of within-gene metrics).
    """
    y = np.asarray(y_true)
    p = np.asarray(prob, dtype=float)
    g = np.asarray(groups)
    if y.ndim != 1 or p.shape != y.shape or g.shape != y.shape or not len(y):
        raise ValueError("y, prob and groups must be nonempty aligned vectors")
    if not np.isin(y, [0, 1]).all() or not np.isfinite(p).all():
        raise ValueError("labels must be binary and scores finite")
    if pd.isna(g).any() or any(not str(v).strip() for v in g):
        raise ValueError("biological gene groups must be present")
    if all(str(v).startswith("row") and str(v)[3:].isdigit() for v in g):
        raise ValueError("row split identifiers are not biological gene groups")
    codes, names = pd.factorize(g, sort=True)
    ng = len(names)
    if ng < 2 or not isinstance(n_boot, (int, np.integer)) or n_boot < 100:
        raise ValueError("need at least two groups and 100 bootstrap draws")
    y = y.astype(int)
    if np.unique(y).size != 2:
        raise ValueError("original sample must contain both label classes")
    pr = None if pred is None else np.asarray(pred)
    if pr is not None and (pr.shape != y.shape or not np.isin(pr, [0, 1]).all()):
        raise ValueError("pred must be an aligned binary vector")

    # Sort once; cluster multiplicities are exact observation weights. Equal
    # score bins retain the half-credit convention for tied case-control pairs.
    order = np.argsort(p, kind="stable")
    ys = y[order]
    starts = np.r_[0, np.flatnonzero(np.diff(p[order])) + 1]

    def measures(w):
        ws = w[order]
        pos = np.add.reduceat(ws * ys, starts)
        neg = np.add.reduceat(ws * (1 - ys), starts)
        npos, nneg = pos.sum(), neg.sum()
        auc = (np.sum(pos * (np.cumsum(neg) - 0.5 * neg)) / (npos * nneg)
               if npos and nneg else np.nan)
        values = {"AUC": float(auc)}
        if pr is not None:
            tp = np.sum(w * y * pr)
            called = np.sum(w * pr)
            values.update(R90_recall=float(tp / npos) if npos else np.nan,
                          R90_precision=float(tp / called) if called else np.nan)
        return values

    point = measures(np.ones(len(y)))
    draws = {key: [] for key in point}
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        counts = np.bincount(rng.integers(0, ng, ng), minlength=ng)
        for key, val in measures(counts[codes]).items():
            if np.isfinite(val):
                draws[key].append(val)
    out = {}
    for key, estimate in point.items():
        values = np.asarray(draws[key])
        enough = len(values) >= max(100, int(np.ceil(0.9 * n_boot)))
        out[key] = {
            "estimate": estimate if np.isfinite(estimate) else None,
            "bootstrap_mean": float(values.mean()) if len(values) else None,
            "lo": float(np.percentile(values, 2.5)) if enough else None,
            "hi": float(np.percentile(values, 97.5)) if enough else None,
            "n_valid": len(values), "n_undefined": n_boot - len(values),
            "status": "ok" if enough else "insufficient_valid_draws",
        }
    out["metadata"] = {
        "method": "gene_group_cluster_percentile", "confidence_level": 0.95,
        "n_groups": ng, "n_variants": len(y), "n_boot": n_boot, "seed": seed,
        "estimand": "variant_weighted_metrics_of_fixed_seed_averaged_oof_predictions",
        "scope": "conditional_on_fitted_models_and_observed_groups; no_refitting",
        "undefined_draws": "AUC: missing label class; recall: no positives; precision: no calls",
    }
    return out


def paired_delta(per_fold: pd.DataFrame, a: str, b: str, value: str = "auc") -> dict:
    """Descriptive b-a on matched folds; no independent-fold SE or t test.

    Fold SD and the range of seed means describe observed variability, not a
    confidence interval. Repeated seeds reuse subjects and fitted training data.
    """
    selected = per_fold[per_fold["feature_set"].isin([a, b])]
    if selected.empty:
        return {}
    if set(selected["feature_set"]) != {a, b}:
        raise ValueError(f"missing model in comparison {b} - {a}")
    keys = ["seed", "outer_fold"]
    if selected[keys].isna().any().any():
        raise ValueError("pairing keys must be present")
    if selected.duplicated(keys + ["feature_set"]).any():
        raise ValueError("duplicate model/fold rows")
    wide = selected.pivot(index=keys, columns="feature_set", values=value)
    presence = selected.assign(present=1).pivot(index=keys, columns="feature_set", values="present")
    if presence.isna().any().any():
        raise ValueError("model comparisons must have matching folds")
    valid = np.isfinite(wide[a]) & np.isfinite(wide[b])
    d = (wide.loc[valid, b] - wide.loc[valid, a])
    seed_means = d.groupby(level="seed").mean()
    return {
        "contrast": f"{b} - {a}", "inference": "descriptive_only",
        "n_folds": len(wide), "n_pairs": len(d), "n_undefined_pairs": int((~valid).sum()),
        "mean_delta": float(d.mean()) if len(d) else None,
        "sd_delta": float(d.std(ddof=1)) if len(d) > 1 else None,
        "n_seeds": len(seed_means),
        "seed_mean_min": float(seed_means.min()) if len(d) else None,
        "seed_mean_max": float(seed_means.max()) if len(d) else None,
    }
