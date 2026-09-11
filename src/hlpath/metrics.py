"""Evaluation metrics: the R90 operating point, confusion at a threshold, a
patient-level bootstrap CI, and paired per-fold comparisons.
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


def bootstrap_ci(y_true, prob, pred=None, n_boot: int = 2000, seed: int = 0) -> dict:
    """Patient-level bootstrap 95% CI on pooled out-of-fold predictions.

    The mean +/- SD across overlapping folds is not a confidence interval
    (Nadeau & Bengio 2003); resampling variants is.
    """
    y = np.asarray(y_true).astype(int)
    p = np.asarray(prob, dtype=float)
    rng = np.random.default_rng(seed)
    n = len(y)
    aucs, recs, precs = [], [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yb = y[idx]
        if yb.sum() in (0, len(yb)):
            continue
        aucs.append(_rank_auc(yb, p[idx]))
        if pred is not None:
            prb = np.asarray(pred)[idx]
            tp = int(((prb == 1) & (yb == 1)).sum())
            fp = int(((prb == 1) & (yb == 0)).sum())
            fn = int(((prb == 0) & (yb == 1)).sum())
            recs.append(tp / (tp + fn) if (tp + fn) else np.nan)
            precs.append(tp / (tp + fp) if (tp + fp) else np.nan)

    def ci(values):
        v = np.asarray(values, dtype=float)
        v = v[np.isfinite(v)]
        if not len(v):
            return None
        return {"mean": float(v.mean()),
                "lo": float(np.percentile(v, 2.5)),
                "hi": float(np.percentile(v, 97.5))}

    out = {"AUC": ci(aucs)}
    if pred is not None:
        out["R90_recall"] = ci(recs)
        out["R90_precision"] = ci(precs)
    return out


def paired_delta(per_fold: pd.DataFrame, a: str, b: str, value: str = "auc") -> dict:
    """Paired per-fold difference b - a with a paired-t statistic.

    Comparing a mean delta to the marginal SD is invalid for paired cross-validation;
    the paired difference has a much smaller standard error.
    """
    wide = per_fold.pivot_table(index=["seed", "outer_fold"], columns="feature_set", values=value)
    if a not in wide or b not in wide:
        return {}
    d = (wide[b] - wide[a]).dropna()
    se = d.std(ddof=1) / np.sqrt(len(d))
    return {
        "contrast": f"{b} - {a}",
        "n_folds": int(len(d)),
        "mean_delta": float(d.mean()),
        "paired_se": float(se),
        "paired_t": float(d.mean() / se) if se else float("nan"),
        "significant_p05": bool(abs(d.mean() / se) > 2.01) if se else False,
    }
