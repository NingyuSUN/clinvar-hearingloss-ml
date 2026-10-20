"""Paired whole-gene-group bootstrap for comparing two models' AUCs.

Both models' seed-averaged OOF scores are evaluated on the *same* resampled
variants in every draw, so the interval is for the AUC difference itself. The
resampling scheme matches ``metrics.bootstrap_ci``: sample G biological gene
groups with replacement and keep every variant of each sampled group. Intervals
are conditional on the fitted models and their fixed predictions; they do not
capture refitting or all cross-validation dependence.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _weighted_auc(y: np.ndarray, p: np.ndarray):
    """Return ``w -> AUC`` with cluster multiplicities as weights; tied scores
    receive half credit, as in ``metrics.bootstrap_ci``."""
    order = np.argsort(p, kind="stable")
    ys = y[order]
    starts = np.r_[0, np.flatnonzero(np.diff(p[order])) + 1]

    def auc(w: np.ndarray) -> float:
        ws = w[order]
        pos = np.add.reduceat(ws * ys, starts)
        neg = np.add.reduceat(ws * (1 - ys), starts)
        npos, nneg = pos.sum(), neg.sum()
        if not npos or not nneg:
            return float("nan")
        return float(np.sum(pos * (np.cumsum(neg) - 0.5 * neg)) / (npos * nneg))

    return auc


def _summary(point: float, values: np.ndarray, n_boot: int) -> dict:
    enough = len(values) >= max(100, int(np.ceil(0.9 * n_boot)))
    return {
        "estimate": point if np.isfinite(point) else None,
        "bootstrap_mean": float(values.mean()) if len(values) else None,
        "lo": float(np.percentile(values, 2.5)) if enough else None,
        "hi": float(np.percentile(values, 97.5)) if enough else None,
        "n_valid": int(len(values)), "n_undefined": int(n_boot - len(values)),
        "status": "ok" if enough else "insufficient_valid_draws",
    }


def paired_group_bootstrap_auc(y_true, score_a, score_b, *, groups, n_boot: int = 2000,
                               seed: int = 0) -> dict:
    """Empirical AUC(b) - AUC(a) with a paired gene-group percentile interval.

    ``superior`` is True only when the interval's lower bound is above zero (the
    prespecified rule); a lower bound at or below zero is "not demonstrated",
    never evidence of equivalence.
    """
    y = np.asarray(y_true)
    a = np.asarray(score_a, dtype=float)
    b = np.asarray(score_b, dtype=float)
    g = np.asarray(groups)
    if y.ndim != 1 or not len(y) or any(v.shape != y.shape for v in (a, b, g)):
        raise ValueError("y, both scores and groups must be nonempty aligned vectors")
    if not np.isin(y, [0, 1]).all() or not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError("labels must be binary and scores finite")
    if pd.isna(g).any() or any(not str(v).strip() for v in g):
        raise ValueError("biological gene groups must be present")
    if all(str(v).startswith("row") and str(v)[3:].isdigit() for v in g):
        raise ValueError("row split identifiers are not biological gene groups")
    y = y.astype(int)
    if np.unique(y).size != 2:
        raise ValueError("original sample must contain both label classes")
    codes, names = pd.factorize(g, sort=True)
    ng = len(names)
    if ng < 2 or not isinstance(n_boot, (int, np.integer)) or n_boot < 100:
        raise ValueError("need at least two groups and 100 bootstrap draws")

    auc_a, auc_b = _weighted_auc(y, a), _weighted_auc(y, b)
    ones = np.ones(len(y))
    pa, pb = auc_a(ones), auc_b(ones)
    da, db, dd = [], [], []
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        w = np.bincount(rng.integers(0, ng, ng), minlength=ng)[codes].astype(float)
        va, vb = auc_a(w), auc_b(w)
        if np.isfinite(va) and np.isfinite(vb):
            da.append(va)
            db.append(vb)
            dd.append(vb - va)
    da, db, dd = map(np.asarray, (da, db, dd))
    delta = _summary(pb - pa, dd, n_boot)
    delta["fraction_draws_le_zero"] = float((dd <= 0).mean()) if len(dd) else None
    return {
        "auc_a": _summary(pa, da, n_boot),
        "auc_b": _summary(pb, db, n_boot),
        "delta_b_minus_a": delta,
        "superior": bool(delta["status"] == "ok" and delta["lo"] > 0),
        "metadata": {
            "method": "paired_gene_group_cluster_percentile", "confidence_level": 0.95,
            "n_groups": ng, "n_variants": int(len(y)), "n_positive": int(y.sum()),
            "n_boot": n_boot, "seed": seed,
            "estimand": "difference_in_variant_weighted_AUC_of_fixed_seed_averaged_oof_predictions",
            "scope": "conditional_on_fitted_models_and_observed_groups; no_refitting",
            "rule": "superior iff 2.5th percentile of the paired delta > 0",
        },
    }
