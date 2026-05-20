"""Evaluation metrics for pathogenicity prediction."""
from __future__ import annotations

import numpy as np
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score


def metrics_at_threshold(y_true, y_prob, threshold: float) -> dict:
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    y_pred = (y_prob >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    return {
        "threshold": float(threshold),
        "TN": int(tn),
        "FP": int(fp),
        "FN": int(fn),
        "TP": int(tp),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "recall": float(recall),
        "precision": float(precision),
    }


def r90_metrics(y_true, y_prob, target_recall: float = 0.90) -> dict:
    """Choose a threshold with recall >= target_recall and minimal FP.

    If several thresholds have the same FP count, the higher threshold is preferred.
    """
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    best = None
    for thr in np.unique(y_prob):
        m = metrics_at_threshold(y_true, y_prob, thr)
        if m["recall"] >= target_recall:
            cand = {
                "thr": m["threshold"],
                "FN": m["FN"],
                "FP": m["FP"],
                "Recall": m["recall"],
                "Precision": m["precision"],
            }
            if best is None or (cand["FP"], -cand["thr"]) < (best["FP"], -best["thr"]):
                best = cand
    if best is None:
        prevalence = float(np.mean(y_true))
        best = {
            "thr": float("-inf"),
            "FN": int(np.sum(y_true == 1)),
            "FP": int(np.sum(y_true == 0)),
            "Recall": 1.0,
            "Precision": prevalence,
        }
    return best


def evaluate_auc_r90(y_true, y_prob, target_recall: float = 0.90) -> dict:
    return {
        "AUC": float(roc_auc_score(y_true, y_prob)),
        **{f"R90_{k}": v for k, v in r90_metrics(y_true, y_prob, target_recall).items()},
    }
