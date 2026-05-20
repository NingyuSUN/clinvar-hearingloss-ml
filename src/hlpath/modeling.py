"""Model training utilities for strict gene-based evaluation."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

from .metrics import r90_metrics


@dataclass
class GeneSplitResult:
    experiment: str
    auc: float
    r90: dict
    best_iteration: int
    booster: xgb.Booster
    test_pred: pd.DataFrame
    train_genes: np.ndarray
    test_genes: np.ndarray


def default_xgb_params(seed: int = 42) -> dict:
    return {
        "objective": "binary:logistic",
        "eval_metric": "auc",
        "eta": 0.03,
        "max_depth": 5,
        "subsample": 0.9,
        "colsample_bytree": 0.9,
        "lambda": 1.0,
        "min_child_weight": 1,
        "tree_method": "hist",
        "seed": seed,
    }


def run_gene_split_xgb(
    df: pd.DataFrame,
    feature_cols: list[str],
    experiment: str = "xgboost_gene_split",
    test_size_genes: float = 0.2,
    valid_size_rows: float = 0.2,
    seed: int = 42,
    num_boost_round: int = 2000,
    early_stopping_rounds: int = 50,
) -> GeneSplitResult:
    """Train XGBoost with held-out genes and strict early stopping.

    Test genes are never used for early stopping.
    Missing values are imputed with the median from the training subset only.
    """
    d = df.dropna(subset=["label", "GeneSymbol"]).copy()
    d["label"] = d["label"].astype(int)

    genes = d["GeneSymbol"].astype(str).unique()
    train_genes, test_genes = train_test_split(genes, test_size=test_size_genes, random_state=seed)

    train_df = d[d["GeneSymbol"].astype(str).isin(train_genes)].copy().reset_index(drop=True)
    test_df = d[d["GeneSymbol"].astype(str).isin(test_genes)].copy().reset_index(drop=True)

    y_train_full = train_df["label"].astype(int).values
    tr_idx, va_idx = train_test_split(
        np.arange(len(train_df)),
        test_size=valid_size_rows,
        random_state=seed,
        stratify=y_train_full,
    )

    tr = train_df.iloc[tr_idx].copy().reset_index(drop=True)
    va = train_df.iloc[va_idx].copy().reset_index(drop=True)

    X_tr = tr[feature_cols].copy()
    y_tr = tr["label"].astype(int).values
    X_va = va[feature_cols].copy()
    y_va = va["label"].astype(int).values
    X_te = test_df[feature_cols].copy()
    y_te = test_df["label"].astype(int).values

    med = X_tr.median(numeric_only=True)
    X_tr = X_tr.fillna(med)
    X_va = X_va.fillna(med)
    X_te = X_te.fillna(med)

    dtrain = xgb.DMatrix(X_tr, label=y_tr, feature_names=feature_cols)
    dvalid = xgb.DMatrix(X_va, label=y_va, feature_names=feature_cols)
    dtest = xgb.DMatrix(X_te, label=y_te, feature_names=feature_cols)

    booster = xgb.train(
        params=default_xgb_params(seed),
        dtrain=dtrain,
        num_boost_round=num_boost_round,
        evals=[(dvalid, "valid")],
        early_stopping_rounds=early_stopping_rounds,
        verbose_eval=False,
    )

    prob = booster.predict(dtest, iteration_range=(0, booster.best_iteration + 1))
    auc = roc_auc_score(y_te, prob)
    r90 = r90_metrics(y_te, prob)

    test_pred = test_df.copy()
    test_pred["pred_prob"] = prob
    test_pred["pred_label_05"] = (prob >= 0.5).astype(int)
    test_pred["pred_label_r90"] = (prob >= r90["thr"]).astype(int)

    return GeneSplitResult(
        experiment=experiment,
        auc=float(auc),
        r90=r90,
        best_iteration=int(booster.best_iteration),
        booster=booster,
        test_pred=test_pred,
        train_genes=train_genes,
        test_genes=test_genes,
    )


def result_row(result: GeneSplitResult) -> dict:
    return {
        "Experiment": result.experiment,
        "AUC": result.auc,
        "R90_thr": result.r90["thr"],
        "R90_FN": result.r90["FN"],
        "R90_FP": result.r90["FP"],
        "R90_Recall": result.r90["Recall"],
        "R90_Precision": result.r90["Precision"],
        "best_iteration": result.best_iteration,
        "n_test_rows": len(result.test_pred),
        "n_train_genes": len(result.train_genes),
        "n_test_genes": len(result.test_genes),
    }


def compare_feature_sets(df: pd.DataFrame, feature_sets: dict[str, list[str]], seed: int = 42) -> pd.DataFrame:
    rows = []
    for name, cols in feature_sets.items():
        res = run_gene_split_xgb(df, cols, experiment=name, seed=seed)
        rows.append(result_row(res))
    return pd.DataFrame(rows)


def run_random_split_xgb(
    df: pd.DataFrame,
    feature_cols: list[str],
    experiment: str = "xgboost_random_split",
    test_size_rows: float = 0.2,
    valid_size_rows: float = 0.2,
    seed: int = 42,
    num_boost_round: int = 2000,
    early_stopping_rounds: int = 50,
) -> GeneSplitResult:
    """Train XGBoost with an ordinary stratified random row split.

    This is included only as a diagnostic baseline to compare with gene-based
    splitting. GeneSymbol is not used as a feature.
    """
    d = df.dropna(subset=["label", "GeneSymbol"]).copy()
    d["label"] = d["label"].astype(int)

    train_df, test_df = train_test_split(
        d,
        test_size=test_size_rows,
        random_state=seed,
        stratify=d["label"].astype(int),
    )
    train_df = train_df.reset_index(drop=True)
    test_df = test_df.reset_index(drop=True)

    tr, va = train_test_split(
        train_df,
        test_size=valid_size_rows,
        random_state=seed,
        stratify=train_df["label"].astype(int),
    )
    tr = tr.reset_index(drop=True)
    va = va.reset_index(drop=True)

    X_tr = tr[feature_cols].copy()
    y_tr = tr["label"].astype(int).values
    X_va = va[feature_cols].copy()
    y_va = va["label"].astype(int).values
    X_te = test_df[feature_cols].copy()
    y_te = test_df["label"].astype(int).values

    med = X_tr.median(numeric_only=True)
    X_tr = X_tr.fillna(med)
    X_va = X_va.fillna(med)
    X_te = X_te.fillna(med)

    dtrain = xgb.DMatrix(X_tr, label=y_tr, feature_names=feature_cols)
    dvalid = xgb.DMatrix(X_va, label=y_va, feature_names=feature_cols)
    dtest = xgb.DMatrix(X_te, label=y_te, feature_names=feature_cols)

    booster = xgb.train(
        params=default_xgb_params(seed),
        dtrain=dtrain,
        num_boost_round=num_boost_round,
        evals=[(dvalid, "valid")],
        early_stopping_rounds=early_stopping_rounds,
        verbose_eval=False,
    )

    prob = booster.predict(dtest, iteration_range=(0, booster.best_iteration + 1))
    auc = roc_auc_score(y_te, prob)
    r90 = r90_metrics(y_te, prob)

    test_pred = test_df.copy()
    test_pred["pred_prob"] = prob
    test_pred["pred_label_05"] = (prob >= 0.5).astype(int)
    test_pred["pred_label_r90"] = (prob >= r90["thr"]).astype(int)

    return GeneSplitResult(
        experiment=experiment,
        auc=float(auc),
        r90=r90,
        best_iteration=int(booster.best_iteration),
        booster=booster,
        test_pred=test_pred,
        train_genes=train_df["GeneSymbol"].astype(str).unique(),
        test_genes=test_df["GeneSymbol"].astype(str).unique(),
    )
