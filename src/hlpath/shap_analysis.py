"""SHAP utilities for model interpretation."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import shap
import xgboost as xgb


def tree_shap_values(booster: xgb.Booster, X: pd.DataFrame):
    explainer = shap.TreeExplainer(booster)
    values = explainer.shap_values(X)
    return explainer, values


def save_summary_plot(shap_values, X: pd.DataFrame, out_path: str | Path, max_display: int = 20):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    shap.summary_plot(shap_values, X, max_display=max_display, show=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def save_dependence_plot(shap_values, X: pd.DataFrame, feature: str, out_path: str | Path, interaction_index=None):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    shap.dependence_plot(feature, shap_values, X, interaction_index=interaction_index, show=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def gray_zone(X: pd.DataFrame, col: str = "ensembl_conservation", low: float = -1, high: float = 1) -> pd.DataFrame:
    return X[(X[col] >= low) & (X[col] <= high)].copy()
