"""Paths, seeds, and cross-validation constants.

The data directory can be overridden with the ``HLPATH_DATA`` environment variable;
by default it is ``<repo>/data``.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("HLPATH_DATA", REPO_ROOT / "data"))
RESULTS_DIR = Path(os.environ.get("HLPATH_RESULTS", REPO_ROOT / "results"))

# The per-variant modelling matrix (features + label + gene group + cohort flags).
MODELING_MATRIX = DATA_DIR / "modeling_matrix.csv"

# Cross-validation.
K_OUTER = 5                       # outer grouped folds
INNER_VALID_FRAC = 0.20          # gene groups held out inside each training fold
SEEDS = list(range(101, 111))    # 10 frozen seeds; each is one balanced fold assignment
TARGET_RECALL = 0.90            # the R90 operating point

# XGBoost.
NUM_BOOST_ROUND = 700           # hard cap; best iteration is chosen post hoc
POSTHOC_MIN_ITER = 10


def xgb_params(seed: int) -> dict:
    return {
        "objective": "binary:logistic",
        "eval_metric": ["auc", "logloss"],
        "eta": 0.03,
        "max_depth": 5,
        "subsample": 0.9,
        "colsample_bytree": 0.9,
        "lambda": 1.0,
        "min_child_weight": 1,
        "tree_method": "hist",
        "seed": seed,
        "nthread": int(os.environ.get("HLPATH_NTHREAD", "8")),
    }
