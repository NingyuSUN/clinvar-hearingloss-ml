#!/usr/bin/env python3
"""Fit the PLM-ESM-2026-10 feature sets under the unchanged gene-held-out protocol.

    HLPATH_NTHREAD=10 python scripts/run_plm_evaluation.py --out runs/plm_esm_2026-10

Two tasks, the same seeds (101-110), outer folds, inner split, XGBoost settings,
iteration selection, R90 rule and imputation as the published evaluation:

  * ``plm_codingnt``  headline cohort restricted to coding non-truncating
                      (trained and evaluated within the subset, like the
                      published 0.815 model)
  * ``plm_headline``  the full 7,125-variant headline cohort

Writes per-fold tables and OOF predictions, then ``prediction_seal.json`` with
the SHA256 of every prediction file. No metric summary is printed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

from hlpath.config import K_OUTER, MODELING_MATRIX, NUM_BOOST_ROUND, REPO_ROOT, SEEDS, xgb_params
from hlpath.features import load_matrix
from hlpath.pipeline import run_experiment
from hlpath.plm import PLM_CONTINUOUS, PLM_FEATURE_SETS, PROTOCOL_ID, attach_plm_features

TASKS = {
    "plm_codingnt": ("p4_headline_cohort", ["coding_nontruncating"]),
    "plm_headline": ("p4_headline_cohort", None),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, check=True,
                              capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--features", type=Path, default=REPO_ROOT / "data/plm/plm_features.csv")
    ap.add_argument("--tasks", nargs="*", default=list(TASKS), choices=list(TASKS))
    args = ap.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit(f"{args.out} must be new or empty")
    args.out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    matrix_path = MODELING_MATRIX.with_suffix(".csv.gz")
    matrix = attach_plm_features(load_matrix(), pd.read_csv(args.features))
    files = {}
    for tag in args.tasks:
        flag, classes = TASKS[tag]
        sub = matrix if classes is None else matrix[matrix["consequence_class"].isin(classes)]
        per_fold, oof = run_experiment(sub.reset_index(drop=True), flag, list(PLM_FEATURE_SETS),
                                       seeds=SEEDS, feature_map=PLM_FEATURE_SETS,
                                       continuous=PLM_CONTINUOUS)
        path = args.out / f"per_fold_{tag}.csv"
        per_fold.to_csv(path, index=False, lineterminator="\n")
        files[path.name] = sha256(path)
        for name, frame in oof.items():
            path = args.out / f"oof_{tag}__{name}.csv"
            frame.to_csv(path, index=False, lineterminator="\n")
            files[path.name] = sha256(path)
        print(f"{tag}: {len(per_fold)} fold-runs, {len(oof)} feature sets", flush=True)

    seal = {
        "protocol": PROTOCOL_ID,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(),
        "tasks": {t: {"cohort_flag": TASKS[t][0], "consequence_classes": TASKS[t][1]} for t in args.tasks},
        "feature_sets": PLM_FEATURE_SETS,
        "k_outer": K_OUTER, "seeds": SEEDS, "num_boost_round_cap": NUM_BOOST_ROUND,
        "xgb_params": xgb_params(SEEDS[0]),
        "versions": {"xgboost": xgb.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        "inputs_sha256": {"data/modeling_matrix.csv.gz": sha256(matrix_path),
                          "data/plm/plm_features.csv": sha256(args.features)},
        "predictions_sha256": files,
        "seconds": round(time.time() - t0, 1),
    }
    (args.out / "prediction_seal.json").write_text(json.dumps(seal, indent=2) + "\n")
    print(f"sealed {len(files)} prediction files in {args.out / 'prediction_seal.json'}")


if __name__ == "__main__":
    main()
