#!/usr/bin/env python3
"""Run gene-held-out cross-validation for every feature set on the headline cohort
and the missense subset, plus the review-status sensitivity cohorts and the
random-split leakage diagnostic.

    python scripts/run_evaluation.py [--boot 2000] [--out results/]

Writes per-fold tables, out-of-fold predictions, and a run manifest.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import xgboost as xgb

from hlpath.config import K_OUTER, MODELING_MATRIX, NUM_BOOST_ROUND, RESULTS_DIR, SEEDS, xgb_params
from hlpath.features import load_matrix
from hlpath.pipeline import run_experiment
from hlpath.protocol import FEATURE_SETS

SENSITIVITY_SETS = ["L2_base", "L6_full"]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=RESULTS_DIR)
    ap.add_argument("--matrix", type=Path, default=None)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    matrix = load_matrix(args.matrix)
    tasks = {
        "headline":            ("p4_headline_cohort",   list(FEATURE_SETS), False, None),
        "headline_missense":   ("p4_headline_cohort",   list(FEATURE_SETS), False, ["coding_nontruncating"]),
        "headline_randomsplit":("p4_headline_cohort",   ["L2_base", "L6_full"], True, None),
        "sensitivity_alltier": ("p4_in_modeling_cohort", SENSITIVITY_SETS, False, None),
        "sensitivity_ge2star": ("p4_sens_ge2star",       SENSITIVITY_SETS, False, None),
    }

    for tag, (flag, sets, random_split, consequence) in tasks.items():
        sub = matrix
        if consequence is not None:
            sub = matrix[matrix["consequence_class"].isin(consequence)].reset_index(drop=True)
        per_fold, oof = run_experiment(sub, flag, sets, seeds=SEEDS, random_split=random_split)
        per_fold.to_csv(args.out / f"per_fold_{tag}.csv", index=False)
        for name, frame in oof.items():
            frame.to_csv(args.out / f"oof_{tag}__{name}.csv", index=False)
        print(f"{tag}: {len(per_fold)} fold-runs, {len(oof)} feature sets")

    matrix_path = args.matrix or MODELING_MATRIX
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "k_outer": K_OUTER, "seeds": SEEDS, "num_boost_round_cap": NUM_BOOST_ROUND,
        "model_selection": "argmin of the inner-validation log-loss curve; no in-loop early stopping",
        "xgb_params": xgb_params(SEEDS[0]),
        "xgboost_version": xgb.__version__,
        "matrix": str(matrix_path),
        "matrix_sha256": _sha(matrix_path) if matrix_path.exists() else None,
    }
    (args.out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {args.out}/run_manifest.json")


if __name__ == "__main__":
    main()
