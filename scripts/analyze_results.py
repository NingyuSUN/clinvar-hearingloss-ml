#!/usr/bin/env python3
"""Build the reported tables from the per-fold output of run_evaluation.py.

    python scripts/analyze_results.py [--dir results/] [--boot 2000]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from hlpath.analysis import (ladder_table, paired_ladder_table, pooled_oof_table,
                             split_gap, stratified_table)
from hlpath.config import RESULTS_DIR


def _load(dir_: Path, tag: str):
    per_fold = pd.read_csv(dir_ / f"per_fold_{tag}.csv")
    oof = {p.stem.split("__", 1)[1]: pd.read_csv(p)
           for p in dir_.glob(f"oof_{tag}__*.csv")}
    return per_fold, oof


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=RESULTS_DIR)
    ap.add_argument("--boot", type=int, default=2000)
    args = ap.parse_args()
    d = args.dir

    summary: dict = {}
    for tag in ["headline", "headline_missense"]:
        per_fold, oof = _load(d, tag)
        ladder_table(per_fold).to_csv(d / f"ladder_{tag}.csv", index=False)
        paired_ladder_table(per_fold).to_csv(d / f"paired_ladder_{tag}.csv", index=False)
        summary[f"pooled_{tag}"] = pooled_oof_table(per_fold, oof, "L6_full", args.boot)
        if tag == "headline":
            stratified_table(per_fold, oof, "L6_full").to_csv(d / "stratified_headline.csv",
                                                             index=False)

    gene_pf, gene_oof = _load(d, "headline")
    rnd_pf, rnd_oof = _load(d, "headline_randomsplit")
    summary["split_gap"] = split_gap(gene_pf, gene_oof, rnd_pf, rnd_oof, "L6_full")

    sens_rows = []
    for tag in ["sensitivity_alltier", "headline", "sensitivity_ge2star"]:
        per_fold, oof = _load(d, tag)
        p = pooled_oof_table(per_fold, oof, "L6_full", args.boot)
        sens_rows.append({"cohort": tag, "n_variants": p["n_variants"],
                          **{f"AUC_{k}": round(v, 4) for k, v in p["bootstrap_ci"]["AUC"].items()}})
    pd.DataFrame(sens_rows).to_csv(d / "review_status_sensitivity.csv", index=False)
    summary["review_status_sensitivity"] = sens_rows

    (d / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
