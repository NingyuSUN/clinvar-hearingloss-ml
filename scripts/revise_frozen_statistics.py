#!/usr/bin/env python3
"""Reanalyse frozen 2026-09-10 OOF files without fitting models.

Only aggregate outputs are written. The legacy CSV adapter is explicit; analysis
code always comes from this repository, not the archived implementation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from hlpath.analysis import _pooled, ladder_table, paired_ladder_table, pooled_oof_table, stratified_table
from hlpath.metrics import paired_delta
from hlpath.config import SEEDS, K_OUTER

ALIASES = {"L2_allele_conservation": "L2_base", "L3_+consequence": "L3_consequence",
           "L4_+constraint": "L4_constraint", "L5_+frequency": "L5_frequency"}
TAGS = {"headline": "headline", "coding_nontruncating": "headline_codingnt",
        "random_split": "headline_randomsplit", "all_tier": "sens_alltier", "ge2star": "sens_ge2star"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--boot", type=int, default=2000)
    args = ap.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        ap.error("output directory must be new or empty")
    args.output.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def read(name):
        path = args.source / name
        inputs[name] = sha(path)
        return pd.read_csv(path)

    summary = {"schema_version": "2.0", "analysis": "conditional_gene_group_bootstrap",
               "training_performed": False, "cohorts": {}}
    estimates, comparisons, descriptive = [], [], []
    for label, tag in TAGS.items():
        pf = read(f"results_{tag}.csv").rename(columns={"thr": "threshold", "thr_ok": "threshold_ok",
                    "strict_recall": "R90_recall", "strict_precision": "R90_precision"})
        pf["feature_set"] = pf["feature_set"].replace(ALIASES)
        expected_folds = {(seed, fold) for seed in SEEDS for fold in range(K_OUTER)}
        for name, rows in pf.groupby("feature_set"):
            if rows.duplicated(["seed", "outer_fold"]).any() or set(map(tuple, rows[["seed", "outer_fold"]].to_numpy())) != expected_folds:
                raise ValueError(f"{tag}/{name}: expected all 10 seeds and 5 folds")
        names = ["L6_full"]
        if label in ["headline", "coding_nontruncating"]:
            names.append("base_freq_only")
            table = paired_ladder_table(pf)
            extra = pd.DataFrame([paired_delta(pf, "base_freq_only", "L6_full")])
            table = pd.concat([table, extra], ignore_index=True)
            table.insert(0, "cohort", label)
            descriptive.append(table)
            ladder_table(pf).to_csv(args.output / f"ladder_{label}.csv", index=False, lineterminator="\n")
        frames = {}
        for name in names:
            frames[name] = read(f"oof_{tag}_{name}.csv").rename(columns={
                "p_strict": "prob", "p4_gene_group": "gene_group", "p4_review_star": "review_star"})
        for name, frame in frames.items():
            if label != "random_split" and (frame.groupby(["seed", "gene_group"])["outer_fold"].nunique() != 1).any():
                raise ValueError(f"{tag}/{name}: gene group crosses test folds")
        if "base_freq_only" in frames:
            keys = ["seed", "VariationID", "outer_fold", "y", "gene_group", "consequence_class"]
            left = frames["L6_full"][keys].sort_values(keys[:2]).reset_index(drop=True)
            right = frames["base_freq_only"][keys].sort_values(keys[:2]).reset_index(drop=True)
            if not left.equals(right):
                raise ValueError(f"{tag}: models do not share identical evaluation rows")
        result = pooled_oof_table(pf, frames, n_boot=args.boot)
        summary["cohorts"][label] = result
        for metric in ["AUC", "R90_recall", "R90_precision"]:
            estimates.append({"cohort": label, "metric": metric,
                              **result["bootstrap_ci"][metric],
                              **{k: result["bootstrap_ci"]["metadata"][k] for k in ["n_variants", "n_groups", "n_boot"]}})
        if label == "headline":
            stratified_table(pf, frames).to_csv(args.output / "stratified_headline.csv", index=False, lineterminator="\n")
        if "base_freq_only" in frames:
            # Different estimand from mean per-fold deltas; no inferential claim.
            base = pooled_oof_table(pf, frames, feature_set="base_freq_only", n_boot=args.boot)
            comparisons.append({"cohort": label, "contrast": "L6_full - base_freq_only",
                                "estimand": "difference_of_seed_averaged_pooled_auc",
                                "full_auc": result["bootstrap_ci"]["AUC"]["estimate"],
                                "frequency_only_auc": base["bootstrap_ci"]["AUC"]["estimate"],
                                "delta": result["bootstrap_ci"]["AUC"]["estimate"] - base["bootstrap_ci"]["AUC"]["estimate"],
                                "inference": "descriptive_only"})
        print(f"{label}: n={result['n_variants']}, groups={result['bootstrap_ci']['metadata']['n_groups']}", flush=True)
    # Fail if the archive does not reproduce the committed reference cohort/calls.
    repo = Path(__file__).resolve().parents[1]
    reference = json.loads((repo / "results/summary.json").read_text(encoding="utf-8"))
    for label, key in [("headline", "pooled_headline"), ("coding_nontruncating", "pooled_missense")]:
        actual, expected = summary["cohorts"][label], reference[key]
        if actual["n_variants"] != expected["n_variants"] or actual["confusion"] != expected["confusion"]:
            raise ValueError(f"{label} does not match published reference counts and confusion")
    gap = reference["split_gap"]
    for label, key in [("headline", "gene_held_out_auc"), ("random_split", "random_split_auc")]:
        if abs(summary["cohorts"][label]["bootstrap_ci"]["AUC"]["estimate"] - gap[key]) > 0.000051:
            raise ValueError(f"{label} pooled AUC does not match published split comparison")
    summary["reference_checks"] = "cohort_counts_confusion_and_split_auc_match"
    summary["scope"] = "Post hoc conditional uncertainty; fixed models and scores, no external validation or calibration. Paired contrasts descriptive; no p-values."
    pd.DataFrame(estimates).to_csv(args.output / "cluster_intervals.csv", index=False, lineterminator="\n")
    pd.concat(descriptive, ignore_index=True).to_csv(args.output / "paired_fold_deltas.csv", index=False, lineterminator="\n")
    pd.DataFrame(comparisons).to_csv(args.output / "pooled_auc_comparisons.csv", index=False, lineterminator="\n")
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    sources = [repo / "src/hlpath/metrics.py", repo / "src/hlpath/analysis.py", repo / "src/hlpath/config.py", Path(__file__).resolve()]
    manifest = {"schema_version": "2.0", "inputs_sha256": inputs,
                "code_sha256": {str(p.relative_to(repo)).replace("\\", "/"): sha(p) for p in sources},
                "reference_summary_sha256": sha(repo / "results/summary.json"),
                "code_base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
                "code_is_local_revision": True, "bootstrap_seed": 0, "n_boot": args.boot,
                "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
                "outputs_sha256": {p.name: sha(p) for p in args.output.iterdir() if p.is_file()}}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
