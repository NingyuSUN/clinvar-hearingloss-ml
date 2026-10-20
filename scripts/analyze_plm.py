#!/usr/bin/env python3
"""Prespecified analysis for PLM-ESM-2026-10 (see docs/plm_protocol.md).

    python scripts/analyze_plm.py --runs runs/plm_esm_2026-10 --out results/plm_esm_2026-10

Verifies the prediction seal, then reports per-model pooled metrics with
gene-group intervals, the prespecified paired contrasts, zero-shot ESM AUCs
on scored missense variants, matched-fold descriptive deltas and the
ORIGINAL18 reproduction check against the published statistics revision.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from hlpath.analysis import _pooled
from hlpath.compare import paired_group_bootstrap_auc
from hlpath.config import REPO_ROOT
from hlpath.metrics import bootstrap_ci, paired_delta
from hlpath.plm import ESM1V_MODELS, ESM2_MODEL, PLM_FEATURE_SETS, PROTOCOL_ID, SCORED_FLAG

N_BOOT = 2000
BOOT_SEED = 0
# (id, tier, task, evaluation set, model a, model b): delta = AUC(b) - AUC(a).
CONTRASTS = [
    ("P1", "primary", "plm_codingnt", "all", "ORIGINAL18", "ORIGINAL18_ESM1V"),
    ("S1", "secondary", "plm_headline", "all", "ORIGINAL18", "ORIGINAL18_ESM1V"),
    ("S2", "secondary", "plm_codingnt", "scored_missense", "ORIGINAL18", "ORIGINAL18_ESM1V"),
    ("S3", "secondary", "plm_codingnt", "all", "ORIGINAL18_SCOREDFLAG", "ORIGINAL18_ESM1V"),
    ("D1", "descriptive", "plm_codingnt", "all", "ORIGINAL18", "ORIGINAL18_ESM2"),
    ("D2", "descriptive", "plm_headline", "all", "ORIGINAL18", "ORIGINAL18_ESM2"),
    ("D3", "descriptive", "plm_codingnt", "scored_missense", "ORIGINAL18", "ORIGINAL18_ESM2"),
    ("D4", "descriptive", "plm_codingnt", "all", "ORIGINAL18", "ORIGINAL18_SCOREDFLAG"),
    ("D5", "descriptive", "plm_codingnt", "all", "ORIGINAL18", "ESM1V_ONLY"),
    ("D6", "descriptive", "plm_codingnt", "scored_missense", "ORIGINAL18", "ESM1V_ONLY"),
    ("D7", "descriptive", "plm_headline", "all", "ORIGINAL18_SCOREDFLAG", "ORIGINAL18_ESM1V"),
]
# Published empirical estimates for the identical ORIGINAL18 (L6_full) protocol.
REFERENCE = {"plm_codingnt": "coding_nontruncating", "plm_headline": "headline"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_seal(runs: Path) -> dict:
    seal = json.loads((runs / "prediction_seal.json").read_text())
    for name, digest in seal["predictions_sha256"].items():
        if sha256(runs / name) != digest:
            raise SystemExit(f"{name} changed after sealing")
    if seal["feature_sets"] != PLM_FEATURE_SETS:
        raise SystemExit("sealed feature sets differ from the protocol")
    return seal


def pooled_sets(runs: Path, tag: str) -> dict[str, pd.DataFrame]:
    per_fold = pd.read_csv(runs / f"per_fold_{tag}.csv")
    out = {}
    for name in PLM_FEATURE_SETS:
        oof = pd.read_csv(runs / f"oof_{tag}__{name}.csv")
        pf = per_fold[per_fold["feature_set"] == name]
        out[name] = _pooled(pf, oof).sort_values("VariationID").reset_index(drop=True)
    ids = {tuple(df["VariationID"]) for df in out.values()}
    if len(ids) != 1:
        raise SystemExit(f"{tag}: models cover different variants")
    return out | {"_per_fold": per_fold}


def evaluation_rows(pooled: pd.DataFrame, subset: str, scored_ids: set) -> pd.DataFrame:
    if subset == "all":
        return pooled
    if subset == "scored_missense":
        return pooled[pooled["VariationID"].isin(scored_ids)]
    raise ValueError(subset)


def model_metrics(task: str, subset: str, name: str, df: pd.DataFrame) -> dict:
    pred = (df["hit"] >= 0.5).to_numpy().astype(int)
    ci = bootstrap_ci(df["y"].to_numpy().astype(int), df["prob"].to_numpy(), pred,
                      n_boot=N_BOOT, seed=BOOT_SEED, groups=df["gene_group"].to_numpy())
    row = {"task": task, "evaluation_set": subset, "model": name,
           "n_variants": len(df), "n_positive": int(df["y"].sum()),
           "n_gene_groups": ci["metadata"]["n_groups"]}
    for key in ["AUC", "R90_recall", "R90_precision"]:
        for part in ["estimate", "lo", "hi", "n_valid"]:
            row[f"{key}_{part}"] = ci[key][part]
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "results/plm_esm_2026-10")
    ap.add_argument("--features", type=Path, default=REPO_ROOT / "data/plm/plm_features.csv")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    seal = verify_seal(args.runs)
    if seal["inputs_sha256"]["data/plm/plm_features.csv"] != sha256(args.features):
        raise SystemExit("features differ from those used for fitting")
    features = pd.read_csv(args.features)
    scored_ids = set(features.loc[features[SCORED_FLAG] == 1, "VariationID"])

    pooled = {tag: pooled_sets(args.runs, tag) for tag in seal["tasks"]}

    metrics = []
    for tag in pooled:
        subsets = ["all", "scored_missense"] if tag == "plm_codingnt" else ["all"]
        for subset in subsets:
            for name in PLM_FEATURE_SETS:
                df = evaluation_rows(pooled[tag][name], subset, scored_ids)
                metrics.append(model_metrics(tag, subset, name, df))
    metrics = pd.DataFrame(metrics)

    contrasts = []
    for cid, tier, tag, subset, a, b in CONTRASTS:
        da = evaluation_rows(pooled[tag][a], subset, scored_ids)
        db = evaluation_rows(pooled[tag][b], subset, scored_ids)
        if not (da["VariationID"].to_numpy() == db["VariationID"].to_numpy()).all():
            raise SystemExit(f"{cid}: unaligned variants")
        res = paired_group_bootstrap_auc(da["y"].to_numpy(), da["prob"].to_numpy(), db["prob"].to_numpy(),
                                         groups=da["gene_group"].to_numpy(), n_boot=N_BOOT, seed=BOOT_SEED)
        d = res["delta_b_minus_a"]
        contrasts.append({
            "contrast": cid, "tier": tier, "task": tag, "evaluation_set": subset,
            "model_a": a, "model_b": b, "n_variants": res["metadata"]["n_variants"],
            "n_positive": res["metadata"]["n_positive"], "n_gene_groups": res["metadata"]["n_groups"],
            "auc_a": res["auc_a"]["estimate"], "auc_b": res["auc_b"]["estimate"],
            "delta": d["estimate"], "delta_lo": d["lo"], "delta_hi": d["hi"],
            "delta_bootstrap_mean": d["bootstrap_mean"], "fraction_draws_le_zero": d["fraction_draws_le_zero"],
            "n_valid_draws": d["n_valid"], "status": d["status"],
            "lower_bound_above_zero": res["superior"],
        })
    contrasts = pd.DataFrame(contrasts)

    # Zero-shot: -LLR ranks variants (lower LLR = less likely than wild type).
    base = evaluation_rows(pooled["plm_codingnt"]["ORIGINAL18"], "scored_missense", scored_ids)
    zs = base.merge(features, on="VariationID", how="left", validate="one_to_one")
    zero_shot = []
    llr_cols = {"esm1v_ensemble": "esm1v_llr", "esm2_650M": "esm2_llr",
                **{m.split("/")[-1]: f"llr__{m.split('/')[-1]}" for m in ESM1V_MODELS}}
    for label, col in llr_cols.items():
        if zs[col].isna().any():
            raise SystemExit(f"zero-shot {col} missing for scored variants")
        ci = bootstrap_ci(zs["y"].to_numpy().astype(int), -zs[col].to_numpy(), None,
                          n_boot=N_BOOT, seed=BOOT_SEED, groups=zs["gene_group"].to_numpy())
        zero_shot.append({"score": label, "column": col, "n_variants": len(zs),
                          "n_positive": int(zs["y"].sum()), "n_gene_groups": ci["metadata"]["n_groups"],
                          "AUC_estimate": ci["AUC"]["estimate"], "AUC_lo": ci["AUC"]["lo"],
                          "AUC_hi": ci["AUC"]["hi"]})
    zs_vs = paired_group_bootstrap_auc(zs["y"].to_numpy(), zs["prob"].to_numpy(), -zs["esm1v_llr"].to_numpy(),
                                       groups=zs["gene_group"].to_numpy(), n_boot=N_BOOT, seed=BOOT_SEED)
    zero_shot = pd.DataFrame(zero_shot)

    fold_deltas = []
    for cid, tier, tag, subset, a, b in CONTRASTS:
        if subset == "all":
            fold_deltas.append({"contrast": cid, "tier": tier, "task": tag,
                                **paired_delta(pooled[tag]["_per_fold"], a, b)})
    fold_deltas = pd.DataFrame(fold_deltas)

    published = json.loads((REPO_ROOT / "results/statistics_20260921/summary.json").read_text())
    reproduction = {}
    for tag, key in REFERENCE.items():
        if tag not in pooled:
            continue
        df = pooled[tag]["ORIGINAL18"]
        pred = (df["hit"] >= 0.5).astype(int)
        y = df["y"].astype(int)
        confusion = {"TP": int(((pred == 1) & (y == 1)).sum()), "FP": int(((pred == 1) & (y == 0)).sum()),
                     "TN": int(((pred == 0) & (y == 0)).sum()), "FN": int(((pred == 0) & (y == 1)).sum())}
        auc = metrics.query("task == @tag and evaluation_set == 'all' and model == 'ORIGINAL18'")["AUC_estimate"].item()
        ref = published["cohorts"][key]
        ref_auc = ref["bootstrap_ci"]["AUC"]["estimate"]
        reproduction[tag] = {
            "published_cohort": key, "n_variants": int(len(df)),
            "published_n_variants": ref["n_variants"],
            "auc_rerun": auc, "auc_published": ref_auc, "abs_difference": abs(auc - ref_auc),
            "confusion_rerun": confusion, "confusion_published": ref["confusion"],
            "exact_match": bool(abs(auc - ref_auc) < 1e-9 and confusion == ref["confusion"]),
        }

    paths = {
        "model_metrics.csv": metrics, "contrasts.csv": contrasts,
        "zero_shot_scored_missense.csv": zero_shot, "matched_fold_deltas.csv": fold_deltas,
    }
    for name, frame in paths.items():
        frame.to_csv(args.out / name, index=False, lineterminator="\n", float_format="%.6f")
    primary = contrasts[contrasts["contrast"] == "P1"].iloc[0]
    s3 = contrasts[contrasts["contrast"] == "S3"].iloc[0]
    if not primary["lower_bound_above_zero"]:
        verdict = "improvement_not_demonstrated"
    elif s3["lower_bound_above_zero"]:
        verdict = "improvement_demonstrated_beyond_scored_status_marker"
    else:
        verdict = "improvement_demonstrated_not_attributable_beyond_scored_status_marker"
    records = json.loads(contrasts.to_json(orient="records", double_precision=15))
    summary = {
        "protocol": PROTOCOL_ID, "training_performed_here": False,
        "n_boot": N_BOOT, "bootstrap_seed": BOOT_SEED,
        "primary": records[0], "primary_verdict": verdict,
        "contrasts": records,
        "zero_shot_vs_original18_scored_missense": {
            "auc_original18": zs_vs["auc_a"], "auc_zero_shot_esm1v": zs_vs["auc_b"],
            "delta_zero_shot_minus_original18": zs_vs["delta_b_minus_a"]},
        "reproduction_of_published_original18": reproduction,
        "n_scored_variants_headline": int(len(scored_ids)),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    code = ["scripts/analyze_plm.py", "scripts/run_plm_evaluation.py", "scripts/plm_build_features.py",
            "src/hlpath/plm.py", "src/hlpath/compare.py", "src/hlpath/pipeline.py",
            "src/hlpath/analysis.py", "src/hlpath/metrics.py", "kaggle/esm_score_kernel.py"]
    manifest = {
        "protocol": PROTOCOL_ID, "prediction_seal": seal,
        "prediction_seal_sha256": sha256(args.runs / "prediction_seal.json"),
        "features_sha256": sha256(args.features),
        "code_sha256": {c: sha256(REPO_ROOT / c) for c in code},
        "outputs_sha256": {n: sha256(args.out / n) for n in [*paths, "summary.json"]},
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"primary_verdict": verdict, "reproduction": {
        k: v["exact_match"] for k, v in reproduction.items()}}, indent=2))


if __name__ == "__main__":
    main()
