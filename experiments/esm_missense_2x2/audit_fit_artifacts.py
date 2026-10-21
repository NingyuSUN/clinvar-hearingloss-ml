#!/usr/bin/env python3
"""Independent fit-artifact audit without new-runner functions or training.

The base commit's balanced_group_folds is the sole reused split algorithm.
The audit independently joins source files, reconstructs roles and medians,
selects the iteration/positive-score R90 quantile and replays saved models."""
import argparse
import hashlib
import importlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

BASE_COMMIT = "d27ff142fed63a110114f8efa928cdd26bedd2df"
INPUT_HASHES = {
    "matrix": "225fbf1c5060baadd157149b2c43424614a0fd88550ce22ecc6920e3ca17aceb",
    "plm_features": "16b463e18e61445213228f924f6188c30549897fd294cc7f76b6f38a7125667b",
}
SEEDS, ARMS = tuple(range(101, 111)), tuple("ABCD")
ROLES = ("inner_train", "inner_val", "outer_test")
CP_FILES = {"test_predictions.csv", "validation_predictions.csv", "fold.json",
            "validation_logloss.json", "model.ubj"}
SPLIT_COLUMNS = ["VariationID", "y", "gene_group", "plm_scored", "seed", "outer_fold", "arm", "role"]


def need(condition, text):
    if not condition:
        raise ValueError(text)


def digest(path):
    with Path(path).open("rb") as stream:
        h = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def read_csv(path):
    return pd.read_csv(path, dtype={"VariationID": str, "gene_group": str, "group": str},
                       float_precision="round_trip")


def counts(frame):
    return {"n": len(frame), "positive": int(frame.y.sum()),
            "negative": int((frame.y == 0).sum()), "gene_groups": int(frame.gene_group.nunique())}


def check_identity_columns(frame, context):
    for key in ("VariationID", "gene_group"):
        need(frame[key].notna().all(), f"{context}: missing {key}")
        value = frame[key].astype(str)
        need((value == value.str.strip()).all() and
             (~value.str.lower().isin(["", "nan", "none", "null"])).all(), f"{context}: invalid {key}")
    need(frame.y.isin([0, 1]).all(), f"{context}: nonbinary label")


def source_modules(repo, runner_dir, manifest):
    """Prove imported base files equal Git's fixed commit, before importing them."""
    module_names = ("config", "protocol", "pipeline", "plm", "features", "metrics")
    hashes = {}
    for name in module_names:
        relative = f"src/hlpath/{name}.py"
        frozen = subprocess.run(["git", "-C", str(repo), "show", f"{BASE_COMMIT}:{relative}"],
                                check=True, capture_output=True).stdout
        local = (repo / relative).read_bytes()
        need(local == frozen, f"base source differs from fixed commit: {relative}")
        hashes[f"hlpath.{name}"] = hashlib.sha256(local).hexdigest()
    for name, filename in (("fit_core", "fit_core.py"), ("run_missense_2x2.py", "run_missense_2x2.py")):
        hashes[name] = digest(runner_dir / filename)
    need(hashes == manifest["code_sha256"], "current code hashes differ from frozen manifest")
    sys.path.insert(0, str(repo / "src"))
    pipeline = importlib.import_module("hlpath.pipeline")
    protocol = importlib.import_module("hlpath.protocol")
    config = importlib.import_module("hlpath.config")
    need(Path(pipeline.__file__).resolve() == (repo / "src/hlpath/pipeline.py").resolve(),
         "wrong base pipeline import")
    return pipeline.balanced_group_folds, protocol, config, hashes


def load_independent_inputs(matrix_path, plm_path, features):
    # Match the frozen base matrix loader's standard pandas float parsing;
    # the PLM feature table uses round_trip in the declared new runner.
    matrix = pd.read_csv(matrix_path, low_memory=False)
    pf = pd.read_csv(plm_path, float_precision="round_trip")
    for frame, name in ((matrix, "matrix"), (pf, "PLM")):
        need(frame.VariationID.notna().all() and not frame.VariationID.duplicated().any(),
             f"{name}: missing or duplicate source ID")
        frame["VariationID"] = frame.VariationID.astype(str)
    matrix["y"] = pd.to_numeric(matrix.y, errors="coerce").astype("Int64")
    need(set(features) <= set(matrix), "missing base features")
    need(pf.plm_scored.isin([0, 1]).all(), "invalid source PLM scored flags")
    need(np.isfinite(pf.loc[pf.plm_scored == 1, "esm1v_llr"]).all(), "missing eligible ESM score")
    need(pf.loc[pf.plm_scored == 0, "esm1v_llr"].isna().all(), "unscored record has ESM1v value")
    joined = matrix.merge(pf[["VariationID", "plm_scored", "esm1v_llr"]], how="left",
                          on="VariationID", validate="one_to_one", sort=False)
    joined["plm_scored"] = joined.plm_scored.fillna(0).astype(int)
    selected = joined.loc[(joined.p4_headline_cohort == 1) &
                          (joined.consequence_class == "coding_nontruncating")].copy()
    selected["gene_group"] = selected.p4_gene_group
    selected["y"] = selected.y.astype(int)
    check_identity_columns(selected, "mixed cohort")
    selected["gene_group"] = selected.gene_group.astype(str)
    for feature in [*features, "esm1v_llr"]:
        selected[feature] = pd.to_numeric(selected[feature], errors="raise").astype(float)
        need(not np.isinf(selected[feature].to_numpy()).any(), f"infinite feature {feature}")
    need(counts(selected) == {"n": 3121, "positive": 989, "negative": 2132, "gene_groups": 148},
         "mixed source denominator mismatch")
    scored = selected[selected.plm_scored == 1]
    need(counts(scored) == {"n": 1328, "positive": 925, "negative": 403, "gene_groups": 135},
         "eligible missense source denominator mismatch")
    return selected.reset_index(drop=True)


def reconstruct_splits(data, assign):
    sizes = data.groupby("gene_group").size().to_dict()
    prevalence = data.groupby("gene_group").y.mean().to_dict()
    reconstructed, outer_record, inner_record, blocks = {}, {}, {}, []
    for seed in SEEDS:
        outer = assign(sizes, prevalence, 5, seed)
        outer_record[str(seed)] = outer
        for fold in range(5):
            candidates = [g for g in outer if outer[g] != fold]
            inner_raw = assign({g: sizes[g] for g in candidates},
                               {g: prevalence[g] for g in candidates}, 5, seed * 97 + fold)
            inner = {g: int(inner_raw[g] == 0) for g in candidates}
            inner_record[f"{seed}:{fold}"] = inner
            outside = data.gene_group.map(outer).eq(fold)
            validation = ~outside & data.gene_group.map(inner).eq(1)
            training = ~outside & data.gene_group.map(inner).eq(0)
            need((outside.astype(int) + validation.astype(int) + training.astype(int) == 1).all(),
                 "reconstructed roles do not partition mixed input")
            raw_roles = [data.loc[m] for m in (training, validation, outside)]
            gs = [set(frame.gene_group) for frame in raw_roles]
            need(not (gs[0] & gs[1] or gs[0] & gs[2] or gs[1] & gs[2]), "source groups leak across roles")
            for arm in ARMS:
                used = {}
                for role, frame in zip(ROLES, raw_roles):
                    chosen = frame if role == "inner_train" and arm in "AB" else frame[frame.plm_scored == 1]
                    need(set(chosen.y) == {0, 1}, f"one-class reconstructed role: {seed}/{fold}/{arm}/{role}")
                    used[role] = chosen.copy()
                    part = chosen[["VariationID", "y", "gene_group", "plm_scored"]].copy()
                    part = part.assign(seed=seed, outer_fold=fold, arm=arm, role=role)
                    blocks.append(part)
                reconstructed[seed, fold, arm] = used
    return reconstructed, pd.concat(blocks, ignore_index=True), outer_record, inner_record


def compare_split_tables(actual, expected):
    need(list(actual.columns) == SPLIT_COLUMNS, "split-manifest columns/order differ")
    need(not actual.duplicated(["seed", "outer_fold", "arm", "role", "VariationID"]).any(),
         "duplicate split membership")
    # Exact row order matters for reproducible training and aligns every role.
    pd.testing.assert_frame_equal(actual.reset_index(drop=True), expected[SPLIT_COLUMNS].reset_index(drop=True),
                                  check_dtype=False, check_exact=True)


def expected_medians(train, features, continuous):
    med = {feature: float(train[feature].median()) for feature in features if feature in continuous}
    need(all(np.isfinite(value) for value in med.values()), "all-missing/nonfinite training median")
    return med


def imputed_matrix(frame, features, medians):
    values = frame[features].copy()
    for feature in features:
        values[feature] = values[feature].fillna(medians.get(feature, 0.0))
    result = values.to_numpy(dtype=np.float64)
    need(np.isfinite(result).all(), "nonfinite value after independent imputation")
    return result


def independent_r90(y, scores):
    y, scores = np.asarray(y), np.asarray(scores, dtype=float)
    need(set(y) == {0, 1} and np.isfinite(scores).all(), "R90 needs finite scores and both classes")
    positives = np.sort(scores[y == 1])[::-1]
    # The largest feasible threshold minimizes FP and resolves ties upward.
    # Compute the needed rank by the original tp/n comparison, avoiding a
    # floating-point ceil(0.9*n) boundary difference at integer products.
    k = next(k for k in range(1, len(positives) + 1) if k / len(positives) >= 0.9)
    return float(positives[k - 1])


def selected_logloss(y, scores):
    p = np.clip(np.asarray(scores, dtype=float), 1e-16, 1 - 1e-16)
    y = np.asarray(y, dtype=int)
    return float(-np.mean(np.where(y == 1, np.log(p), np.log1p(-p))))


def validate_prediction_table(frame, reference, identity, threshold):
    expected = ["VariationID", "y", "gene_group", "score", "threshold", "decision", "seed", "outer_fold", "arm"]
    need(list(frame.columns) == expected, "checkpoint prediction schema mismatch")
    need(not frame.VariationID.duplicated().any(), "duplicate checkpoint prediction ID")
    for col in ("VariationID", "y", "gene_group"):
        need(list(frame[col]) == list(reference[col]), f"checkpoint source alignment mismatch: {col}")
    for col, value in identity.items():
        need(frame[col].eq(value).all(), f"checkpoint identity mismatch: {col}")
    need(np.isfinite(frame.score).all() and frame.score.between(0, 1).all(), "invalid checkpoint score")
    need(frame.threshold.eq(threshold).all(), "checkpoint threshold mismatch")
    need(frame.decision.isin([0, 1]).all() and
         frame.decision.eq(frame.score.ge(threshold).astype(int)).all(), "checkpoint decision mismatch")


def run_audit(args, result):
    manifest_path = args.run_dir / "frozen_manifest.json"
    manifest, manifest_hash = read_json(manifest_path), digest(manifest_path)
    need(manifest["source_commit"] == BASE_COMMIT, "wrong source commit in manifest")
    actual_inputs = {"matrix": digest(args.matrix), "plm_features": digest(args.plm_features)}
    need(actual_inputs == INPUT_HASHES == manifest["inputs_sha256"], "fixed input hash mismatch")
    assign, protocol, config, code_hashes = source_modules(args.repo_root, args.runner_dir, manifest)
    base_features = list(protocol.ALL_FEATURES)
    need(len(base_features) == 18 and len(set(base_features)) == 18, "base features not original18")
    arm_settings = {a: {"training": "mixed" if a in "AB" else "scored_missense",
                        "features": base_features + (["esm1v_llr"] if a in "BD" else [])} for a in ARMS}
    need(manifest["arms"] == arm_settings, "arm feature/training specification mismatch")
    cfg = manifest["config"]
    need(1 <= cfg["threads"] <= 10, "training threads outside cap")
    params101 = config.xgb_params(101)
    params101["nthread"] = cfg["threads"]
    need(cfg == {"seeds": list(SEEDS), "k_outer": 5, "num_boost_round": 700,
                 "posthoc_min_iter": 10, "inner_valid_fraction": 0.2, "target_recall": 0.9,
                 "threads": cfg["threads"], "xgb_params_seed101": params101,
                 "validation_population": "scored_missense", "refit": False}, "frozen scientific config differs")
    versions = {"python": platform.python_version(), **{name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "xgboost", "scikit-learn", "scipy")}}
    need(versions == manifest["versions"], "audit runtime differs from frozen fit runtime")
    data = load_independent_inputs(args.matrix, args.plm_features, base_features)
    need(manifest["cohort_counts"] == {"mixed": counts(data), "missense": counts(data[data.plm_scored == 1])},
         "manifest cohort counts differ from source")
    roles, split, outer, inner = reconstruct_splits(data, assign)
    split_path = args.run_dir / "split_manifest.csv"
    need(digest(split_path) == manifest["split_manifest_sha256"], "split file hash mismatch")
    compare_split_tables(read_csv(split_path), split)
    need(outer == manifest["outer_assignments"] and inner == manifest["inner_assignments"],
         "manifest assignments differ from fixed-base reconstruction")
    continuous = set(protocol.CONTINUOUS) | {"esm1v_llr"}
    prepared = []
    for (seed, fold, arm), used in roles.items():
        fs = arm_settings[arm]["features"]
        medians = expected_medians(used["inner_train"], fs, continuous)
        # Exercise all role transforms before any fit, without retaining raw arrays.
        for role in ROLES:
            imputed_matrix(used[role], fs, medians)
        prepared.append({"seed": seed, "outer_fold": fold, "arm": arm,
                         "roles_counts": {r: counts(used[r]) for r in ROLES},
                         "expected_imputation_medians": medians})
    result.update(source_commit=BASE_COMMIT, inputs_sha256=actual_inputs, code_sha256=code_hashes,
                  versions=versions, frozen_manifest_sha256=manifest_hash, split_manifest_sha256=digest(split_path),
                  cohort_counts=manifest["cohort_counts"], split_rows=len(split), assignments=50,
                  planned_fits=200, prepared_roles_and_medians=prepared,
                  preparation_checks_passed=True,
                  preparation_scope="independent source join, exact split row order/identity, group isolation, per-arm medians")
    if args.stage == "splits":
        result.update(passed=True, scope="PREPARE_ONLY; no fitted artifacts checked", replayed_models=0)
        return
    checkpoint_root = args.run_dir / "checkpoints"
    expected_names = {f"s{s}_f{f}_{a}" for s in SEEDS for f in range(5) for a in ARMS}
    # Check only completed checkpoint-shaped directories; unrelated hidden temp
    # files/locks are neither changed nor interpreted as permission to remove.
    actual_names = {p.name for p in checkpoint_root.iterdir() if p.is_dir() and p.name.startswith("s")}
    need(actual_names == expected_names, "expected exactly 200 named checkpoint directories")
    completed_hashes, all_test, fit_audits = {}, [], []
    max_prediction_error, max_logloss_error = 0.0, 0.0
    for item in prepared:
        seed, fold, arm = item["seed"], item["outer_fold"], item["arm"]
        name = f"s{seed}_f{fold}_{arm}"
        cp = checkpoint_root / name
        identity = {"seed": seed, "outer_fold": fold, "arm": arm}
        completion = read_json(cp / "COMPLETE.json")
        need(completion["identity"] == identity and completion["frozen_manifest_sha256"] == manifest_hash,
             f"{name}: checkpoint identity/binding mismatch")
        need(set(completion["files_sha256"]) == CP_FILES, f"{name}: incomplete checkpoint file manifest")
        for filename, expected_hash in completion["files_sha256"].items():
            need(digest(cp / filename) == expected_hash, f"{name}: changed checkpoint file {filename}")
        completed_hashes[name] = digest(cp / "COMPLETE.json")
        row = read_json(cp / "fold.json")
        for key, value in identity.items():
            need(row[key] == value, f"{name}: fold metadata identity mismatch")
        need(row["roles_counts"] == item["roles_counts"], f"{name}: role counts mismatch")
        fs = arm_settings[arm]["features"]
        need(row["features"] == fs, f"{name}: fitted feature whitelist/order differs")
        med = item["expected_imputation_medians"]
        need(set(row["imputation_medians"]) == set(med), f"{name}: imputation feature keys differ")
        need(all(np.isclose(row["imputation_medians"][f], v, rtol=0, atol=1e-12) for f, v in med.items()),
             f"{name}: imputation medians differ from source inner train")
        params = config.xgb_params(seed)
        params["nthread"] = cfg["threads"]
        need(row["xgb_params"] == params, f"{name}: declared XGBoost params differ")
        curve = np.asarray(read_json(cp / "validation_logloss.json"), dtype=float)
        need(curve.shape == (700,) and np.isfinite(curve).all() and (curve >= 0).all(),
             f"{name}: invalid 700-round validation logloss curve")
        best = int(np.argmin(curve[10:])) + 10
        need(row["best_iteration"] == best, f"{name}: selected iteration not first argmin after index10")
        threshold = float(row["threshold"])
        need(np.isfinite(threshold) and 0 <= threshold <= 1 and row["threshold_reason"] == "ok",
             f"{name}: invalid R90 metadata")
        used = roles[seed, fold, arm]
        validation, test = read_csv(cp / "validation_predictions.csv"), read_csv(cp / "test_predictions.csv")
        validate_prediction_table(validation, used["inner_val"], identity, threshold)
        validate_prediction_table(test, used["outer_test"], identity, threshold)
        need(independent_r90(validation.y.to_numpy(), validation.score.to_numpy()) == threshold,
             f"{name}: threshold differs from independent positive-order-statistic R90")
        loss_error = abs(selected_logloss(validation.y, validation.score) - curve[best])
        need(loss_error <= 1e-6, f"{name}: selected validation logloss differs from saved prediction scores")
        max_logloss_error = max(max_logloss_error, loss_error)
        booster = xgb.Booster()
        booster.load_model(str(cp / "model.ubj"))
        need(booster.feature_names == fs and booster.num_features() == len(fs), f"{name}: model feature schema differs")
        need(booster.num_boosted_rounds() == 700, f"{name}: saved model is not complete 700-round fit")
        booster.set_param({"nthread": args.replay_threads})
        errors = {}
        for role, stored in (("inner_val", validation), ("outer_test", test)):
            dm = xgb.DMatrix(imputed_matrix(used[role], fs, med), feature_names=fs)
            predicted = booster.predict(dm, iteration_range=(0, best + 1))
            error = float(np.max(np.abs(predicted.astype(float) - stored.score.to_numpy())))
            need(error <= 1e-7, f"{name}: saved model prediction replay mismatch on {role}")
            need(np.array_equal((predicted >= threshold).astype(int), stored.decision.to_numpy()),
                 f"{name}: replay changes a threshold decision on {role}")
            errors[role] = error
            max_prediction_error = max(max_prediction_error, error)
        all_test.append(test.rename(columns={"gene_group": "group", "outer_fold": "fold"}))
        fit_audits.append({**identity, "best_iteration": best, "threshold": threshold,
                           "validation_logloss_abs_error": loss_error, "replay_max_abs_error": errors})
        if len(fit_audits) % 20 == 0:
            print(json.dumps({"stage": "audit_replay", "checked_models": len(fit_audits), "total": 200}), flush=True)
    oof = read_csv(args.run_dir / "oof_predictions.csv")
    expected_oof = pd.concat(all_test, ignore_index=True)
    need(len(oof) == 53120 and not oof.duplicated(["arm", "VariationID", "seed"]).any(), "final OOF incomplete/duplicate")
    pd.testing.assert_frame_equal(oof, expected_oof, check_dtype=False, check_exact=True)
    seal = read_json(args.run_dir / "prediction_seal.json")
    need(seal["total_fits"] == 200 and seal["oof_rows"] == 53120 and
         seal["frozen_manifest_sha256"] == manifest_hash and
         seal["checkpoint_complete_sha256"] == completed_hashes, "final seal binding mismatch")
    need(set(seal["files_sha256"]) == {"oof_predictions.csv", "per_fold.csv"}, "final seal files differ")
    for filename, expected_hash in seal["files_sha256"].items():
        need(digest(args.run_dir / filename) == expected_hash, f"final seal hash differs: {filename}")
    result.update(passed=True, scope="COMPLETE_FIT_ARTIFACT_AUDIT", replayed_models=200,
                  validation_and_test_prediction_replays=400, final_oof_rows=len(oof),
                  prediction_seal_sha256=digest(args.run_dir / "prediction_seal.json"),
                  max_prediction_abs_error=max_prediction_error, max_selected_logloss_abs_error=max_logloss_error,
                  fit_audits=fit_audits,
                  limitations=["Checks engineering consistency and declared training isolation; no external/clinical validation",
                               "No training rerun; uses saved model predictions, metadata and audited source code",
                               "All 700 stored logloss values/argmin checked; only selected-round logloss independently replayed",
                               "Narrow inner-validation gene diversity and fixed-prediction uncertainty remain scientific limits"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--runner-dir", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--plm-features", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--stage", choices=("splits", "complete"), default="complete")
    parser.add_argument("--replay-threads", type=int, default=10)
    args = parser.parse_args()
    need(not args.out.exists(), "preserve prior audit receipt; choose a new output filename")
    need(1 <= args.replay_threads <= 10, "replay thread cap is 10")
    result = {"checked_utc": datetime.now(timezone.utc).isoformat(), "stage": args.stage,
              "script_sha256": digest(__file__), "passed": False}
    started = time.monotonic()
    try:
        run_audit(args, result)
    except Exception as exc:
        result["error"] = type(exc).__name__ + ": " + str(exc)
    result["elapsed_seconds"] = time.monotonic() - started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: result.get(k) for k in ("passed", "stage", "scope", "replayed_models", "error")},
                     ensure_ascii=False), flush=True)
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
