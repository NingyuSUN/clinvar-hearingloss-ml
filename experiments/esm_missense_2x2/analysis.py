#!/usr/bin/env python3
"""Paired, group-bootstrap analysis of four arms' fixed repeated OOF predictions.

Per-variant inputs and *_LOCAL_ONLY.csv outputs are not public aggregate results."""
import argparse
import csv
import hashlib
import json
import platform
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ARMS = ("A", "B", "C", "D")
SEEDS = tuple(range(101, 111))
CONTRASTS = (("D", "C", "primary"), ("B", "A", "exploratory"),
             ("C", "A", "exploratory"), ("D", "B", "exploratory"))
METRICS = ("auc", "recall", "fpr", "specificity")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_rows(rows, expected_n=1328, expected_positive=925, expected_groups=135,
                  seeds=SEEDS, folds=tuple(range(5))):
    require(0 < expected_positive < expected_n and expected_groups >= 2,
            "Need both expected classes and at least two groups")
    required = {"VariationID", "y", "group", "score", "threshold", "decision", "seed", "fold", "arm"}
    data, meta, assignments, thresholds = {}, {}, {}, {}
    group_folds = {}
    for line, row in enumerate(rows, 2):
        require(required <= row.keys() and all(row[k] is not None for k in required),
                f"Line {line}: missing required columns or values")
        vid, group, arm = str(row["VariationID"]), str(row["group"]), str(row["arm"])
        require(all(value == value.strip() for value in (vid, group, arm)),
                f"Line {line}: whitespace in ID/group/arm; repair source explicitly")
        require(vid.strip().lower() not in ("", "nan", "none", "null") and
                group.strip().lower() not in ("", "nan", "none", "null"), f"Line {line}: empty ID/group")
        require(arm in ARMS, f"Line {line}: unexpected arm {arm!r}")
        for field in ("y", "decision", "seed", "fold"):
            try:
                parsed = int(row[field])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Line {line}: {field} must be an integer") from exc
            require(str(row[field]).strip() == str(parsed), f"Line {line}: {field} must be an integer")
        y, decision, seed, fold = (int(row[k]) for k in ("y", "decision", "seed", "fold"))
        score, threshold = float(row["score"]), float(row["threshold"])
        require(y in (0, 1) and decision in (0, 1), f"Line {line}: non-binary label/decision")
        require(seed in seeds and fold in folds, f"Line {line}: unexpected seed/fold")
        require(np.isfinite(score) and 0 <= score <= 1, f"Line {line}: invalid probability")
        require(np.isfinite(threshold) and 0 <= threshold <= 1, f"Line {line}: invalid threshold")
        require(decision == int(score >= threshold), f"Line {line}: decision != score >= threshold")
        key = (arm, vid, seed)
        require(key not in data, f"Line {line}: duplicate arm/ID/seed: {key}")
        require(vid not in meta or meta[vid] == (y, group), f"Line {line}: label/group inconsistency")
        require((vid, seed) not in assignments or assignments[vid, seed] == fold,
                f"Line {line}: outer fold differs across arms")
        require((group, seed) not in group_folds or group_folds[group, seed] == fold,
                f"Line {line}: group split across outer test folds")
        require((arm, seed, fold) not in thresholds or thresholds[arm, seed, fold] == threshold,
                f"Line {line}: threshold differs within one arm/seed/fold")
        meta[vid] = y, group
        assignments[vid, seed] = fold
        group_folds[group, seed] = fold
        thresholds[arm, seed, fold] = threshold
        data[key] = score, decision, fold, threshold
    ids = sorted(meta)
    require(len(ids) == expected_n, f"Expected {expected_n} unique IDs, got {len(ids)}")
    require(sum(meta[v][0] for v in ids) == expected_positive, "Positive count mismatch")
    require(len({meta[v][1] for v in ids}) == expected_groups, "Group count mismatch")
    require(len(data) == len(ARMS) * len(ids) * len(seeds), "Incomplete arm/ID/seed panel")
    for arm in ARMS:
        for seed in seeds:
            require({data[arm, v, seed][2] for v in ids if (arm, v, seed) in data} == set(folds),
                    f"Missing outer fold in {arm}/{seed}")
            for vid in ids:
                require((arm, vid, seed) in data, f"Missing prediction {arm}/{vid}/{seed}")
    y = np.array([meta[v][0] for v in ids], dtype=int)
    groups = np.array([meta[v][1] for v in ids], dtype=str)
    scores, decisions = {}, {}
    for arm in ARMS:
        scores[arm] = np.array([[data[arm, v, seed][0] for seed in seeds] for v in ids])
        decisions[arm] = np.array([[data[arm, v, seed][1] for seed in seeds] for v in ids])
    return {"ids": ids, "y": y, "groups": groups, "scores": scores, "decisions": decisions,
            "seeds": tuple(seeds), "rows": len(data), "folds": tuple(folds)}


def weighted_auc(y, score, weight=None):
    """Mann-Whitney AUC with half credit for ties; weights replicate group draws."""
    y, score = np.asarray(y), np.asarray(score)
    weight = np.ones(len(y), dtype=float) if weight is None else np.asarray(weight, dtype=float)
    require(y.ndim == score.ndim == weight.ndim == 1 and len(y) == len(score) == len(weight),
            "AUC array mismatch")
    require(np.all(np.isfinite(weight)) and np.all(weight >= 0), "Invalid AUC weights")
    _, bins = np.unique(score, return_inverse=True)  # ascending score, ties share a bin
    pos = np.bincount(bins, weights=weight * y)
    neg = np.bincount(bins, weights=weight * (1 - y))
    p, n = float(pos.sum()), float(neg.sum())
    if p == 0 or n == 0:
        return None
    lower_neg = np.cumsum(neg) - neg
    return float(np.sum(pos * (lower_neg + 0.5 * neg)) / (p * n))


def metrics(y, score, decision, weight=None):
    y, decision = np.asarray(y), np.asarray(decision)
    w = np.ones(len(y), dtype=int) if weight is None else np.asarray(weight)
    tp = int(w[(y == 1) & (decision == 1)].sum())
    fn = int(w[(y == 1) & (decision == 0)].sum())
    fp = int(w[(y == 0) & (decision == 1)].sum())
    tn = int(w[(y == 0) & (decision == 0)].sum())
    return {"n": tp + fn + fp + tn, "positive": tp + fn, "negative": fp + tn,
            "TP": tp, "FN": fn, "FP": fp, "TN": tn,
            "recall": tp / (tp + fn) if tp + fn else None,
            "precision": tp / (tp + fp) if tp + fp else None,
            "fpr": fp / (fp + tn) if fp + tn else None,
            "specificity": tn / (fp + tn) if fp + tn else None,
            "auc": weighted_auc(y, score, w)}


def aggregate(panel):
    scores = {a: panel["scores"][a].mean(axis=1) for a in ARMS}
    # Ten votes with five or more positive: a 5/5 tie is positive by frozen rule.
    decisions = {a: (panel["decisions"][a].sum(axis=1) >= len(panel["seeds"]) / 2).astype(int)
                 for a in ARMS}
    summary = {a: metrics(panel["y"], scores[a], decisions[a]) for a in ARMS}
    return scores, decisions, summary


def bootstrap(panel, scores, decisions, summary, iterations=2000, seed=0):
    require(iterations >= 2, "Need at least two bootstrap draws")
    unique_groups, inverse = np.unique(panel["groups"], return_inverse=True)
    rng = np.random.default_rng(seed)
    draws = {f"{a}-{b}": {m: [] for m in METRICS} for a, b, _ in CONTRASTS}
    draw_hash = hashlib.sha256()
    invalid = 0
    for _ in range(iterations):
        selected = rng.integers(0, len(unique_groups), size=len(unique_groups))
        multiplicity = np.bincount(selected, minlength=len(unique_groups)).astype("<i8")
        draw_hash.update(multiplicity.tobytes())
        weight = multiplicity[inverse]
        # One common group draw feeds all arms and all contrasts.
        resampled = {a: metrics(panel["y"], scores[a], decisions[a], weight) for a in ARMS}
        if resampled["A"]["auc"] is None:
            invalid += 1
            continue
        for a, b, _ in CONTRASTS:
            for metric in METRICS:
                draws[f"{a}-{b}"][metric].append(resampled[a][metric] - resampled[b][metric])
    valid = iterations - invalid
    intervals_available = valid >= max(100, int(np.ceil(0.9 * iterations)))
    contrasts = []
    for a, b, role in CONTRASTS:
        for metric in METRICS:
            low, high = (np.quantile(draws[f"{a}-{b}"][metric], [0.025, 0.975], method="linear")
                         if intervals_available else (None, None))
            contrasts.append({"contrast": f"{a}-{b}", "role": role, "metric": metric,
                              "difference": summary[a][metric] - summary[b][metric],
                              "ci95_low": float(low) if low is not None else None,
                              "ci95_high": float(high) if high is not None else None,
                              "interval_status": "available" if intervals_available else "withheld_insufficient_valid_draws",
                              "valid_draws": valid})
    return contrasts, {"iterations_requested": iterations, "seed": seed,
                       "groups": len(unique_groups), "invalid_single_class_draws": invalid,
                       "valid_draws": valid, "intervals_available": intervals_available,
                       "group_multiplicities_sha256": draw_hash.hexdigest(),
                       "method": "paired gene-group percentile bootstrap; common draws across arms",
                       "scope": "conditional on fixed OOF predictions; no refitting uncertainty",
                       "multiplicity": "no multiplicity adjustment; AUC and decision effects reported together"}


def discordance(panel, scores, decisions):
    counts, rows = {}, []
    for label, title in ((1, "pathogenic"), (0, "benign")):
        ix = panel["y"] == label
        counts[title] = {"n": int(ix.sum()), "C0_D0": int((ix & (decisions["C"] == 0) & (decisions["D"] == 0)).sum()),
                         "C0_D1": int((ix & (decisions["C"] == 0) & (decisions["D"] == 1)).sum()),
                         "C1_D0": int((ix & (decisions["C"] == 1) & (decisions["D"] == 0)).sum()),
                         "C1_D1": int((ix & (decisions["C"] == 1) & (decisions["D"] == 1)).sum())}
    for i, vid in enumerate(panel["ids"]):
        if decisions["C"][i] != decisions["D"][i]:
            rows.append({"VariationID": vid, "y": int(panel["y"][i]), "group": str(panel["groups"][i]),
                         "C_score": float(scores["C"][i]), "D_score": float(scores["D"][i]),
                         "C_decision": int(decisions["C"][i]), "D_decision": int(decisions["D"][i]),
                         "change": "corrected" if decisions["D"][i] == panel["y"][i] else "new_error"})
    return counts, rows


def write_csv(path, rows, fields=None):
    fields = fields if fields is not None else (list(rows[0]) if rows else None)
    require(bool(fields), "Cannot write an empty CSV without fieldnames")
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--expected-n", type=int, default=1328)
    parser.add_argument("--expected-positive", type=int, default=925)
    parser.add_argument("--expected-groups", type=int, default=135)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=0)
    args = parser.parse_args()
    require(not args.out.exists() or not any(args.out.iterdir()), "Output directory must be new or empty")
    with args.input.open(encoding="utf-8-sig", newline="") as stream:
        panel = validate_rows(csv.DictReader(stream), args.expected_n, args.expected_positive, args.expected_groups)
    scores, decisions, summary = aggregate(panel)
    contrasts, bootmeta = bootstrap(panel, scores, decisions, summary, args.bootstrap, args.bootstrap_seed)
    changed_counts, changed_rows = discordance(panel, scores, decisions)
    seed_rows = [{"arm": a, "seed": seed, **metrics(panel["y"], panel["scores"][a][:, j],
                                                    panel["decisions"][a][:, j])}
                 for a in ARMS for j, seed in enumerate(panel["seeds"])]
    seed_differences = [{"contrast": f"{a}-{b}", "role": role, "seed": seed,
                         **{m: seed_rows[ARMS.index(a) * len(panel["seeds"]) + j][m] -
                               seed_rows[ARMS.index(b) * len(panel["seeds"]) + j][m] for m in METRICS}}
                        for a, b, role in CONTRASTS for j, seed in enumerate(panel["seeds"])]
    aggregate_rows = [{"VariationID": vid, "y": int(panel["y"][i]), "group": str(panel["groups"][i]),
                       **{f"{a}_{m}": (float(scores[a][i]) if m == "score" else int(decisions[a][i]))
                          for a in ARMS for m in ("score", "decision")}}
                      for i, vid in enumerate(panel["ids"])]
    result = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "input_sha256": sha256(args.input), "analysis_sha256": sha256(__file__),
              "runtime": {"python": platform.python_version(), "numpy": np.__version__},
              "validation": {"passed": True, "rows": panel["rows"], "n": len(panel["ids"]),
                             "positive": int(panel["y"].sum()), "groups": len(set(panel["groups"])),
                             "seeds": list(panel["seeds"]), "folds": list(panel["folds"]),
                             "scope": "OOF panel integrity; cannot establish training isolation without split audit"},
              "aggregation": {"auc": "per-ID mean score over ten seeds; ties count half",
                              "decision": "at least five positive fold-threshold decisions out of ten",
                              "seed_metrics": "descriptive stability only; not ten independent datasets"},
              "summary": summary, "contrasts": contrasts, "bootstrap": bootmeta,
              "C_D_discordance": changed_counts,
              "limitations": ["post hoc exploratory internal validation, not external or clinical validation",
                              "confidence intervals condition on the stored predictions and do not include refitting",
                              "group bootstrap preserves within-group dependence, not all cross-group homology",
                              "aggregate AUC uses score means; errors use majority votes, different summaries",
                              "finite test performance cannot guarantee R90 recall for new variants"]}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    write_csv(args.out / "summary.csv", [{"arm": a, **summary[a]} for a in ARMS])
    write_csv(args.out / "contrasts.csv", contrasts)
    write_csv(args.out / "seed_metrics.csv", seed_rows)
    write_csv(args.out / "seed_contrasts.csv", seed_differences)
    write_csv(args.out / "aggregated_predictions_LOCAL_ONLY.csv", aggregate_rows)
    write_csv(args.out / "discordant_C_D_LOCAL_ONLY.csv", changed_rows,
              ["VariationID", "y", "group", "C_score", "D_score", "C_decision", "D_decision", "change"])
    print(json.dumps({"status": "PASS", "n": len(panel["ids"]), "summary": summary,
                      "primary_contrasts": [r for r in contrasts if r["role"] == "primary"]}, allow_nan=False))


if __name__ == "__main__":
    main()
