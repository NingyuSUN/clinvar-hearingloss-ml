"""CSV batch entry point for the five-model prediction core.

Scores are separate, uncalibrated classifier outputs — never averaged across
models, never a calibrated probability. See MODEL_CARD.md and docs/predict.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import time
from pathlib import Path

import numpy as np
import xgboost as xgb

from prediction_core import Bundle, FastaReference, INPUT_FIELDS, MODELS, predict_batch, sha, write_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--annotations")
    ap.add_argument("--annotation-manifest")
    ap.add_argument(
        "--reference",
        help="Optional GRCh38 Ensembl115 FASTA for additional REF checking; "
        "novel annotations remain unsupported in v1",
    )
    args = ap.parse_args()
    started = time.time()
    out = Path(args.output)

    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a fresh output directory")
    if bool(args.annotations) != bool(args.annotation_manifest):
        raise ValueError("External annotations and manifest must be provided together")

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []
        if set(INPUT_FIELDS) - set(columns):
            raise ValueError("Required columns: " + ",".join(INPUT_FIELDS))
        rows = list(reader)
    if not rows:
        raise ValueError("No input rows")

    bundle = Bundle(args.bundle)
    external = bundle.external_annotations(args.annotations, args.annotation_manifest) if args.annotations else {}
    reference = FastaReference(args.reference, bundle.annotation_manifest["reference_manifest"]) if args.reference else None

    try:
        results = predict_batch(rows, bundle, external, reference)
    finally:
        if reference:
            reference.close()

    out.mkdir(parents=True, exist_ok=True)
    (out / "predictions.jsonl").write_text(
        "".join(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n" for r in results)
    )

    longs = []
    wides = []
    common_fields = [
        "input_index", "variant_id", "variant_key", "assembly", "chrom", "pos", "ref", "alt",
        "variant_type", "annotation_status", "group_context", "model_fold",
        "in_development_7125", "in_common_4050",
    ]
    for r in results:
        common = {k: r.get(k) for k in common_fields}
        wide = common.copy()
        for name, m in r["models"].items():
            one = {**common, "model": name, **m}
            for k, v in one.items():
                if isinstance(v, (dict, list)):
                    one[k] = json.dumps(v, ensure_ascii=False, allow_nan=False)
            longs.append(one)
            wide[name + "_score"] = m["score"]
            wide[name + "_status"] = m["status"]
            wide[name + "_decision"] = m["research_decision"]
        wides.append(wide)

    for name, data in [("predictions.csv", longs), ("predictions_wide.csv", wides)]:
        fields = list(dict.fromkeys(k for r in data for k in r))
        with (out / name).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(data)

    counts = {n: {} for n in MODELS}
    for r in results:
        for n, m in r["models"].items():
            counts[n][m["status"]] = counts[n].get(m["status"], 0) + 1

    manifest = {
        "status": "completed_research_prediction",
        "input_rows": len(rows),
        "model_rows": len(longs),
        "input_sha256": sha(args.input),
        "bundle_manifest_sha256": sha(Path(args.bundle) / "bundle_manifest.json"),
        "runtime_code_sha256": {p.name: sha(p) for p in [Path(__file__), Path(__file__).with_name("prediction_core.py")]},
        "annotation_sha256": sha(args.annotations) if args.annotations else None,
        "annotation_manifest_sha256": sha(args.annotation_manifest) if args.annotation_manifest else None,
        "reference_fasta_sha256": (
            bundle.annotation_manifest["reference_manifest"]["files"][Path(args.reference).name]["sha256"]
            if args.reference
            else None
        ),
        "ignored_input_columns": sorted(set(columns) - set(INPUT_FIELDS)),
        "models": counts,
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "xgboost": xgb.__version__},
        "seconds": time.time() - started,
        "limits": [
            "No calibrated probabilities.",
            "This is a fixed seed101 research profile; do not interpret outputs as independent external performance.",
            "Unsupported or missing annotations yield null, not benign predictions.",
            "Novel external annotations are unverified and yield null in v1; source hashes alone do not authorize scoring.",
        ],
        "output_hashes": {p.name: sha(p) for p in out.iterdir() if p.is_file()},
    }
    write_json(out / "run_manifest.json", manifest)
    print(json.dumps({"rows": len(rows), "model_rows": len(longs), "models": counts, "seconds": manifest["seconds"]}))


if __name__ == "__main__":
    main()
