"""Score a genuinely novel variant — one NOT already in predict/bundle's frozen cache.

This is a separate, weaker-guarantee path from predict_variants.py. That
script only ever replays a hash-verified, pre-reviewed cache; this script
calls five live public/keyed data sources at request time (Ensembl VEP,
gnomAD, a UCSC conservation track, a Hugging Face-hosted GPN-Star lookup, and
optionally the AlphaGenome API) and scores whatever comes back. Read
docs/predict_novel.md before trusting its output for anything.

Usage:
    export ALPHAGENOME_API_KEY=...          # optional; omit to skip AVI/splice_sites
    python predict_novel.py --input variants.csv --bundle bundle --output out/
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import xgboost as xgb

sys.path.insert(0, str(Path(__file__).parent))
from prediction_core import (  # noqa: E402
    INPUT_FIELDS,
    MODELS,
    Bundle,
    extract_features,
    normalize_variant,
    prediction_set,
    select_fold,
    sha,
    write_json,
)
import live_sources as ls  # noqa: E402

ORIGIN = "live_annotated_v1"


def annotate_variant(row: dict, api_key: str | None) -> dict:
    """Build one variant's feature dict from live sources, or an explained failure."""
    r = normalize_variant(row)
    result = {"variant_key": r["variant_key"], "error": r["error"], "variant_type": None, "gene_symbols": None, "features": {}, "sources": {}}
    if r["error"]:
        return result

    chrom, pos, ref, alt = r["chrom"], r["pos"], r["ref"], r["alt"]

    try:
        vep = ls.fetch_vep(chrom, pos, ref, alt)
    except ls.AnnotationError as exc:
        result["error"] = f"vep_failed:{exc}"
        return result

    if vep["variant_type"] is None:
        result["error"] = "unsupported_variant_type"
        return result
    gene_symbol = vep["gene_symbol"]
    if not gene_symbol:
        result["error"] = "unresolved_gene_identity"
        return result

    features = dict(vep["features"])
    features.update({
        "variant_length": {"value": 0.0, "status": "observed"},
        "has_del": {"value": 0.0, "status": "observed"},
        "has_ins": {"value": 0.0, "status": "observed"},
        "has_dup": {"value": 0.0, "status": "observed"},
    })
    result["sources"]["transcript_id"] = vep["transcript_id"]
    result["sources"]["transcript_policy"] = vep["transcript_policy"]

    try:
        features.update(ls.fetch_gnomad_af(chrom, pos, ref, alt))
    except ls.AnnotationError as exc:
        result["error"] = f"gnomad_af_failed:{exc}"
        return result

    try:
        features.update(ls.fetch_gnomad_constraint(gene_symbol))
    except ls.AnnotationError as exc:
        result["error"] = f"gnomad_constraint_failed:{exc}"
        return result

    try:
        features.update(ls.fetch_conservation(chrom, pos))
    except ls.AnnotationError as exc:
        features["ensembl_conservation"] = {"value": None, "status": "missing", "error": str(exc)}

    features.update(ls.fetch_gpn_star(chrom, pos, ref, alt))

    if api_key:
        try:
            features.update(ls.fetch_alphagenome(chrom, pos, ref, alt, gene_symbol, api_key))
        except Exception as exc:  # noqa: BLE001 — never let AlphaGenome failure block the other four models
            features["avi"] = {"value": None, "status": "missing", "error": str(exc)[:200]}
            features["splice_sites"] = {"value": None, "status": "missing", "error": str(exc)[:200]}
    else:
        features["avi"] = {"value": None, "status": "not_requested"}
        features["splice_sites"] = {"value": None, "status": "not_requested"}

    result.update(variant_type=vep["variant_type"], gene_symbols=[gene_symbol], features=features)
    return result


def score_annotation(annotation: dict, bundle: Bundle) -> dict:
    """Score one already-annotated variant against all five models.

    Mirrors prediction_core.predict_batch's per-variant scoring logic, but
    for a single live-annotated record instead of a cache hit — deliberately
    reimplemented here rather than calling predict_batch, which is hard-locked
    to the verified cache by design (see docs/predict.md's "closed failure"
    note) and must stay that way.
    """
    models_out = {n: {"status": "not_scored", "score": None, "score_type": "uncalibrated_classifier_score", "probability": None, "probability_status": "not_calibrated", "research_decision": "uncertain_unavailable", "raw_set": None, "reasons": []} for n in MODELS}

    if annotation["error"]:
        for mo in models_out.values():
            mo.update(status=annotation["error"], reasons=[annotation["error"]])
        return models_out

    try:
        fold, context, groups = select_fold(annotation["gene_symbols"], bundle.registry)
    except ValueError as exc:
        for mo in models_out.values():
            mo.update(status=str(exc), reasons=[str(exc)])
        return models_out
    annotation["model_fold"], annotation["group_context"] = fold, context

    route = "protein" if annotation["variant_type"] == "missense" else "RNA"
    fake_ann = {"features": annotation["features"], "origin": ORIGIN}
    for name, mo in models_out.items():
        booster, meta, thresholds = bundle.model(name, fold, route)
        values, errors = extract_features(fake_ann, meta["features"], bundle.continuous)
        mo.update(model_route=meta.get("route", route if name == "FULL_STRATIFIED" else "pooled"), model_weight_sha256=meta["weight_sha256"], annotation_origin=ORIGIN)
        if errors:
            mo.update(status="missing_or_invalid_features", reasons=errors)
            continue
        x = np.asarray([values], dtype=float)
        for j, feat in enumerate(meta["features"]):
            fill = meta["medians"].get(feat, 0.0)
            if fill is not None and np.isnan(x[0, j]):
                x[0, j] = fill
        score = float(booster.predict(xgb.DMatrix(x, feature_names=meta["features"]))[0])
        mo.update(status="scored_research", score=score, **prediction_set(score, annotation["variant_type"], thresholds))
    return models_out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--api-key", help="AlphaGenome API key; falls back to $ALPHAGENOME_API_KEY. Omit both to skip AVI/splice_sites.")
    args = ap.parse_args()
    started = time.time()
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a fresh output directory")

    api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("Warning: no AlphaGenome API key given — avi/splice_sites will be null, and FULL_UNIFIED/FULL_STRATIFIED will report missing_or_invalid_features.", file=sys.stderr)

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []
        if set(INPUT_FIELDS) - set(columns):
            raise ValueError("Required columns: " + ",".join(INPUT_FIELDS))
        rows = list(reader)
    if not rows:
        raise ValueError("No input rows")

    bundle = Bundle(args.bundle)
    results = []
    for index, row in enumerate(rows):
        annotation = annotate_variant(row, api_key)
        models = score_annotation(annotation, bundle)
        results.append({
            "input_index": index,
            "variant_id": row.get("variant_id"),
            "variant_key": annotation["variant_key"],
            "assembly": row.get("assembly"),
            "chrom": row.get("chrom"),
            "pos": row.get("pos"),
            "ref": row.get("ref"),
            "alt": row.get("alt"),
            "annotation_status": ORIGIN if not annotation["error"] else annotation["error"],
            "variant_type": annotation["variant_type"],
            "gene_symbols": annotation["gene_symbols"],
            "model_fold": annotation.get("model_fold"),
            "group_context": annotation.get("group_context"),
            "in_development_7125": False,
            "in_common_4050": False,
            "components": {k: v for k, v in annotation["features"].items() if k in ("gpn_v100", "gpn_m447", "gpn_p243", "avi", "splice_sites")},
            "sources": annotation.get("sources", {}),
            "models": models,
        })

    out.mkdir(parents=True, exist_ok=True)
    (out / "predictions.jsonl").write_text("".join(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n" for r in results))

    manifest = {
        "status": "completed_live_annotated_prediction",
        "origin": ORIGIN,
        "input_rows": len(rows),
        "input_sha256": sha(args.input),
        "bundle_manifest_sha256": sha(Path(args.bundle) / "bundle_manifest.json"),
        "alphagenome_requested": bool(api_key),
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "xgboost": xgb.__version__},
        "seconds": time.time() - started,
        "limits": [
            "This path calls live external services at request time and is far less verified than predict_variants.py's frozen-cache path.",
            "ensembl_conservation uses UCSC phyloP100way as a documented substitute for the original (unaudited) Ensembl GERP-style score — different scale, not a reproduction.",
            "The VEP transcript-selection 'target-gene' preference tier is not reproduced (the original gene panel list was not preserved); MANE Select > canonical > protein-coding > any is used instead.",
            "No calibrated probabilities. No external validation of this live path exists yet — see PROJECT_STATUS.md.",
        ],
    }
    write_json(out / "run_manifest.json", manifest)
    print(json.dumps({"rows": len(rows), "seconds": manifest["seconds"], "alphagenome_requested": manifest["alphagenome_requested"]}))


if __name__ == "__main__":
    main()
