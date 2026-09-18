"""Versioned, label-free variant scoring against a frozen, hash-verified bundle.

No model fitting and no remote calls happen at inference time. Every score
comes from replaying an already-trained XGBoost booster on an already-verified
cached annotation; a variant that is not already in the frozen cache (or that
arrives via an unverified external annotation) scores as null on all five
models, by design (see MODEL_CARD.md and docs/predict.md).
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import xgboost as xgb

MODELS = ["ORIGINAL18", "GPN3", "ORIGINAL18_GPN", "FULL_UNIFIED", "FULL_STRATIFIED"]
TYPES = {"missense", "canonical_splice", "splice_region", "synonymous"}
INPUT_FIELDS = ["variant_id", "assembly", "chrom", "pos", "ref", "alt"]
FORBIDDEN = {
    "y",
    "label",
    "clinical_label",
    "clinical_significance",
    "clinvar_classification",
    "ReviewStatus",
}


def sha(path):
    """Return the sha256 hex digest of a file, read in 1 MiB chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(p):
    return json.loads(Path(p).read_text())


def write_json(p, obj):
    Path(p).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def normalize_variant(row):
    """Parse and validate one input row into a canonical variant record.

    Sets `error` to a short machine-readable reason and returns early as soon
    as a field fails to parse; `variant_key` stays None until the identity is
    fully resolved.
    """
    r = {k: str(row.get(k, "") or "").strip() for k in INPUT_FIELDS}
    r["error"] = None
    r["variant_key"] = None

    if not r["variant_id"]:
        r["error"] = "missing_variant_id"
        return r
    if r["assembly"] != "GRCh38":
        r["error"] = "unsupported_assembly"
        return r

    chrom = re.sub("^chr", "", r["chrom"], flags=re.I).upper()
    if chrom == "M":
        chrom = "MT"
    if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
        r["error"] = "unsupported_chromosome"
        return r
    if not re.fullmatch("[1-9][0-9]*", r["pos"]):
        r["error"] = "invalid_position"
        return r

    r.update(chrom=chrom, pos=int(r["pos"]), ref=r["ref"].upper(), alt=r["alt"].upper())
    if not re.fullmatch("[ACGT]+", r["ref"]) or not re.fullmatch("[ACGT]+", r["alt"]):
        r["error"] = "invalid_allele"
        return r

    r["variant_key"] = f"{chrom}-{r['pos']}-{r['ref']}-{r['alt']}"
    if r["ref"] == r["alt"]:
        r["error"] = "identical_alleles"
    elif len(r["ref"]) != 1 or len(r["alt"]) != 1:
        r["error"] = "unsupported_non_snv"
    return r


def normalize_rows(rows):
    """Normalize every row and flag duplicate variant_ids / variant_keys."""
    out = [normalize_variant(r) for r in rows]
    ids = Counter(r["variant_id"] for r in out if r["variant_id"])
    keys = Counter(r["variant_key"] for r in out if r["variant_key"])
    for r in out:
        if ids[r["variant_id"]] > 1:
            r["error"] = "duplicate_variant_id"
        elif r["variant_key"] and keys[r["variant_key"]] > 1:
            r["error"] = "duplicate_variant_key"
    return out


def reject_labels(obj):
    """Raise if a ClinVar-style outcome label is present anywhere in `obj`.

    Applied to every annotation record (cached or externally supplied) so
    that a clinical label can never leak into a runtime feature vector.
    """
    if isinstance(obj, dict):
        if FORBIDDEN.intersection(obj):
            raise ValueError("Outcome labels are forbidden in annotations")
        for value in obj.values():
            reject_labels(value)
    elif isinstance(obj, list):
        for value in obj:
            reject_labels(value)


def extract_features(annotation, names, continuous):
    """Pull a model's declared feature vector out of one annotation record.

    Returns `(values, [])` on success or `(None, errors)` listing which named
    features were missing/invalid. A feature frozen as legitimately missing in
    the verified cache becomes NaN (later filled by the model's frozen median);
    anything else missing or malformed is a hard error, not a silent zero.
    """
    result = []
    errors = []
    for name in names:
        feature = annotation.get("features", {}).get(name)
        if not isinstance(feature, dict):
            errors.append(name + ":missing_feature")
            continue
        value, status = feature.get("value"), feature.get("status")
        if value is None and status == "frozen_missing" and annotation.get("origin") == "verified_cache":
            result.append(float("nan"))
            continue
        if status != "observed" or value is None:
            errors.append(name + ":" + str(status or "missing_status"))
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            errors.append(name + ":invalid_numeric")
            continue
        if name not in continuous and value not in (0, 1):
            errors.append(name + ":invalid_binary")
            continue
        result.append(float(value))
    return (None, errors) if errors else (result, [])


def select_fold(genes, registry):
    """Pick the gene-held-out fold that never trained/selected/calibrated on `genes`.

    Returns `(fold, context, groups)`. Unknown genes get fold 0 under
    `novel_gene_context` (there is no held-out guarantee to give them);
    genes that resolve to more than one historical gene-group, or that mix
    known and unrepresented genes, are rejected rather than guessed at.
    """
    if not isinstance(genes, list) or not genes or any(not isinstance(g, str) or not g.strip() for g in genes):
        raise ValueError("unresolved_gene_identity")

    known = {registry["gene_to_group"][g] for g in genes if g in registry["gene_to_group"]}
    if len(known) > 1:
        raise ValueError("unsupported_group_mapping")
    if not known:
        return 0, "novel_gene_context", []

    group = next(iter(known))
    if group not in registry["group_to_fold"]:
        raise ValueError("unsupported_group_mapping")
    # Combining known and unrepresented genes changes historical connected-group semantics.
    if any(g not in registry["gene_to_group"] for g in genes):
        raise ValueError("unsupported_group_mapping")
    return int(registry["group_to_fold"][group]), "known_gene_held_out", [group]


def prediction_set(score, variant_type, thresholds):
    """Turn a raw score into a frozen-threshold research decision.

    `thresholds` holds two rows (label 0 and label 1) per consequence class,
    each a false-positive-rate-style quantile at alpha=0.05 (or `unbounded`
    if that class never had enough calibration data). A directional research
    call is only ever reported for missense; every other supported subtype
    is deliberately reported as `uncertain_subtype_validation` regardless of
    its raw set, because those subtypes' directional calls are not validated.
    """
    if not math.isfinite(score) or not 0 <= score <= 1:
        raise ValueError("Invalid score")

    recs = [r for r in thresholds if r["cell"] == variant_type]
    if len(recs) != 2 or {r["label"] for r in recs} != {0, 1}:
        raise ValueError("Incomplete calibration cell")

    p = float(np.float32(score))
    include = {}
    unbounded = False
    for r in recs:
        q, u = r["q"], r["unbounded"]
        if type(u) is not bool or u != (q is None):
            raise ValueError("Invalid unbounded threshold")
        if q is not None and (not math.isfinite(q) or not 0 <= q <= 1):
            raise ValueError("Invalid threshold")
        include[r["label"]] = True if u else (p if r["label"] == 0 else 1.0 - p) <= q
        unbounded |= u

    raw_set = [label for label in (0, 1) if include[label]]
    decision = {
        (0,): "lower_risk_support",
        (1,): "higher_risk_support",
        (0, 1): "uncertain_both",
        (): "uncertain_empty",
    }[tuple(raw_set)]
    reason = (
        "insufficient_calibration_groups"
        if unbounded
        else (
            "overlapping_sets"
            if len(raw_set) == 2
            else ("empty_conflict" if not raw_set else "single_label_set")
        )
    )
    return {
        "raw_set": raw_set,
        "raw_set_decision": decision,
        "research_decision": decision if variant_type == "missense" else "uncertain_subtype_validation",
        "decision_reason": reason if variant_type == "missense" else "subtype_directional_calls_not_validated",
        "calibration_alpha": 0.05,
    }


class FastaReference:
    """Optional hash-verified GRCh38 FASTA lookup for extra REF checking.

    Never used to derive features for novel variants — v1 has no independently
    reconstructed feature-extraction pipeline for anything outside the frozen
    cache, so this class only ever confirms or contradicts a REF allele.
    """

    def __init__(self, path, manifest):
        self.path = Path(path)
        fai = Path(str(path) + ".fai")
        self.manifest = manifest
        for p in (self.path, fai):
            entry = manifest["files"][p.name]
            if p.stat().st_size != entry["bytes"] or sha(p) != entry["sha256"]:
                raise ValueError("Reference file hash mismatch")
        self.index = {}
        for line in fai.read_text().splitlines():
            contig, length, offset, bases, width = line.split("\t")[:5]
            self.index[contig] = tuple(map(int, (length, offset, bases, width)))
        self.handle = self.path.open("rb")

    def base(self, chrom, pos):
        if chrom not in self.index:
            raise ValueError("unsupported_reference_contig")
        length, offset, bases, width = self.index[chrom]
        if not 1 <= pos <= length:
            raise ValueError("position_out_of_range")
        n = pos - 1
        self.handle.seek(offset + (n // bases) * width + n % bases)
        return self.handle.read(1).decode().upper()

    def close(self):
        self.handle.close()


class Bundle:
    """A loaded, hash-verified inference bundle: models + frozen annotation cache.

    Every file listed in `bundle_manifest.json` is re-hashed against its
    recorded digest at load time, and every cached annotation is checked for
    a `verified_cache` origin and rejected if it contains an outcome label.
    A tampered or partial bundle fails loudly here, not silently at scoring
    time.
    """

    def __init__(self, root):
        self.root = Path(root).resolve()
        manifest = read_json(self.root / "bundle_manifest.json")
        for rel, expected in manifest["files"].items():
            p = self.root / rel
            if not p.resolve().is_relative_to(self.root) or sha(p) != expected:
                raise ValueError("Bundle integrity mismatch: " + rel)

        self.registry = read_json(self.root / "model_registry.json")
        self.annotation_manifest = read_json(self.root / "annotation_manifest.json")
        if set(self.registry["models"]) != set(MODELS):
            raise ValueError("Unexpected model roster")

        self.cache = {}
        self.reference_by_locus = {}
        self.boosters = {}
        for line in (self.root / "annotation_cache.jsonl").read_text().splitlines():
            ann = json.loads(line)
            reject_labels(ann)
            key = ann["variant_key"]
            if key in self.cache:
                raise ValueError("Duplicate cache key")
            if ann.get("origin") != "verified_cache":
                raise ValueError("Wrong cache origin")
            self.cache[key] = ann

            chrom, pos, ref, alt = key.split("-")
            if len(ref) == 1 and len(alt) == 1 and ann["reference_status"] == "matched":
                locus = (chrom, int(pos))
                if locus in self.reference_by_locus and self.reference_by_locus[locus] != ref:
                    raise ValueError("Conflicting reference cache")
                self.reference_by_locus[locus] = ref

        self.continuous = set(self.registry["feature_contract"]["continuous_names"])

    def model(self, name, fold, route):
        folds = [f for f in self.registry["models"][name]["folds"] if f["fold"] == fold]
        if len(folds) != 1:
            raise ValueError("Missing/duplicate model fold")
        fold_meta = folds[0]
        key = route if name == "FULL_STRATIFIED" else "pooled"
        meta = fold_meta["routes"][key]

        rel = meta["weight_path"]
        path = self.root / rel
        if rel not in self.boosters:
            if not path.resolve().is_relative_to(self.root) or sha(path) != meta["weight_sha256"]:
                raise ValueError("Weight hash mismatch")
            booster = xgb.Booster()
            booster.load_model(path)
            booster.set_param({"nthread": 1})
            self.boosters[rel] = booster
        return self.boosters[rel], meta, fold_meta["thresholds"]

    def external_annotations(self, path, manifest_path):
        """Load caller-supplied annotations, but only as a cross-check.

        A caller-supplied source name and hash establish internal consistency,
        not provenance: any record that doesn't already match a frozen cache
        entry is marked `unverified_external_annotation` and scores null on
        every model. There is no independently reconstructed adapter in this
        release for scoring a genuinely novel variant end to end.
        """
        manifest = read_json(manifest_path)
        if manifest["annotation_file_sha256"] != sha(path):
            raise ValueError("External annotation hash mismatch")
        if manifest["feature_contract_id"] != self.annotation_manifest["feature_contract_id"]:
            raise ValueError("Incompatible external feature contract")

        sources = manifest.get("source_artifacts", {})
        base = Path(manifest_path).resolve().parent
        for source, rec in sources.items():
            if sha(base / rec["path"]) != rec["sha256"]:
                raise ValueError("External source artifact hash mismatch")

        out = {}
        for line in Path(path).read_text().splitlines():
            a = json.loads(line)
            reject_labels(a)
            key = a["variant_key"]
            parts = key.split("-")
            if len(parts) != 4:
                raise ValueError("Malformed annotation key")

            r = normalize_variant(dict(zip(INPUT_FIELDS, ["annotation", a.get("assembly", ""), *parts])))
            if r["error"] or r["variant_key"] != key:
                raise ValueError("Noncanonical external key")
            if key in out:
                raise ValueError("Duplicate external annotation")
            if a.get("annotation_source_id") not in sources:
                raise ValueError("Missing annotation identity source")
            if not a.get("gene_symbols") or a.get("variant_type") not in TYPES:
                raise ValueError("Invalid annotation context")

            if key in self.cache:
                cached = self.cache[key]
                if a["gene_symbols"] != cached["gene_symbols"] or a["variant_type"] != cached["variant_type"]:
                    raise ValueError("External annotation conflicts with frozen cache")
                for name, rec in a.get("features", {}).items():
                    if name not in cached["features"] or rec.get("value") != cached["features"][name]["value"]:
                        raise ValueError("External feature conflicts with frozen cache")
                out[key] = cached
                continue

            a["origin"] = "unverified_external_annotation"
            a["in_development_7125"] = False
            a["in_common_4050"] = False
            out[key] = a
        return out


def empty_model():
    return {
        "status": "not_scored",
        "score": None,
        "score_type": "uncalibrated_classifier_score",
        "probability": None,
        "probability_status": "not_calibrated",
        "research_decision": "uncertain_unavailable",
        "raw_set": None,
        "reasons": [],
    }


def predict_batch(rows, bundle, external=None, reference=None):
    """Score every input row against all five models.

    One row's failure (bad identity, missing annotation, unresolved gene
    group, ...) never drops the row: it comes back with a per-model status
    explaining why, never a silent zero or a fabricated benign call.
    """
    requests = normalize_rows(rows)
    outputs = []
    pending = defaultdict(list)
    external = external or {}

    for index, r in enumerate(requests):
        result = {k: v for k, v in r.items() if k != "error"}
        result.update(
            input_index=index,
            profile_id=bundle.registry.get("profile_id", "seed101_research_v1"),
            annotation_status=None,
            variant_type=None,
            gene_symbols=None,
            model_fold=None,
            group_context=None,
            in_development_7125=None,
            in_common_4050=None,
            components={},
            models={n: empty_model() for n in MODELS},
        )
        outputs.append(result)
        error = r["error"]
        a = None

        if not error:
            expected = bundle.reference_by_locus.get((r["chrom"], r["pos"]))
            if reference:
                try:
                    expected = reference.base(r["chrom"], r["pos"])
                except ValueError as e:
                    error = str(e)
            if not error and expected is not None and expected != r["ref"]:
                error = "reference_mismatch"
            a = bundle.cache.get(r["variant_key"]) or external.get(r["variant_key"])
            if not error and a is None:
                error = "missing_annotation"
            if not error and a is not bundle.cache.get(r["variant_key"]):
                error = "unverified_external_annotation"

        if not error and a.get("origin") == "verified_cache" and a.get("reference_status") != "matched":
            error = "reference_" + str(a.get("reference_status", "unverified"))
        if not error and a["variant_type"] not in TYPES:
            error = "unsupported_variant_type"

        if error:
            result["annotation_status"] = error
            for mo in result["models"].values():
                mo.update(status=error, reasons=[error])
            continue

        result.update(
            annotation_status=a["origin"],
            variant_type=a["variant_type"],
            gene_symbols=a["gene_symbols"],
            in_development_7125=a["in_development_7125"],
            in_common_4050=a["in_common_4050"],
        )
        result["components"] = {
            n: a["features"][n]
            for n in ("gpn_v100", "gpn_m447", "gpn_p243", "avi", "splice_sites")
            if n in a["features"]
        }

        try:
            fold, context, groups = select_fold(a["gene_symbols"], bundle.registry)
        except ValueError as e:
            for mo in result["models"].values():
                mo.update(status=str(e), reasons=[str(e)])
            continue
        result.update(model_fold=fold, group_context=context)

        if a.get("gene_group") and groups and a["gene_group"] != groups[0]:
            raise ValueError("Cache group identity mismatch")

        route = "protein" if a["variant_type"] == "missense" else "RNA"
        for name, mo in result["models"].items():
            booster, meta, thresholds = bundle.model(name, fold, route)
            if any(g not in meta["excluded_groups"] for g in groups):
                raise ValueError("Selected model saw held-out gene group")
            values, errors = extract_features(a, meta["features"], bundle.continuous)
            mo.update(
                model_route=meta.get("route", route if name == "FULL_STRATIFIED" else "pooled"),
                model_weight_sha256=meta["weight_sha256"],
                annotation_origin=a["origin"],
            )
            if errors:
                mo.update(status="missing_or_invalid_features", reasons=errors)
                continue
            mo["imputed_features"] = [n for n, v in zip(meta["features"], values) if math.isnan(v)]
            pending[(name, fold, route)].append((values, mo, a["variant_type"]))

    for (name, fold, route), jobs in pending.items():
        booster, meta, thresholds = bundle.model(name, fold, route)
        x = np.asarray([v for v, _, _ in jobs], dtype=float)
        for j, feat in enumerate(meta["features"]):
            fill = meta["medians"].get(feat, 0.0)
            if fill is not None:
                x[np.isnan(x[:, j]), j] = fill
        scores = booster.predict(xgb.DMatrix(x, feature_names=meta["features"])).astype(np.float32)
        for (_, mo, variant_type), score in zip(jobs, scores):
            mo.update(status="scored_research", score=float(score), **prediction_set(float(score), variant_type, thresholds))

    return outputs
