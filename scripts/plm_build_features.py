#!/usr/bin/env python3
"""Turn the Kaggle ESM log-probabilities into per-variant PLM features.

    python scripts/plm_build_features.py --scores data/plm/esm_scores

Label-blind. Checks that the scoring output came from the committed scoring
input, that every scored position has all six pinned models, and that the
kernel's reference residue matches the mapping. Writes
``data/plm/plm_features.csv`` and ``data/plm/features_manifest.json``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from hlpath.config import REPO_ROOT
from hlpath.plm import (ESM1V_MODELS, ESM2_MODEL, MODEL_REVISIONS, PROTOCOL_ID, SCORED_FLAG,
                        STATUS_SCORED, stable_id, variant_llrs)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path, default=REPO_ROOT / "data/plm/esm_scores")
    ap.add_argument("--plm-dir", type=Path, default=REPO_ROOT / "data/plm")
    args = ap.parse_args()

    scoring_input = args.plm_dir / "esm_scoring_input.json"
    mapping_path = args.plm_dir / "variant_protein_map.csv"
    logprob_path = args.scores / "esm_logprobs.csv.gz"
    kernel_manifest = json.loads((args.scores / "esm_scoring_manifest.json").read_text())

    if kernel_manifest["input_sha256"] != sha256(scoring_input):
        raise SystemExit("ESM output was not produced from the committed scoring input")
    if kernel_manifest["output_sha256"] != sha256(logprob_path):
        raise SystemExit("ESM log-probability file does not match the kernel manifest")
    if kernel_manifest.get("limit") is not None:
        raise SystemExit("smoke-test output cannot be used as protocol features")
    resolved = {m["repo"]: m["resolved_revision"] for m in kernel_manifest["models"]}
    if resolved != MODEL_REVISIONS:
        raise SystemExit(f"model revisions differ from the protocol: {resolved}")

    mapping = pd.read_csv(mapping_path, low_memory=False)
    logprobs = pd.read_csv(logprob_path)
    scored = mapping[mapping["map_status"] == STATUS_SCORED]
    want = set(zip(scored["cached_transcript_id"].map(stable_id), scored["protein_position"].astype(int)))
    for model in (*ESM1V_MODELS, ESM2_MODEL):
        got = logprobs[logprobs["model"] == model]
        have = set(zip(got["transcript_stable_id"], got["protein_position"].astype(int)))
        if have != want:
            raise SystemExit(f"{model}: scored positions differ from the mapping "
                             f"({len(want - have)} missing, {len(have - want)} extra)")
    ref_check = scored.assign(transcript_stable_id=scored["cached_transcript_id"].map(stable_id)).merge(
        logprobs[["transcript_stable_id", "protein_position", "ref_aa"]].drop_duplicates(),
        on=["transcript_stable_id", "protein_position"], suffixes=("", "_kernel"))
    if (ref_check["ref_aa"] != ref_check["ref_aa_kernel"]).any():
        raise SystemExit("kernel reference residue differs from the mapping")

    features = variant_llrs(mapping, logprobs)
    out_path = args.plm_dir / "plm_features.csv"
    features.to_csv(out_path, index=False, lineterminator="\n", float_format="%.6f")
    manifest = {
        "protocol": PROTOCOL_ID,
        "label_column_read": False,
        "inputs_sha256": {
            "data/plm/variant_protein_map.csv": sha256(mapping_path),
            "data/plm/esm_scoring_input.json": sha256(scoring_input),
            "data/plm/esm_scores/esm_logprobs.csv.gz": sha256(logprob_path),
            "data/plm/esm_scores/esm_scoring_manifest.json": sha256(args.scores / "esm_scoring_manifest.json"),
        },
        "outputs_sha256": {"data/plm/plm_features.csv": sha256(out_path)},
        "counts": {
            "variants": int(len(features)),
            "plm_scored": int(features[SCORED_FLAG].sum()),
            "esm1v_llr_present": int(features["esm1v_llr"].notna().sum()),
            "esm2_llr_present": int(features["esm2_llr"].notna().sum()),
        },
        "kernel": {k: kernel_manifest[k] for k in ["generated_utc", "strategy", "versions"]},
    }
    (args.plm_dir / "features_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
