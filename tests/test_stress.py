"""Load/stress tests for the predict/ batch pipeline at scale.

tests/test_predict.py checks correctness on a handful of example rows. This
file checks that the same contracts (every row accounted for, no cross-row
contamination, no silently dropped errors) still hold when the input is tens
of thousands of rows instead of a handful, and that runtime stays roughly
linear rather than blowing up.

Needs xgboost/numpy (predict/requirements-inference.txt), same as
test_predict.py; skipped if that's not installed.
"""

from __future__ import annotations

import csv
import json
import random
import subprocess
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "predict"))

xgboost = pytest.importorskip("xgboost")
pytest.importorskip("numpy")

import prediction_core as pc  # noqa: E402

BUNDLE_DIR = REPO_ROOT / "predict" / "bundle"

# Real GRCh38 chromosome positions never reach this; used to build synthetic
# rows that are guaranteed not to collide with an actual frozen-cache key.
POS_OFFSET = 900_000_000


@pytest.fixture(scope="module")
def bundle():
    return pc.Bundle(BUNDLE_DIR)


def _valid_rows(bundle):
    """One row per distinct cached variant_key, each with a unique variant_id."""
    rows = []
    for i, key in enumerate(bundle.cache):
        chrom, pos, ref, alt = key.split("-")
        rows.append(
            {"variant_id": f"valid{i}", "assembly": "GRCh38", "chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
        )
    return rows


@pytest.fixture(scope="module")
def scorable_ids(bundle):
    """variant_ids of cached rows that actually score (the cache also holds non-SNVs etc.)."""
    return {r["variant_id"] for r in pc.predict_batch(_valid_rows(bundle), bundle) if r["annotation_status"] == "verified_cache"}


def _malformed_rows(n, seed=0):
    """A cycle of invalid/edge-case rows, none of which can collide with a cached key."""
    kinds = [
        "bad_assembly",
        "bad_chrom",
        "bad_pos_text",
        "bad_ref_base",
        "non_snv",
        "missing_id",
        "huge_pos",
        "identical_alleles",
        "well_formed_but_absent",
    ]
    rows = []
    for i in range(n):
        pos = str(POS_OFFSET + i)
        kind = kinds[i % len(kinds)]
        if kind == "bad_assembly":
            row = {"variant_id": f"bad{i}", "assembly": "GRCh37", "chrom": "1", "pos": pos, "ref": "A", "alt": "G"}
        elif kind == "bad_chrom":
            row = {"variant_id": f"bad{i}", "assembly": "GRCh38", "chrom": "ZZ", "pos": pos, "ref": "A", "alt": "G"}
        elif kind == "bad_pos_text":
            row = {"variant_id": f"bad{i}", "assembly": "GRCh38", "chrom": "1", "pos": "not_a_number", "ref": "A", "alt": "G"}
        elif kind == "bad_ref_base":
            row = {"variant_id": f"bad{i}", "assembly": "GRCh38", "chrom": "1", "pos": pos, "ref": "N", "alt": "G"}
        elif kind == "non_snv":
            row = {"variant_id": f"bad{i}", "assembly": "GRCh38", "chrom": "1", "pos": pos, "ref": "AT", "alt": "GC"}
        elif kind == "missing_id":
            row = {"variant_id": "", "assembly": "GRCh38", "chrom": "1", "pos": pos, "ref": "A", "alt": "G"}
        elif kind == "huge_pos":
            row = {"variant_id": f"bad{i}", "assembly": "GRCh38", "chrom": "1", "pos": str(10**18 + i), "ref": "A", "alt": "G"}
        elif kind == "identical_alleles":
            row = {"variant_id": f"bad{i}", "assembly": "GRCh38", "chrom": "1", "pos": pos, "ref": "A", "alt": "A"}
        else:  # well-formed and normalizable, but nowhere in the frozen cache
            row = {"variant_id": f"bad{i}", "assembly": "GRCh38", "chrom": "chr2", "pos": pos, "ref": "a", "alt": "g"}
        rows.append(row)
    random.Random(seed).shuffle(rows)
    return rows


def test_large_mixed_batch_accounts_for_every_row(bundle, scorable_ids):
    valid = _valid_rows(bundle)
    malformed = _malformed_rows(15_000)
    rows = valid + malformed
    random.Random(42).shuffle(rows)

    started = time.time()
    results = pc.predict_batch(rows, bundle)
    elapsed = time.time() - started

    assert len(results) == len(rows), "every input row must come back, none dropped under volume"
    assert elapsed < 30, f"predict_batch took {elapsed:.1f}s for {len(rows)} rows; expected roughly linear scaling"

    by_status = {}
    for r in results:
        assert set(r["models"]) == set(pc.MODELS)
        for mo in r["models"].values():
            assert mo["status"] != "not_scored", "no row's models may be left unaccounted for"
        by_status[r["annotation_status"]] = by_status.get(r["annotation_status"], 0) + 1

    assert by_status.get("verified_cache") == len(scorable_ids)
    assert sum(by_status.values()) == len(rows)


def test_scores_identical_whether_scored_alone_or_buried_in_a_large_batch(bundle, scorable_ids):
    valid = [r for r in _valid_rows(bundle) if r["variant_id"] in scorable_ids]
    sample = valid[::300]
    assert len(sample) >= 10

    padding = _malformed_rows(12_000)
    interleaved = padding[:6_000] + sample + padding[6_000:]

    batch_results = {r["variant_id"]: r for r in pc.predict_batch(interleaved, bundle)}
    solo_results = {r["variant_id"]: r for r in pc.predict_batch(sample, bundle)}

    checked = 0
    for row in sample:
        vid = row["variant_id"]
        for name in pc.MODELS:
            batched = batch_results[vid]["models"][name]
            solo = solo_results[vid]["models"][name]
            assert batched["status"] == solo["status"] == "scored_research"
            assert batched["score"] == pytest.approx(solo["score"], abs=1e-9)
            checked += 1
    assert checked == len(sample) * len(pc.MODELS)


def test_duplicate_ids_and_keys_flagged_at_scale(bundle):
    valid = _valid_rows(bundle)
    dup_id_rows = [dict(r, variant_id="dup_id") for r in valid[:50]]
    dup_key_rows = [dict(r, variant_id=f"dup_key_{i}_{copy}") for copy in (0, 1) for i, r in enumerate(valid[50:100])]
    rows = valid[100:] + dup_id_rows + dup_key_rows

    results = pc.predict_batch(rows, bundle)
    assert len(results) == len(rows)

    dup_id_statuses = [r["annotation_status"] for r in results if r["variant_id"] == "dup_id"]
    assert dup_id_statuses == ["duplicate_variant_id"] * 50

    dup_key_statuses = [r["annotation_status"] for r in results if r["variant_id"].startswith("dup_key_")]
    assert dup_key_statuses == ["duplicate_variant_key"] * 100


def test_cli_handles_large_batch_end_to_end(bundle, scorable_ids, tmp_path):
    valid = _valid_rows(bundle)[:3_500]
    malformed = _malformed_rows(1_500)
    rows = valid + malformed
    random.Random(7).shuffle(rows)

    input_csv = tmp_path / "stress_input.csv"
    with open(input_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=pc.INPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    output_dir = tmp_path / "out"

    started = time.time()
    result = subprocess.run(
        [
            sys.executable, str(REPO_ROOT / "predict" / "predict_variants.py"),
            "--input", str(input_csv), "--bundle", str(BUNDLE_DIR), "--output", str(output_dir),
        ],
        capture_output=True, text=True, timeout=120,
    )
    elapsed = time.time() - started
    assert result.returncode == 0, result.stderr
    assert elapsed < 60, f"predict_variants.py CLI took {elapsed:.1f}s for {len(rows)} rows"

    predictions = (output_dir / "predictions.jsonl").read_text().splitlines()
    assert len(predictions) == len(rows)

    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["input_rows"] == len(rows)
    assert manifest["model_rows"] == len(rows) * len(pc.MODELS)
    expected_scored = sum(r["variant_id"] in scorable_ids for r in valid)
    for name in pc.MODELS:
        assert sum(manifest["models"][name].values()) == len(rows)
        assert manifest["models"][name].get("scored_research", 0) == expected_scored
