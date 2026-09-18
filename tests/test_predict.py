"""Tests for the predict/ inference bundle (prediction_core.py + bundle/).

These need xgboost/numpy (predict/requirements-inference.txt) and are skipped
if that's not installed — CI installs it explicitly for this module; local
`make test` runs the dependency-free suite only unless you've installed those
extras yourself.
"""

import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "predict"))

xgboost = pytest.importorskip("xgboost")
pytest.importorskip("numpy")

import prediction_core as pc  # noqa: E402

BUNDLE_DIR = REPO_ROOT / "predict" / "bundle"


@pytest.fixture(scope="module")
def bundle():
    # Bundle() re-hashes every file in bundle_manifest.json against its
    # recorded digest, so this alone is a full integrity check of what got
    # committed to the repo.
    return pc.Bundle(BUNDLE_DIR)


def _read_examples():
    path = REPO_ROOT / "predict" / "examples" / "input.csv"
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def test_bundle_has_the_five_documented_models(bundle):
    assert set(bundle.registry["models"]) == set(pc.MODELS)
    assert set(pc.MODELS) == {
        "ORIGINAL18",
        "GPN3",
        "ORIGINAL18_GPN",
        "FULL_UNIFIED",
        "FULL_STRATIFIED",
    }


def test_known_variant_scores_on_all_five_models(bundle):
    rows = _read_examples()
    results = pc.predict_batch(rows, bundle)
    by_id = {r["variant_id"]: r for r in results}
    assert "48" in by_id, "example input.csv is expected to include variant_id 48"
    for name in pc.MODELS:
        model_out = by_id["48"]["models"][name]
        assert model_out["status"] == "scored_research"
        assert model_out["score"] is not None
        assert 0.0 <= model_out["score"] <= 1.0
        # The tool must never claim a calibrated probability.
        assert model_out["probability"] is None
        assert model_out["probability_status"] == "not_calibrated"


def test_regression_scores_match_the_verified_2026_09_15_run(bundle):
    """Pins known scores so a future bundle/code change can't silently drift.

    Values captured from a verified run of the original (pre-refactor)
    prediction_core.py against this exact bundle; the refactor that produced
    the current readable prediction_core.py was checked to give byte-identical
    predictions.jsonl/predictions.csv output on all of examples/input.csv,
    examples/mixed_inputs.csv and the full 7,125-variant development cohort
    before this test was written.
    """
    rows = _read_examples()
    results = pc.predict_batch(rows, bundle)
    by_id = {r["variant_id"]: r for r in results}
    scores = by_id["48"]["models"]
    assert scores["ORIGINAL18"]["score"] == pytest.approx(0.620797336101532, abs=1e-6)
    assert scores["GPN3"]["score"] == pytest.approx(0.8927870392799377, abs=1e-6)
    assert scores["ORIGINAL18_GPN"]["score"] == pytest.approx(0.9997738003730774, abs=1e-6)
    assert scores["FULL_UNIFIED"]["score"] == pytest.approx(0.9991288781166077, abs=1e-6)
    assert scores["FULL_STRATIFIED"]["score"] == pytest.approx(0.998528242111206, abs=1e-6)


def test_variant_not_in_frozen_cache_returns_null_not_a_fabricated_score(bundle):
    rows = [{"variant_id": "novel1", "assembly": "GRCh38", "chrom": "1", "pos": "123456", "ref": "A", "alt": "G"}]
    results = pc.predict_batch(rows, bundle)
    assert results[0]["annotation_status"] == "missing_annotation"
    for name in pc.MODELS:
        assert results[0]["models"][name]["score"] is None
        assert results[0]["models"][name]["status"] == "missing_annotation"


def test_mixed_edge_cases_all_get_accounted_for_not_dropped(bundle):
    path = REPO_ROOT / "predict" / "examples" / "mixed_inputs.csv"
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    results = pc.predict_batch(rows, bundle)
    # Every input row must come back with a row, even the malformed ones.
    assert len(results) == len(rows)
    statuses = {r["annotation_status"] for r in results}
    # At least one row should hit a real edge case, proving this fixture
    # actually exercises error handling and not just the happy path.
    assert statuses - {"verified_cache"}
