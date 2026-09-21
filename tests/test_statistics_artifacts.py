"""Integrity and reference consistency of the versioned statistical revision."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REVISION = ROOT / "results/statistics_20260921"


def test_versioned_outputs_match_manifest():
    manifest = json.loads((REVISION / "manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["outputs_sha256"].items():
        assert hashlib.sha256((REVISION / name).read_bytes()).hexdigest() == expected, name
    assert hashlib.sha256((ROOT / "results/summary.json").read_bytes()).hexdigest() == manifest["reference_summary_sha256"]
    assert manifest["n_boot"] == 2000


def test_revision_preserves_reference_predictions_and_cohorts():
    summary = json.loads((REVISION / "summary.json").read_text(encoding="utf-8"))
    old = json.loads((ROOT / "results/summary.json").read_text(encoding="utf-8"))
    assert summary["training_performed"] is False
    for label, key in [("headline", "pooled_headline"), ("coding_nontruncating", "pooled_missense")]:
        current = summary["cohorts"][label]
        assert current["n_variants"] == old[key]["n_variants"]
        assert current["confusion"] == old[key]["confusion"]
    for label, key in [("headline", "gene_held_out_auc"), ("random_split", "random_split_auc")]:
        point = summary["cohorts"][label]["bootstrap_ci"]["AUC"]["estimate"]
        assert abs(point - old["split_gap"][key]) <= 0.000051


def test_interval_metadata_and_valid_draw_accounting():
    summary = json.loads((REVISION / "summary.json").read_text(encoding="utf-8"))
    for cohort in summary["cohorts"].values():
        ci = cohort["bootstrap_ci"]
        meta = ci["metadata"]
        assert meta["method"] == "gene_group_cluster_percentile"
        assert 1 < meta["n_groups"] <= cohort["n_variants"]
        for key in ["AUC", "R90_recall", "R90_precision"]:
            m = ci[key]
            assert m["n_valid"] + m["n_undefined"] == meta["n_boot"]
            assert m["status"] == "ok"
            assert 0 <= m["estimate"] <= 1
            assert 0 <= m["lo"] <= m["hi"] <= 1
