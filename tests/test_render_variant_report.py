"""Tests for predict/render_variant_report.py — the human-readable report layer.

These only need stdlib (json/argparse/pathlib), so they run in every
environment, unlike tests/test_predict.py which needs xgboost/numpy installed.
"""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "predict" / "render_variant_report.py"
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "sample_predictions.jsonl"


def render(variant_id=None):
    args = [sys.executable, str(SCRIPT), "--predictions", str(FIXTURE)]
    if variant_id:
        args += ["--variant-id", variant_id]
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    return result.stdout


def test_scored_variant_reports_all_five_models():
    text = render("48")
    for label in [
        "18 features only",
        "GPN only (3 scores)",
        "18 features + GPN",
        "18 features + GPN + AVI (pooled)",
        "18 features + GPN + AVI (stratified)",
    ]:
        assert label in text


def test_report_never_calls_the_score_a_probability():
    text = render("48")
    assert "not probabilities of pathogenicity" in text
    # The literal word "probability" should only appear in that disclaimer,
    # never attached to a percentage-style claim about this variant.
    assert "% chance" not in text.replace('"X% chance', "")


def test_avi_is_labelled_as_composite_not_pure_alphagenome():
    text = render("48")
    assert "NOT a pure AlphaGenome score" in text
    assert "composite AlphaGenome/AlphaMissense-family score" in text


def test_missing_annotation_variant_explains_why_not_scored():
    text = render("novel_demo")
    assert "missing_annotation" in text
    assert "only replays known, pre-annotated ClinVar hearing-loss variants" in text
    # No model should show a fabricated score for an unscored variant.
    assert "0.62" not in text


def test_render_all_variants_when_no_filter_given():
    text = render()
    assert "Variant `48`" in text
    assert "Variant `novel_demo`" in text
