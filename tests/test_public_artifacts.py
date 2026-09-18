"""Deterministic, fast checks on the committed public artifacts.

These do not retrain anything and do not require the optional `interpret`
extras (shap/matplotlib) — they only check that what is already committed
under results/ and docs/ is well-formed. See tools/validate_public_artifacts.py
for the checks themselves; this module just exercises them under pytest.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

import validate_public_artifacts as vpa  # noqa: E402


def test_no_crlf_in_tracked_text_artifacts():
    problems: list[str] = []
    vpa.check_line_endings(problems)
    assert problems == []


def test_summary_json_has_required_keys():
    problems: list[str] = []
    summary = vpa.check_json_file(
        problems, "results/summary.json", vpa.REQUIRED_SUMMARY_KEYS
    )
    assert problems == []
    assert summary is not None


def test_run_manifest_has_required_keys_and_valid_hash():
    problems: list[str] = []
    manifest = vpa.check_json_file(
        problems, "results/run_manifest.json", vpa.REQUIRED_MANIFEST_KEYS
    )
    vpa.check_manifest_hash(problems, manifest)
    assert problems == []


def test_result_csvs_are_well_formed():
    problems: list[str] = []
    vpa.check_result_csvs(problems)
    assert problems == []


def test_no_pycache_tracked_by_git():
    problems: list[str] = []
    vpa.check_no_stray_bytecode(problems)
    assert problems == []


def test_full_validator_passes_end_to_end():
    assert vpa.main() == 0
