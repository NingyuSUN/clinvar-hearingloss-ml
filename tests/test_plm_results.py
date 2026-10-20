"""Integrity and protocol conformance of the committed PLM-ESM-2026-10 results."""
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from hlpath import plm  # noqa: E402

RESULTS = ROOT / "results/plm_esm_2026-10"
PLM = ROOT / "data/plm"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_esm_scores_and_features_match_their_manifests():
    kernel = load(PLM / "esm_scores/esm_scoring_manifest.json")
    assert kernel["input_sha256"] == sha(PLM / "esm_scoring_input.json")
    assert kernel["output_sha256"] == sha(PLM / "esm_scores/esm_logprobs.csv.gz")
    assert kernel["limit"] is None and kernel["embedded_input"] is True
    assert {m["repo"]: m["resolved_revision"] for m in kernel["models"]} == plm.MODEL_REVISIONS
    for m in kernel["models"]:
        assert m["weight_sha256"] == m["hub_lfs_sha256"]
    feats = load(PLM / "features_manifest.json")
    for rel, digest in {**feats["inputs_sha256"], **feats["outputs_sha256"]}.items():
        assert sha(ROOT / rel) == digest, rel


def test_results_match_manifest_and_seal():
    manifest = load(RESULTS / "manifest.json")
    for name, digest in manifest["outputs_sha256"].items():
        assert sha(RESULTS / name) == digest, name
    assert manifest["features_sha256"] == sha(PLM / "plm_features.csv")
    seal = manifest["prediction_seal"]
    assert seal["protocol"] == plm.PROTOCOL_ID
    assert seal["feature_sets"] == plm.PLM_FEATURE_SETS
    assert seal["seeds"] == list(range(101, 111)) and seal["k_outer"] == 5
    assert seal["inputs_sha256"]["data/plm/plm_features.csv"] == sha(PLM / "plm_features.csv")
    assert seal["inputs_sha256"]["data/modeling_matrix.csv.gz"] == sha(ROOT / "data/modeling_matrix.csv.gz")


def test_prespecified_contrasts_and_denominators():
    c = pd.read_csv(RESULTS / "contrasts.csv").set_index("contrast")
    assert list(c.index) == ["P1", "S1", "S2", "S3", "D1", "D2", "D3", "D4", "D5", "D6", "D7"]
    assert c.loc["P1", "tier"] == "primary" and (c["tier"] == "primary").sum() == 1
    assert (c.loc["P1", "model_a"], c.loc["P1", "model_b"]) == ("ORIGINAL18", "ORIGINAL18_ESM1V")
    expected = {"plm_codingnt|all": (3121, 148), "plm_headline|all": (7125, 167)}
    for cid, row in c.iterrows():
        key = f"{row['task']}|{row['evaluation_set']}"
        if key in expected:
            assert (row["n_variants"], row["n_gene_groups"]) == expected[key], cid
        else:
            assert key == "plm_codingnt|scored_missense" and row["n_variants"] == 1328, cid
        assert row["n_valid_draws"] >= 1800 and row["status"] == "ok"
        assert row["delta_lo"] <= row["delta_hi"]
        assert row["lower_bound_above_zero"] == (row["delta_lo"] > 0)


def test_summary_verdict_follows_the_prespecified_rule():
    s = load(RESULTS / "summary.json")
    c = {r["contrast"]: r for r in s["contrasts"]}
    if not c["P1"]["lower_bound_above_zero"]:
        assert s["primary_verdict"] == "improvement_not_demonstrated"
    elif c["S3"]["lower_bound_above_zero"]:
        assert s["primary_verdict"] == "improvement_demonstrated_beyond_scored_status_marker"
    else:
        assert s["primary_verdict"] == "improvement_demonstrated_not_attributable_beyond_scored_status_marker"
    assert s["n_boot"] == 2000 and s["bootstrap_seed"] == 0
    assert s["n_scored_variants_headline"] == 1328


def test_model_metrics_cover_every_model_and_set():
    m = pd.read_csv(RESULTS / "model_metrics.csv")
    sets = {("plm_codingnt", "all"), ("plm_codingnt", "scored_missense"), ("plm_headline", "all")}
    assert set(zip(m["task"], m["evaluation_set"])) == sets
    for (task, subset), g in m.groupby(["task", "evaluation_set"]):
        assert set(g["model"]) == set(plm.PLM_FEATURE_SETS), (task, subset)
        assert g["n_variants"].nunique() == 1 and g["AUC_n_valid"].min() >= 1800
        assert ((g["AUC_lo"] <= g["AUC_estimate"]) & (g["AUC_estimate"] <= g["AUC_hi"])).all()


def test_report_quotes_the_primary_result():
    s = load(RESULTS / "summary.json")
    report = (ROOT / "docs/plm_results.md").read_text(encoding="utf-8")
    p = s["primary"]
    for value in (p["auc_a"], p["auc_b"]):
        assert f"{value:.3f}" in report
    for value in (p["delta"], p["delta_lo"], p["delta_hi"]):
        assert f"{value:+.3f}" in report
    assert s["primary_verdict"] in report
