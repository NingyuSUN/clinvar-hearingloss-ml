"""Tests for predict/live_sources.py and predict/predict_novel.py.

These make REAL network calls to Ensembl, gnomAD, UCSC and Hugging Face — no
mocks — because the whole point of this module is "does the live pipeline
actually reproduce the frozen cache's values." That makes them unsuitable for
required CI (external services can be slow, rate-limited, or briefly down for
reasons that have nothing to do with this repository), so they are skipped
unless RUN_NETWORK_TESTS=1 is set:

    RUN_NETWORK_TESTS=1 python -m pytest tests/test_predict_novel.py -v

They were passing, with every value below matching the frozen cache exactly
except ensembl_conservation (a deliberate, documented substitution — see
docs/predict_novel.md), when this file was written.
"""

import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "predict"))

pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_NETWORK_TESTS") != "1",
    reason="Set RUN_NETWORK_TESTS=1 to run tests that call live external services.",
)

pytest.importorskip("requests")
pytest.importorskip("polars")

import live_sources as ls  # noqa: E402

# A variant already in predict/bundle's frozen cache (16-2496587-G-C /
# TBC1D24), used as ground truth: every reproducible feature below should
# match the cache's recorded value exactly.
CHROM, POS, REF, ALT = "16", 2496587, "G", "C"
GENE = "TBC1D24"


def test_vep_picks_the_same_transcript_and_flags_as_the_frozen_cache():
    vep = ls.fetch_vep(CHROM, POS, REF, ALT)
    assert vep["variant_type"] == "missense"
    assert vep["gene_symbol"] == GENE
    assert vep["transcript_id"] == "ENST00000646147"
    assert vep["transcript_policy"] == "mane_select"
    flags = {k: v["value"] for k, v in vep["features"].items() if k.startswith("p1_")}
    assert flags == {
        "p1_has_frameshift": 0.0,
        "p1_has_stop_gained": 0.0,
        "p1_has_canonical_splice": 0.0,
        "p1_has_splice_region": 0.0,
    }
    assert vep["features"]["p3_in_specific_domain"]["value"] == 1.0
    assert vep["features"]["p3_specific_domain_count"]["value"] == 2.0
    assert vep["features"]["p3_variant_is_coding"]["value"] == 1.0


def test_vep_surfaces_a_useful_message_on_a_wrong_ref_allele():
    with pytest.raises(ls.AnnotationError) as exc_info:
        # chr1:1,000,000 reference base is G, not A — this should fail loudly.
        ls.fetch_vep("1", 1_000_000, "A", "G")
    assert "vep_http_400" in str(exc_info.value)


def test_gnomad_af_matches_the_frozen_cache():
    af = ls.fetch_gnomad_af(CHROM, POS, REF, ALT)
    assert af["p4_freq_log10"]["value"] == pytest.approx(-5.684291721999604, abs=1e-4)
    assert af["p2_af_known_missing"]["value"] == 0.0
    assert af["p2_not_observed_in_gnomad_r4"]["value"] == 0.0


def test_gnomad_constraint_matches_the_frozen_cache():
    constraint = ls.fetch_gnomad_constraint(GENE)
    assert constraint["p3_gene_constraint_oe_lof_upper"]["value"] == pytest.approx(1.123, abs=1e-3)
    assert constraint["p3_gene_constraint_oe_mis_upper"]["value"] == pytest.approx(0.974, abs=1e-3)


def test_gpn_star_matches_the_frozen_cache_exactly():
    gpn = ls.fetch_gpn_star(CHROM, POS, REF, ALT)
    assert gpn["gpn_v100"]["value"] == pytest.approx(-8.375, abs=1e-3)
    assert gpn["gpn_m447"]["value"] == pytest.approx(-6.839, abs=1e-3)
    assert gpn["gpn_p243"]["value"] == pytest.approx(-4.198, abs=1e-3)


def test_conservation_returns_a_finite_value_but_is_not_the_original_source():
    cons = ls.fetch_conservation(CHROM, POS)
    value = cons["ensembl_conservation"]["value"]
    assert value is not None
    # This is the UCSC phyloP substitute, not the original 1.71 GERP-style
    # value in the frozen cache — different source, different scale, on
    # purpose. This test only asserts it returns *something* usable.
    assert isinstance(value, float)


def test_end_to_end_novel_scoring_without_alphagenome_key(tmp_path):
    """predict_novel.py should score the three non-AVI models and cleanly
    report the two AVI-dependent models as missing when no key is given."""
    import subprocess

    input_csv = tmp_path / "variants.csv"
    input_csv.write_text(f"variant_id,assembly,chrom,pos,ref,alt\n48,GRCh38,{CHROM},{POS},{REF},{ALT}\n")
    output_dir = tmp_path / "out"

    env = dict(os.environ)
    env.pop("ALPHAGENOME_API_KEY", None)
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "predict" / "predict_novel.py"),
            "--input", str(input_csv),
            "--bundle", str(REPO_ROOT / "predict" / "bundle"),
            "--output", str(output_dir),
        ],
        capture_output=True, text=True, env=env, timeout=120,
    )
    assert result.returncode == 0, result.stderr

    import json

    record = json.loads((output_dir / "predictions.jsonl").read_text().splitlines()[0])
    assert record["variant_type"] == "missense"
    assert record["gene_symbols"] == [GENE]
    for name in ("ORIGINAL18", "GPN3", "ORIGINAL18_GPN"):
        assert record["models"][name]["status"] == "scored_research"
        assert record["models"][name]["probability"] is None
    for name in ("FULL_UNIFIED", "FULL_STRATIFIED"):
        assert record["models"][name]["status"] == "missing_or_invalid_features"
        assert record["models"][name]["score"] is None

    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["alphagenome_requested"] is False
    assert any("substitute" in limit for limit in manifest["limits"])
