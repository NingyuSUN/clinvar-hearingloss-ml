"""Integrity of the committed, label-free PLM mapping and scoring inputs."""
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from hlpath import plm  # noqa: E402

PLM = ROOT / "data/plm"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name):
    return json.loads((PLM / name).read_text(encoding="utf-8"))


def test_mapping_outputs_and_inputs_match_manifest():
    manifest = load("mapping_manifest.json")
    assert manifest["protocol"] == plm.PROTOCOL_ID and manifest["label_column_read"] is False
    for name, digest in manifest["outputs_sha256"].items():
        assert sha(PLM / name) == digest, name
    for rel, digest in manifest["inputs_sha256"].items():
        assert sha(ROOT / rel) == digest, rel


def test_mapping_table_is_label_free_and_complete():
    mapping = pd.read_csv(PLM / "variant_protein_map.csv", low_memory=False)
    assert len(mapping) == 7125 and mapping["VariationID"].is_unique
    assert "y" not in mapping.columns and "ClinicalSignificance" not in mapping.columns
    counts = load("mapping_manifest.json")["counts"]
    assert mapping["map_status"].value_counts().to_dict() == counts["map_status"]
    scored = mapping[mapping["map_status"] == plm.STATUS_SCORED]
    assert set(scored["consequence_class"]) <= set(plm.CANDIDATE_CLASSES)
    noncandidates = mapping[~mapping["consequence_class"].isin(plm.CANDIDATE_CLASSES)]
    assert noncandidates["map_status"].eq(plm.STATUS_NOT_CANDIDATE).all()


def test_scoring_input_matches_mapping_and_window_rule():
    mapping = pd.read_csv(PLM / "variant_protein_map.csv", low_memory=False)
    spec = load("esm_scoring_input.json")
    proteins = load("proteins.json")
    assert spec["protocol"] == plm.PROTOCOL_ID and spec["max_residues"] == plm.MAX_RESIDUES
    scored = mapping[mapping["map_status"] == plm.STATUS_SCORED]
    want = {(plm.stable_id(t), int(p)) for t, p in zip(scored["cached_transcript_id"], scored["protein_position"])}
    have = set()
    for prot in spec["proteins"]:
        seq = prot["sequence"]
        assert seq == proteins[prot["transcript_stable_id"]]
        for pos in prot["positions"]:
            assert tuple(pos["window"]) == plm.esm_window(len(seq), pos["position"])
            have.add((prot["transcript_stable_id"], pos["position"]))
    assert have == want
    for _, row in scored.iterrows():
        seq = proteins[plm.stable_id(row["cached_transcript_id"])]
        assert seq[int(row["protein_position"]) - 1] == row["ref_aa"]
        assert len(seq) == row["protein_length"]
