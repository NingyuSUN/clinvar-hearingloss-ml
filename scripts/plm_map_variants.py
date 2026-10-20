#!/usr/bin/env python3
"""Map headline-cohort variants to protein substitutions for ESM scoring.

    python scripts/plm_map_variants.py --out data/plm

Label-blind: the label column is never read. Uses the frozen cached transcript
(``predict/bundle/annotation_cache.jsonl``), the Ensembl VEP REST API and the
Ensembl translation of that transcript. Writes

  * ``variant_protein_map.csv``   one row per headline variant with a map status
  * ``vep_cached_transcript.jsonl.gz``  the VEP entry for the cached transcript
  * ``proteins.json``             transcript stable ID -> translation
  * ``esm_scoring_input.json``    the (protein, position) pairs to score
  * ``mapping_manifest.json``     Ensembl release, input/output SHA256, counts
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

from hlpath.config import REPO_ROOT
from hlpath.plm import (CANDIDATE_CLASSES, MAX_RESIDUES, PROTOCOL_ID, STATUS_NO_TRANSCRIPT,
                        STATUS_NOT_CANDIDATE, STATUS_SCORED, STATUS_VEP_NO_RESULT, esm_window,
                        map_consequence, select_transcript_consequence, stable_id, vcf_input)

REST = "https://rest.ensembl.org"
VEP_BATCH = 200
SEQ_BATCH = 50
VEP_PARAMS = {"protein": 1, "mane": 1, "canonical": 1, "transcript_version": 1}
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
MATRIX_COLUMNS = ["VariationID", "p4_headline_cohort", "consequence_class", "GeneSymbol",
                  "p4_gene_group"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def post(session: requests.Session, path: str, payload: dict, params: dict) -> list:
    for attempt in range(6):
        resp = session.post(f"{REST}{path}", params=params, json=payload, headers=HEADERS,
                            timeout=300)
        if resp.status_code == 200:
            return resp.json()
        if resp.status_code in (429, 500, 502, 503, 504):
            time.sleep(float(resp.headers.get("Retry-After", 2 ** attempt)))
            continue
        raise RuntimeError(f"{path}: HTTP {resp.status_code}: {resp.text[:300]}")
    raise RuntimeError(f"{path}: retries exhausted")


def get_json(session: requests.Session, path: str) -> dict:
    resp = session.get(f"{REST}{path}", headers=HEADERS, timeout=60)
    resp.raise_for_status()
    return resp.json()


def write_gz_jsonl(path: Path, rows: list[dict]) -> None:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        for row in rows:
            gz.write((json.dumps(row, sort_keys=True) + "\n").encode())
    path.write_bytes(buf.getvalue())


def load_inputs(root: Path) -> tuple[pd.DataFrame, dict[Path, str]]:
    paths = {
        "matrix": root / "data/modeling_matrix.csv.gz",
        "keys": root / "predict/examples/known_variants_7125.csv",
        "cache": root / "predict/bundle/annotation_cache.jsonl",
    }
    matrix = pd.read_csv(paths["matrix"], usecols=MATRIX_COLUMNS, low_memory=False)
    headline = matrix[matrix["p4_headline_cohort"] == 1].copy()
    keys = pd.read_csv(paths["keys"], dtype=str)
    keys["VariationID"] = keys["variant_id"].astype(int)
    keys["variant_key"] = keys["chrom"] + "-" + keys["pos"] + "-" + keys["ref"] + "-" + keys["alt"]
    cache = pd.DataFrame([
        {"variant_key": r["variant_key"], "cached_transcript_id": r["transcript_ids"],
         "cached_transcript_policy": r["transcript_policy"],
         "cached_variant_type": r["variant_type"]}
        for r in map(json.loads, paths["cache"].read_text().splitlines())
    ])
    df = (headline.merge(keys[["VariationID", "variant_key"]], on="VariationID", how="left",
                         validate="one_to_one")
          .merge(cache, on="variant_key", how="left", validate="many_to_one"))
    if df["variant_key"].isna().any() or df["cached_transcript_id"].isna().any():
        raise SystemExit("every headline variant needs a variant key and a cached transcript")
    if len(df) != 7125:
        raise SystemExit(f"expected 7,125 headline variants, found {len(df)}")
    return df, {p: sha256(p) for p in paths.values()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "data/plm")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    df, input_hashes = load_inputs(REPO_ROOT)
    session = requests.Session()
    release = get_json(session, "/info/data")["releases"]
    software = get_json(session, "/info/software")

    single = ~df["cached_transcript_id"].str.contains(r"[,;|]")
    candidate = df["consequence_class"].isin(CANDIDATE_CLASSES)
    query = df[candidate & single]
    inputs = sorted({vcf_input(k) for k in query["variant_key"]})
    vep: dict[str, dict] = {}
    for i in range(0, len(inputs), VEP_BATCH):
        for rec in post(session, "/vep/human/region", {"variants": inputs[i:i + VEP_BATCH]},
                        VEP_PARAMS):
            if rec["input"] in vep:
                raise SystemExit(f"duplicate VEP record for {rec['input']}")
            vep[rec["input"]] = rec
        time.sleep(0.2)

    entries, vep_rows = {}, []
    for _, row in query.iterrows():
        rec = vep.get(vcf_input(row["variant_key"]))
        entry = select_transcript_consequence(rec, row["cached_transcript_id"])
        entries[row["VariationID"]] = (rec, entry)
        vep_rows.append({"VariationID": int(row["VariationID"]), "input": vcf_input(row["variant_key"]),
                         "cached_transcript_id": row["cached_transcript_id"],
                         "vep_found": rec is not None, "entry": entry})

    need_seq = sorted({stable_id(row["cached_transcript_id"]) for _, row in query.iterrows()
                       if (e := entries[row["VariationID"]][1]) is not None
                       and "missense_variant" in e.get("consequence_terms", [])})
    proteins: dict[str, str] = {}
    for i in range(0, len(need_seq), SEQ_BATCH):
        for rec in post(session, "/sequence/id", {"ids": need_seq[i:i + SEQ_BATCH]},
                        {"type": "protein", "object_type": "transcript"}):
            proteins[stable_id(rec["query"])] = rec["seq"]
        time.sleep(0.2)

    rows = []
    for _, row in df.iterrows():
        base = {k: row[k] for k in ["VariationID", "variant_key", "GeneSymbol", "p4_gene_group",
                                    "consequence_class", "cached_transcript_id",
                                    "cached_transcript_policy", "cached_variant_type"]}
        if row["consequence_class"] not in CANDIDATE_CLASSES:
            rows.append(base | {"map_status": STATUS_NOT_CANDIDATE})
            continue
        if not single[row.name]:
            rows.append(base | {"map_status": STATUS_NO_TRANSCRIPT})
            continue
        rec, entry = entries[row["VariationID"]]
        if rec is None:
            rows.append(base | {"map_status": STATUS_VEP_NO_RESULT})
            continue
        seq = proteins.get(stable_id(row["cached_transcript_id"]))
        rows.append(base | map_consequence(entry, seq))
    mapping = pd.DataFrame(rows).sort_values("VariationID").reset_index(drop=True)
    for col in ["protein_position", "protein_length"]:
        mapping[col] = mapping[col].astype("Int64")

    scored = mapping[mapping["map_status"] == STATUS_SCORED]
    positions = (scored.assign(transcript_stable_id=scored["cached_transcript_id"].map(stable_id))
                 .groupby("transcript_stable_id")["protein_position"]
                 .agg(lambda s: sorted({int(v) for v in s})))
    scoring_input = {
        "protocol": PROTOCOL_ID, "ensembl_release": release, "max_residues": MAX_RESIDUES,
        "proteins": [
            {"transcript_stable_id": tid, "sequence": proteins[tid], "positions": [
                {"position": p, "window": list(esm_window(len(proteins[tid]), p))} for p in pos]}
            for tid, pos in positions.items()
        ],
    }

    paths = {
        "map": args.out / "variant_protein_map.csv",
        "vep": args.out / "vep_cached_transcript.jsonl.gz",
        "proteins": args.out / "proteins.json",
        "scoring_input": args.out / "esm_scoring_input.json",
    }
    mapping.to_csv(paths["map"], index=False, lineterminator="\n")
    write_gz_jsonl(paths["vep"], sorted(vep_rows, key=lambda r: r["VariationID"]))
    paths["proteins"].write_text(json.dumps(dict(sorted(proteins.items())), indent=1) + "\n")
    paths["scoring_input"].write_text(json.dumps(scoring_input, indent=1) + "\n")

    status = mapping["map_status"].value_counts().to_dict()
    manifest = {
        "protocol": PROTOCOL_ID,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "label_column_read": False,
        "ensembl_rest": REST, "ensembl_release": release, "ensembl_software": software,
        "vep_endpoint": "POST /vep/human/region", "vep_params": VEP_PARAMS,
        "sequence_endpoint": "POST /sequence/id (type=protein, object_type=transcript)",
        "inputs_sha256": {str(p.relative_to(REPO_ROOT)).replace("\\", "/"): h
                          for p, h in input_hashes.items()},
        "outputs_sha256": {p.name: sha256(p) for p in paths.values()},
        "counts": {
            "headline_variants": int(len(mapping)),
            "vep_queried_variants": int(len(query)),
            "map_status": {k: int(v) for k, v in sorted(status.items())},
            "scored_by_consequence_class": {k: int(v) for k, v in
                                            scored["consequence_class"].value_counts().items()},
            "scored_by_cached_variant_type": {k: int(v) for k, v in
                                              scored["cached_variant_type"].value_counts().items()},
            "vep_transcript_version_changed": int(
                (scored["vep_transcript_id"] != scored["cached_transcript_id"]).sum()),
            "proteins_scored": int(len(positions)),
            "positions_scored": int(sum(len(v) for v in positions)),
            "proteins_longer_than_window": int(sum(len(proteins[t]) > MAX_RESIDUES
                                                   for t in positions.index)),
            "status_by_cached_variant_type": {
                f"{a}|{b}": int(n) for (a, b), n in
                Counter(zip(mapping["cached_variant_type"], mapping["map_status"])).items()},
        },
    }
    (args.out / "mapping_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
