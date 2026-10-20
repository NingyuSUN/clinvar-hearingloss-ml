"""Protein language model (PLM) feature rules for the PLM-2026-10 protocol.

Everything here is label-blind and deterministic. ``docs/plm_protocol.md`` is
the prespecified protocol; this module is its executable form:

  * which variants receive an ESM score (missense on the frozen cached
    transcript, reference residue verified against the Ensembl translation);
  * the residue window shown to a 1,022-residue ESM model;
  * the masked-marginal log-likelihood ratio and the ESM-1v ensemble mean;
  * the feature sets compared (ORIGINAL18 is the existing 18-feature model).

Unscored variants (synonymous, in-frame indels, start/stop changes, failed
mapping) keep ``NaN`` scores and are imputed with the inner-train median inside
every split, exactly like the existing continuous features.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .protocol import ALL_FEATURES, CONTINUOUS

PROTOCOL_ID = "PLM-ESM-2026-10"

AA20 = "ACDEFGHIKLMNPQRSTVWY"

# Frozen matrix classes whose variants are sent to VEP. Truncating and
# non-coding classes cannot carry a missense change on the cached transcript.
CANDIDATE_CLASSES = ("coding_nontruncating", "canonical_splice")

# ESM-1v / ESM-2 positional limit: 1,024 tokens including <cls> and <eos>.
MAX_RESIDUES = 1022

ESM1V_MODELS = tuple(f"facebook/esm1v_t33_650M_UR90S_{i}" for i in range(1, 6))
ESM2_MODEL = "facebook/esm2_t33_650M_UR50D"
# Hugging Face commit hashes pinned before scoring (queried 2026-10-19).
MODEL_REVISIONS = {
    "facebook/esm1v_t33_650M_UR90S_1": "8bfdb1892536cc77bd0760b9c25ddced2cd0b4c8",
    "facebook/esm1v_t33_650M_UR90S_2": "3c1e9e64480f069b163e456c361bfd80a8bab04c",
    "facebook/esm1v_t33_650M_UR90S_3": "0b00fd112e63f6b5e70a9cd8484d4e660312ce70",
    "facebook/esm1v_t33_650M_UR90S_4": "443968f644da132d323bbc6321a6d443149fa57b",
    "facebook/esm1v_t33_650M_UR90S_5": "fb2e51cb0f605cbad2c4ca3bb784be7bc8e2f4a8",
    "facebook/esm2_t33_650M_UR50D": "08e4846e537177426273712802403f7ba8261b6c",
}

SCORE_COLUMNS = ("esm1v_llr", "esm2_llr")
SCORED_FLAG = "plm_scored"

PLM_FEATURE_SETS = {
    "ORIGINAL18": list(ALL_FEATURES),
    "ESM1V_ONLY": ["esm1v_llr"],
    "ORIGINAL18_ESM1V": [*ALL_FEATURES, "esm1v_llr"],
    "ORIGINAL18_SCOREDFLAG": [*ALL_FEATURES, SCORED_FLAG],
    "ESM2_ONLY": ["esm2_llr"],
    "ORIGINAL18_ESM2": [*ALL_FEATURES, "esm2_llr"],
}
PLM_CONTINUOUS = set(CONTINUOUS) | set(SCORE_COLUMNS)

# Mapping statuses. Only "scored" receives an ESM score.
STATUS_SCORED = "scored"
STATUS_NOT_CANDIDATE = "not_candidate_class"
STATUS_NO_TRANSCRIPT = "no_single_cached_transcript"
STATUS_VEP_NO_RESULT = "vep_no_result"
STATUS_TRANSCRIPT_ABSENT = "cached_transcript_not_in_vep"
STATUS_NOT_MISSENSE = "not_missense"
STATUS_MULTI_RESIDUE = "not_single_residue_substitution"
STATUS_NONSTANDARD = "nonstandard_residue"
STATUS_NO_SEQUENCE = "no_protein_sequence"
STATUS_REF_MISMATCH = "reference_residue_mismatch"


def stable_id(transcript_id: str) -> str:
    """``ENST00000646147.1`` -> ``ENST00000646147``."""
    return str(transcript_id).split(".", 1)[0]


def vcf_input(variant_key: str) -> str:
    """``16-2496587-G-C`` -> the VEP REST VCF-style region string."""
    chrom, pos, ref, alt = str(variant_key).split("-")
    if not pos.isdigit() or not ref or not alt:
        raise ValueError(f"malformed variant key: {variant_key!r}")
    return f"{chrom} {pos} . {ref} {alt} . . ."


def select_transcript_consequence(vep_record: dict | None, transcript_id: str) -> dict | None:
    """Return the VEP transcript consequence for the cached transcript's stable ID."""
    if not vep_record:
        return None
    target = stable_id(transcript_id)
    hits = [tc for tc in vep_record.get("transcript_consequences", [])
            if stable_id(tc.get("transcript_id", "")) == target]
    if len(hits) > 1:
        raise ValueError(f"VEP returned {len(hits)} entries for {target}")
    return hits[0] if hits else None


def map_consequence(entry: dict | None, sequence: str | None) -> dict:
    """Apply the protocol's scoring rule to one VEP transcript consequence."""
    if entry is None:
        return {"map_status": STATUS_TRANSCRIPT_ABSENT}
    out = {
        "vep_transcript_id": entry.get("transcript_id"),
        "protein_id": entry.get("protein_id"),
        "consequence_terms": "&".join(entry.get("consequence_terms", [])),
    }
    if "missense_variant" not in entry.get("consequence_terms", []):
        return out | {"map_status": STATUS_NOT_MISSENSE}
    aa = str(entry.get("amino_acids", ""))
    start, end = entry.get("protein_start"), entry.get("protein_end")
    parts = aa.split("/")
    if len(parts) != 2 or len(parts[0]) != 1 or len(parts[1]) != 1 or start is None or start != end:
        return out | {"map_status": STATUS_MULTI_RESIDUE, "amino_acids": aa}
    ref, alt = parts
    if ref not in AA20 or alt not in AA20 or ref == alt:
        return out | {"map_status": STATUS_NONSTANDARD, "amino_acids": aa}
    out |= {"protein_position": int(start), "ref_aa": ref, "alt_aa": alt}
    if not sequence:
        return out | {"map_status": STATUS_NO_SEQUENCE}
    out["protein_length"] = len(sequence)
    if int(start) > len(sequence) or sequence[int(start) - 1] != ref:
        return out | {"map_status": STATUS_REF_MISMATCH}
    return out | {"map_status": STATUS_SCORED}


def esm_window(length: int, position: int, max_residues: int = MAX_RESIDUES) -> tuple[int, int]:
    """Residue window ``[start, end)`` (0-based) shown to the model for a
    1-based ``position``. Whole protein when it fits; otherwise a
    ``max_residues`` window centred on the position, clamped to the ends."""
    if length < 1 or not 1 <= position <= length:
        raise ValueError("position must lie inside the protein")
    if length <= max_residues:
        return 0, length
    start = min(max(position - 1 - max_residues // 2, 0), length - max_residues)
    return start, start + max_residues


def masked_marginal_llr(logprobs: dict, ref: str, alt: str) -> float:
    """log p(alt) - log p(ref) at the masked position. Negative = less likely
    than wild type, the direction ESM-1v associates with damaging variants."""
    lp_ref, lp_alt = float(logprobs[ref]), float(logprobs[alt])
    if not (math.isfinite(lp_ref) and math.isfinite(lp_alt)):
        raise ValueError("log-probabilities must be finite")
    return lp_alt - lp_ref


def variant_llrs(mapping: pd.DataFrame, logprobs: pd.DataFrame) -> pd.DataFrame:
    """Combine the mapping table with per-position log-probabilities.

    ``logprobs`` has one row per (transcript stable ID, protein_position, model)
    with columns ``lp_<AA>`` for the 20 standard residues. Returns one row per
    VariationID with ``esm1v_llr`` (mean of the five ESM-1v LLRs; NaN unless all
    five are present), ``esm2_llr``, ``plm_scored`` and the per-model LLRs.
    """
    lp_cols = [f"lp_{a}" for a in AA20]
    need = {"transcript_stable_id", "protein_position", "model", *lp_cols}
    if missing := need - set(logprobs.columns):
        raise ValueError(f"log-probability table lacks {sorted(missing)}")
    if logprobs.duplicated(["transcript_stable_id", "protein_position", "model"]).any():
        raise ValueError("duplicate (transcript, position, model) log-probabilities")
    if not np.isfinite(logprobs[lp_cols].to_numpy(float)).all():
        raise ValueError("log-probabilities must be finite")
    if mapping["VariationID"].duplicated().any():
        raise ValueError("duplicate VariationID in mapping")

    scored = mapping[mapping["map_status"] == STATUS_SCORED].copy()
    scored["transcript_stable_id"] = scored["cached_transcript_id"].map(stable_id)
    scored["protein_position"] = scored["protein_position"].astype(int)
    # Several variants can share a position (different alternate residues).
    merged = scored.merge(logprobs, on=["transcript_stable_id", "protein_position"], how="left")
    merged["llr"] = [
        masked_marginal_llr({a: r[f"lp_{a}"] for a in AA20}, r["ref_aa"], r["alt_aa"])
        if pd.notna(r["model"]) else np.nan
        for _, r in merged.iterrows()
    ]
    wide = merged.pivot_table(index="VariationID", columns="model", values="llr", aggfunc="first")
    out = pd.DataFrame({"VariationID": mapping["VariationID"].to_numpy()})
    per_model = {}
    for model in (*ESM1V_MODELS, ESM2_MODEL):
        col = f"llr__{model.split('/')[-1]}"
        per_model[col] = out["VariationID"].map(wide[model]) if model in wide else np.nan
    out = out.assign(**per_model)
    esm1v_cols = [f"llr__{m.split('/')[-1]}" for m in ESM1V_MODELS]
    complete = out[esm1v_cols].notna().all(axis=1)
    out["esm1v_llr"] = out[esm1v_cols].mean(axis=1).where(complete)
    out["esm2_llr"] = out[f"llr__{ESM2_MODEL.split('/')[-1]}"]
    out[SCORED_FLAG] = complete.astype(int)
    return out


def attach_plm_features(matrix: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """Left-join PLM features onto the modelling matrix by VariationID; variants
    absent from ``features`` are unscored (NaN scores, ``plm_scored`` = 0)."""
    cols = ["VariationID", *SCORE_COLUMNS, SCORED_FLAG]
    if features["VariationID"].duplicated().any():
        raise ValueError("duplicate VariationID in PLM features")
    if overlap := (set(cols) - {"VariationID"}) & set(matrix.columns):
        raise ValueError(f"matrix already has PLM columns: {sorted(overlap)}")
    out = matrix.merge(features[cols], on="VariationID", how="left", validate="one_to_one")
    out[SCORED_FLAG] = out[SCORED_FLAG].fillna(0).astype(int)
    if (out.loc[out[SCORED_FLAG] == 1, "esm1v_llr"].isna()).any():
        raise ValueError("scored variants must have an ESM-1v score")
    return out
