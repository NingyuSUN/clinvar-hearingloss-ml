"""The frozen label / cohort / gene-split / feature protocol.

Every decision here is fixed before any model is trained: the label mapping, the
cohort gates, the connected-component gene grouping, the review-status tiers, and
the feature whitelist. The modelling matrix shipped in ``data/`` already carries
the resulting columns (``y``, ``p4_gene_group``, ``p4_review_star``, the cohort
flags, ``consequence_class``); :func:`build_matrix` reconstructs them from a raw
annotated feature table for anyone regenerating the inputs.
"""
from __future__ import annotations

import collections

import numpy as np
import pandas as pd

# --- label ------------------------------------------------------------------
# ClinVar aggregate ClinicalSignificance, not ClinSigSimple. "Likely" is kept.
PATHOGENIC = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
BENIGN = {"benign", "likely benign", "benign/likely benign"}


def canon_label(clinical_significance: object):
    """Aggregate P/LP -> 1, aggregate B/LB -> 0, everything else -> NA."""
    primary = str(clinical_significance).strip().lower().split(";")[0].strip()
    if primary in PATHOGENIC:
        return 1
    if primary in BENIGN:
        return 0
    return pd.NA


# --- review status (ClinVar gold stars) -----------------------------------
STAR = {
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}


def review_star(review_status: object) -> int:
    return STAR.get(str(review_status).strip().lower(), -1)


# --- gene grouping --------------------------------------------------------
def split_genes(gene_symbol: object) -> list[str]:
    out = []
    for part in str(gene_symbol).split(";"):
        part = part.strip()
        if part and part not in {"-", "nan"} and not part.startswith("covers "):
            out.append(part)
    return out


def build_gene_groups(all_symbols) -> dict[str, str]:
    """Connected components of gene co-occurrence: any two genes that ever appear
    together in a multi-gene annotation are held out together."""
    adj = collections.defaultdict(set)
    genes: set[str] = set()
    for symbol in all_symbols:
        parts = split_genes(symbol)
        genes.update(parts)
        for i in range(len(parts)):
            for j in range(i + 1, len(parts)):
                adj[parts[i]].add(parts[j])
                adj[parts[j]].add(parts[i])
    gid, seen, k = {}, set(), 0
    for g in sorted(genes):
        if g in seen:
            continue
        stack, comp = [g], set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.add(x)
            stack.extend(y for y in adj[x] if y not in seen)
        for x in comp:
            gid[x] = f"GG{k:03d}"
        k += 1
    return gid


def row_group(gene_symbol: object, gid: dict[str, str]) -> str:
    groups = sorted({gid[p] for p in split_genes(gene_symbol) if p in gid})
    if len(groups) == 1:
        return groups[0]
    if len(groups) > 1:
        return "|".join(groups)
    return "NA"


# --- features -----------------------------------------------------------
FEATURE_GROUPS = {
    "consequence": ["p1_has_frameshift", "p1_has_stop_gained", "p1_has_canonical_splice",
                    "p1_has_splice_region"],
    "allele": ["variant_length", "has_del", "has_ins", "has_dup"],
    "conservation": ["ensembl_conservation"],
    "constraint": ["p3_gene_constraint_oe_lof_upper", "p3_gene_constraint_oe_mis_upper"],
    "frequency": ["p4_freq_log10", "p2_af_known_missing", "p2_not_observed_in_gnomad_r4",
                  "p2_low_an_flag"],
    "domain": ["p3_in_specific_domain", "p3_specific_domain_count", "p3_variant_is_coding"],
}
ALL_FEATURES = [c for cols in FEATURE_GROUPS.values() for c in cols]
CONTINUOUS = {"variant_length", "ensembl_conservation", "p3_gene_constraint_oe_lof_upper",
              "p3_gene_constraint_oe_mis_upper", "p4_freq_log10", "p3_specific_domain_count"}


def _g(*names):
    return [c for n in names for c in FEATURE_GROUPS[n]]


FEATURE_SETS = {
    "base_freq_only": _g("frequency"),
    "base_consequence_only": _g("consequence"),
    "base_conservation_only": _g("conservation"),
    "L1_allele": _g("allele"),
    "L2_base": _g("allele", "conservation"),
    "L3_consequence": _g("allele", "conservation", "consequence"),
    "L4_constraint": _g("allele", "conservation", "consequence", "constraint"),
    "L5_frequency": _g("allele", "conservation", "consequence", "constraint", "frequency"),
    "L6_full": ALL_FEATURES,
    "F_minus_consequence": [c for c in ALL_FEATURES if c not in _g("consequence")],
    "F_minus_allele": [c for c in ALL_FEATURES if c not in _g("allele")],
    "F_minus_conservation": [c for c in ALL_FEATURES if c not in _g("conservation")],
    "F_minus_constraint": [c for c in ALL_FEATURES if c not in _g("constraint")],
    "F_minus_frequency": [c for c in ALL_FEATURES if c not in _g("frequency")],
    "F_minus_domain": [c for c in ALL_FEATURES if c not in _g("domain")],
}

# gnomAD v4.1.1: variants confirmed absent get a floor allele frequency ~ 1 / 2N_max
NOT_OBSERVED_FLOOR_AF = 1.0 / 1_600_000


def freq_log10(row) -> float:
    """log10(max PASS AF) for observed variants; a floor for confirmed-absent; NaN
    for filtered / unresolved (imputed per fold, paired with the two indicator flags)."""
    status = row.get("p2_af_status")
    af = row.get("p2_max_alt_af_pass")
    if status == "LIVE_PASS" and pd.notna(af):
        value = float(af)
    elif status == "NOT_IN_GNOMAD_R4":
        value = NOT_OBSERVED_FLOOR_AF
    else:
        return np.nan
    return float(np.log10(value + 1e-8))


def consequence_class(row) -> str:
    if row["p1_has_frameshift"] == 1 or row["p1_has_stop_gained"] == 1:
        return "truncating"
    if row["p1_has_canonical_splice"] == 1:
        return "canonical_splice"
    if row.get("p3_variant_is_coding", 0) == 1:
        return "coding_nontruncating"
    return "noncoding_or_other"


COHORT_GATE = [
    "y is not NA",
    "consequence annotation is reliable (single or consensus transcript)",
    "gnomAD gene-constraint present",
    "gnomAD frequency has a definitive status",
    "protein-domain call made (coding / non-coding decided)",
    "ensembl_conservation present (rows without it are dropped, not imputed)",
    "single-gene GeneSymbol",
]


def build_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct the modelling matrix from a raw annotated feature table.

    ``df`` must contain ``VariationID``, ``GeneSymbol``, ``ClinicalSignificance``,
    ``ReviewStatus``, the six feature groups, and the annotation-status columns
    (``p1_status``, ``p2_af_status`` / ``p2_max_alt_af_pass``, ``p3_variant_is_coding``).
    """
    gid = build_gene_groups(df["GeneSymbol"].astype(str).unique())
    out = pd.DataFrame({"VariationID": df["VariationID"]})
    out["y"] = df["ClinicalSignificance"].map(canon_label).astype("Int64")
    out["p4_gene_group"] = df["GeneSymbol"].map(lambda s: row_group(s, gid))
    out["p4_review_star"] = df["ReviewStatus"].map(review_star).astype("Int64")
    out["p4_freq_log10"] = df.apply(freq_log10, axis=1)
    for col in [c for c in ALL_FEATURES if c != "p4_freq_log10"]:
        out[col] = pd.to_numeric(df[col], errors="coerce")
    single_gene = ~df["GeneSymbol"].astype(str).str.contains(";")
    gate = (
        out["y"].notna()
        & df["p1_status"].isin(["ok_single", "ok_consensus"])
        & df["p3_gene_constraint_oe_lof_upper"].notna()
        & ~df["p2_af_status"].isin(["UNRESOLVED", "NO_VALID_KEY", "MT_NOT_SUPPORTED"])
        & df["p3_variant_is_coding"].notna()
        & df["ensembl_conservation"].notna()
        & single_gene
    )
    out["p4_in_modeling_cohort"] = gate.astype("Int64")           # all review tiers
    out["p4_headline_cohort"] = (gate & (out["p4_review_star"] >= 1)).astype("Int64")
    out["p4_sens_ge2star"] = (gate & (out["p4_review_star"] >= 2)).astype("Int64")
    out["p4_expert_panel"] = (gate & (out["p4_review_star"] >= 3)).astype("Int64")
    out["GeneSymbol"] = df["GeneSymbol"].astype(str)
    out["ReviewStatus"] = df["ReviewStatus"].astype(str)
    out["consequence_class"] = df.apply(consequence_class, axis=1)
    return out
