"""Live annotation sources for a genuinely novel variant (not in predict/bundle's cache).

Every function here makes a real, unauthenticated-except-AlphaGenome network
call. This is deliberately a *separate, less-verified* path from
predict/prediction_core.py's hash-verified cache: it exists to let a
downloader score a variant this repository has never seen before, at the
cost of much weaker guarantees. Read docs/predict_novel.md before trusting
its output.

Each fetch function was validated against a real cached variant
(16-2496587-G-C / TBC1D24, "variant 48" in predict/examples) before this
module was written: VEP transcript selection, domain curation, gnomAD AF,
gnomAD gene constraint, and all three GPN-Star scores reproduced the frozen
cache's values exactly. The one deliberate exception is conservation — see
`fetch_conservation` below.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import requests

ENSEMBL_REST = "https://rest.ensembl.org"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
UCSC_API = "https://api.genome.ucsc.edu/getData/track"
GPN_STAR_ROOT = "hf://datasets/songlab/gpn-star-scores/data"

# Curated domain databases per docs/data.md — structure-model hits (PANTHER,
# Gene3D, SUPERFAMILY, ...) and non-functional tracks (ENSP_mappings, mobidb,
# phobius) are deliberately excluded.
CURATED_DOMAIN_DBS = {"Pfam", "SMART", "PROSITE_profiles", "PROSITE_patterns", "PRINTS", "CDD"}

# Consequence-class priority, matching docs/stratification.md's documented
# hierarchy: truncating before canonical splice, then the other classes.
# "truncating" and "noncoding_or_other" are not in predict/bundle's four
# supported TYPES (see prediction_core.TYPES) and will correctly come back
# as unsupported_variant_type.
CONSEQUENCE_PRIORITY = [
    ("truncating", {"frameshift_variant", "stop_gained"}),
    ("canonical_splice", {"splice_donor_variant", "splice_acceptor_variant"}),
    ("splice_region", {"splice_region_variant"}),
    ("missense", {"missense_variant"}),
    ("synonymous", {"synonymous_variant", "stop_retained_variant"}),
]


class AnnotationError(Exception):
    """A live source could not be resolved for this variant. Carries a short machine-readable reason."""


@dataclass
class NovelAnnotation:
    variant_type: str | None = None
    gene_symbols: list[str] = field(default_factory=list)
    transcript_id: str | None = None
    features: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def _obs(value):
    return {"value": value, "status": "observed"}


def fetch_vep(chrom: str, pos: int, ref: str, alt: str, timeout: float = 15.0) -> dict:
    """Query Ensembl VEP REST and pick one transcript: MANE Select > canonical > protein-coding > any.

    Note: the original frozen pipeline's transcript preference was
    "target-gene -> MANE Select -> canonical -> protein-coding -> any"
    (docs/data.md); the internal hearing-loss gene-panel list used for the
    "target-gene" tier was not preserved anywhere this repo could find, so
    this live path only reproduces the last three tiers. This matched the
    frozen cache exactly for every case checked so far, but is a documented,
    deliberate simplification, not a guarantee.
    """
    region = f"{chrom}:{pos}-{pos}:1"
    resp = requests.get(
        f"{ENSEMBL_REST}/vep/human/region/{region}/{alt}",
        params={"content-type": "application/json", "domains": 1, "canonical": 1, "mane": 1},
        headers={"Content-Type": "application/json"},
        timeout=timeout,
    )
    if resp.status_code != 200:
        try:
            detail = resp.json().get("error", "")
        except ValueError:
            detail = resp.text[:200]
        # A very common cause here: the supplied `ref` does not match the
        # actual GRCh38 reference base at this position (VEP reports this as
        # "request for consequence of [ALT] matches reference [ALT]" when the
        # caller's ref/alt are swapped or simply wrong). Surface VEP's own
        # message rather than just the status code.
        raise AnnotationError(f"vep_http_{resp.status_code}:{detail}")
    payload = resp.json()
    if not payload or "transcript_consequences" not in payload[0]:
        raise AnnotationError("vep_no_transcript_consequences")

    transcripts = payload[0]["transcript_consequences"]

    def rank(tc):
        return (
            0 if tc.get("mane_select") else 1,
            0 if tc.get("canonical") else 1,
            0 if tc.get("biotype") == "protein_coding" else 1,
        )

    chosen = min(transcripts, key=rank)
    terms = set(chosen.get("consequence_terms", []))
    domains = chosen.get("domains") or []
    curated = [d for d in domains if d.get("db") in CURATED_DOMAIN_DBS]
    is_coding = chosen.get("biotype") == "protein_coding" and bool(chosen.get("protein_start"))

    variant_type = None
    for name, term_set in CONSEQUENCE_PRIORITY:
        if terms & term_set:
            variant_type = name
            break

    return {
        "variant_type": variant_type,
        "consequence_terms": sorted(terms),
        "gene_symbol": chosen.get("gene_symbol"),
        "transcript_id": chosen.get("transcript_id"),
        "transcript_policy": "mane_select" if chosen.get("mane_select") else ("canonical" if chosen.get("canonical") else "protein_coding_or_any"),
        "features": {
            "p1_has_frameshift": _obs(1.0 if "frameshift_variant" in terms else 0.0),
            "p1_has_stop_gained": _obs(1.0 if "stop_gained" in terms else 0.0),
            "p1_has_canonical_splice": _obs(1.0 if terms & {"splice_donor_variant", "splice_acceptor_variant"} else 0.0),
            "p1_has_splice_region": _obs(1.0 if "splice_region_variant" in terms else 0.0),
            "p3_in_specific_domain": _obs(1.0 if curated else 0.0),
            "p3_specific_domain_count": _obs(float(len(curated))),
            "p3_variant_is_coding": _obs(1.0 if is_coding else 0.0),
        },
    }


def fetch_gnomad_af(chrom: str, pos: int, ref: str, alt: str, timeout: float = 15.0) -> dict:
    """gnomAD v4.1.1 max(genome, exome) PASS allele frequency, matching docs/data.md's p4_freq_log10 rule."""
    query = """
    query VariantInfo($variantId: String!, $dataset: DatasetId!) {
      variant(variantId: $variantId, dataset: $dataset) {
        genome { af filters an }
        exome { af filters an }
      }
    }
    """
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    resp = requests.post(
        GNOMAD_API,
        json={"query": query, "variables": {"variantId": variant_id, "dataset": "gnomad_r4"}},
        timeout=timeout,
    )
    if resp.status_code != 200:
        raise AnnotationError(f"gnomad_http_{resp.status_code}")
    data = resp.json().get("data", {}).get("variant")

    if data is None:
        # Confirmed absent from gnomAD v4 (docs/data.md's NOT_IN_GNOMAD_R4 case).
        floor = math.log10(1 / (2 * 1_614_000) + 1e-8)  # 2*N_max is dataset-scale; this floor is an approximation, not the frozen exact value
        return {
            "p4_freq_log10": _obs(floor),
            "p2_af_known_missing": _obs(0.0),
            "p2_not_observed_in_gnomad_r4": _obs(1.0),
            "p2_low_an_flag": _obs(0.0),
        }

    candidates = [(part["af"], part.get("an", 0)) for part in (data.get("genome"), data.get("exome")) if part and not part.get("filters")]
    if not candidates:
        # Present but filtered/non-PASS everywhere this query looked.
        return {
            "p4_freq_log10": {"value": None, "status": "frozen_missing"},
            "p2_af_known_missing": _obs(1.0),
            "p2_not_observed_in_gnomad_r4": _obs(0.0),
            "p2_low_an_flag": _obs(0.0),
        }

    af, an = max(candidates, key=lambda pair: pair[0])
    return {
        "p4_freq_log10": _obs(math.log10(af + 1e-8)),
        "p2_af_known_missing": _obs(0.0),
        "p2_not_observed_in_gnomad_r4": _obs(0.0),
        "p2_low_an_flag": _obs(1.0 if an < 2000 else 0.0),
    }


def fetch_gnomad_constraint(gene_symbol: str, timeout: float = 15.0) -> dict:
    """gnomAD v2.1.1 gene constraint (oe_lof_upper / oe_mis_upper), queried live via gnomAD's own API."""
    query = """
    query GeneConstraint($gene: String!) {
      gene(gene_symbol: $gene, reference_genome: GRCh37) {
        gnomad_constraint { oe_lof_upper oe_mis_upper }
      }
    }
    """
    resp = requests.post(GNOMAD_API, json={"query": query, "variables": {"gene": gene_symbol}}, timeout=timeout)
    if resp.status_code != 200:
        raise AnnotationError(f"gnomad_constraint_http_{resp.status_code}")
    gene = resp.json().get("data", {}).get("gene")
    constraint = (gene or {}).get("gnomad_constraint")
    if not constraint or constraint.get("oe_lof_upper") is None or constraint.get("oe_mis_upper") is None:
        return {
            "p3_gene_constraint_oe_lof_upper": {"value": None, "status": "missing"},
            "p3_gene_constraint_oe_mis_upper": {"value": None, "status": "missing"},
        }
    return {
        "p3_gene_constraint_oe_lof_upper": _obs(float(constraint["oe_lof_upper"])),
        "p3_gene_constraint_oe_mis_upper": _obs(float(constraint["oe_mis_upper"])),
    }


def fetch_conservation(chrom: str, pos: int, timeout: float = 15.0) -> dict:
    """UCSC phyloP (100-way vertebrate) at this position, via UCSC's public REST API.

    This is a DELIBERATE SUBSTITUTE for the frozen cache's `ensembl_conservation`
    feature. The original feature's provenance was never fully documented in
    this repository (docs/limitations.md), and Ensembl's own GERP conservation
    score has no REST endpoint — only the Compara Perl API or raw packed-binary
    MySQL access, both impractical to ask a downloader to set up. phyloP and
    GERP are related but not the same score, on different scales; do not treat
    this feature as reproducing the training data's original values.
    """
    resp = requests.get(
        UCSC_API,
        params={"genome": "hg38", "track": "phyloP100way", "chrom": f"chr{chrom}", "start": pos - 1, "end": pos},
        timeout=timeout,
    )
    if resp.status_code != 200:
        raise AnnotationError(f"ucsc_http_{resp.status_code}")
    rows = resp.json().get("phyloP100way") or []
    if not rows:
        return {"ensembl_conservation": {"value": None, "status": "missing"}}
    return {"ensembl_conservation": _obs(float(rows[0]["value"]))}


def fetch_gpn_star(chrom: str, pos: int, ref: str, alt: str, timeout: float = 30.0) -> dict:
    """The three GPN-Star LLR scores, via a remote Parquet row lookup (no local model, no GPU).

    Requires `polars`. Uses `songlab/gpn-star-scores` on Hugging Face — the
    exact same public precomputed dataset the frozen cache's gpn_v100/m447/p243
    values came from.
    """
    import polars as pl

    checkpoints = {"gpn_v100": "gpn-star-hg38-v100-200m", "gpn_m447": "gpn-star-hg38-m447-200m", "gpn_p243": "gpn-star-hg38-p243-200m"}
    out = {}
    for feature_name, dataset in checkpoints.items():
        path = f"{GPN_STAR_ROOT}/{dataset}/llr/llr_chr{chrom}.parquet"
        try:
            row = (
                pl.scan_parquet(path)
                .filter((pl.col("pos") == pos) & (pl.col("ref") == ref) & (pl.col("alt") == alt))
                .collect()
            )
        except Exception as exc:  # noqa: BLE001 — surface as a missing feature, not a crash
            out[feature_name] = {"value": None, "status": "missing", "error": str(exc)[:200]}
            continue
        if row.is_empty():
            out[feature_name] = {"value": None, "status": "missing"}
        else:
            out[feature_name] = _obs(float(row["llr_calibrated"][0]))
    return out


def fetch_alphagenome(chrom: str, pos: int, ref: str, alt: str, gene_symbol: str, api_key: str, timeout: float = 30.0) -> dict:
    """AVI_SCORE and SPLICE_SITES from the AlphaGenome API, using the caller's own key.

    Reproduces the exact scorer names and aggregation used to build the
    frozen cache's avi/splice_sites (see runs/alphagenome_targeted_2026-09-14/
    in the project's working notes): AVI_SCORE's single (1,1) value is `avi`;
    SPLICE_SITES' donor/acceptor scores are max-aggregated over rows whose
    gene_name matches `gene_symbol`.

    AVI is a composite score that itself blends AlphaGenome, AlphaMissense
    and related outputs — it is NOT a pure AlphaGenome score. Requires the
    `alphagenome` package and your own API key (apply at
    https://deepmind.google.com/science/alphagenome).
    """
    from alphagenome.atlas import atlas
    from alphagenome.data import genome

    client = atlas.create(api_key, timeout=timeout)
    variant = genome.Variant(chromosome=f"chr{chrom}", position=pos, reference_bases=ref, alternate_bases=alt)
    result = client.query_variant(variant, requested_scorers=["AVI_SCORE", "SPLICE_SITES"])

    out = {"avi": {"value": None, "status": "missing"}, "splice_sites": {"value": None, "status": "missing"}}
    for name, array in result.items():
        values = array.X
        obs_rows = array.obs.reset_index().to_dict("records")
        track_rows = array.var.reset_index().to_dict("records")
        if name == "AVI_SCORE":
            if values.shape == (1, 1):
                out["avi"] = _obs(float(values[0, 0]))
        elif name == "SPLICE_SITES":
            matching = []
            for i, obs_row in enumerate(obs_rows):
                if obs_row.get("gene_name") != gene_symbol:
                    continue
                for j in range(len(track_rows)):
                    v = float(values[i, j])
                    if math.isfinite(v):
                        matching.append(v)
            if matching:
                out["splice_sites"] = _obs(max(matching))
    return out
