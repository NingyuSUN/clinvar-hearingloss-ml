"""Live annotation sources for a genuinely novel variant (not in predict/bundle's cache).

Every function here makes a real, unauthenticated-except-AlphaGenome network
call. This is deliberately a *separate, less-verified* path from
predict/prediction_core.py's hash-verified cache: it exists to let a
downloader score a variant this repository has never seen before, at the
cost of much weaker guarantees. Read docs/predict_novel.md before trusting
its output.

Historical single-variant spot checks are not external validation. The live
conservation track is reference-only and cannot fill the training feature.
"""

from __future__ import annotations

import math

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


# Frozen transform from src/hlpath/protocol.py; checked against the bundle.
NOT_OBSERVED_FLOOR_AF = 1.0 / 1_600_000


def _obs(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise AnnotationError("invalid_numeric")
    return {"value": value, "status": "observed"}


def _request_json(method, source, url, **kwargs):
    try:
        response = getattr(requests, method)(url, **kwargs)
        if response.status_code != 200:
            raise AnnotationError(f"{source}_http_{response.status_code}")
        return response.json()
    except requests.RequestException as exc:
        raise AnnotationError(f"{source}_transport_error") from exc
    except ValueError as exc:
        raise AnnotationError(f"{source}_invalid_json") from exc


def _graphql(payload, key):
    # HTTP 200 may still carry execution errors, including partial data.
    if not isinstance(payload, dict) or ("errors" in payload):
        raise AnnotationError("gnomad_graphql_error")
    data = payload.get("data")
    if not isinstance(data, dict) or key not in data:
        raise AnnotationError("gnomad_invalid_response")
    value = data[key]
    if value is not None and not isinstance(value, dict):
        raise AnnotationError("gnomad_invalid_response")
    return value


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
    payload = _request_json("get", "vep",
        f"{ENSEMBL_REST}/vep/human/region/{region}/{alt}",
        params={"content-type": "application/json", "domains": 1, "canonical": 1, "mane": 1},
        headers={"Content-Type": "application/json"},
        timeout=timeout,
    )
    if not isinstance(payload, list) or len(payload) != 1 or not isinstance(payload[0], dict):
        raise AnnotationError("vep_invalid_response")
    record = payload[0]
    expected = {"assembly_name": "GRCh38", "seq_region_name": chrom,
                "start": pos, "end": pos, "strand": 1, "allele_string": f"{ref}/{alt}"}
    if any(record.get(key) != value for key, value in expected.items()):
        raise AnnotationError("vep_variant_identity_mismatch")
    transcripts = record.get("transcript_consequences")
    if not isinstance(transcripts, list) or not transcripts:
        raise AnnotationError("vep_no_transcript_consequences")
    for tc in transcripts:
        if not isinstance(tc, dict):
            raise AnnotationError("vep_invalid_transcript")
        terms = tc.get("consequence_terms")
        domains = tc.get("domains", [])
        if (not isinstance(terms, list) or not terms or any(not isinstance(t, str) for t in terms)
                or not isinstance(tc.get("transcript_id"), str)
                or not isinstance(domains, list)
                or any(not isinstance(d, dict) or not isinstance(d.get("db"), str) or not d["db"] for d in domains)
                or (tc.get("gene_symbol") is not None and not isinstance(tc["gene_symbol"], str))):
            raise AnnotationError("vep_invalid_transcript")

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
    payload = _request_json("post", "gnomad",
        GNOMAD_API,
        json={"query": query, "variables": {"variantId": variant_id, "dataset": "gnomad_r4"}},
        timeout=timeout,
    )
    data = _graphql(payload, "variant")

    if data is None:
        # Confirmed absent from gnomAD v4 (docs/data.md's NOT_IN_GNOMAD_R4 case).
        floor = math.log10(NOT_OBSERVED_FLOOR_AF + 1e-8)
        return {
            "p4_freq_log10": _obs(floor),
            "p2_af_known_missing": _obs(0.0),
            "p2_not_observed_in_gnomad_r4": _obs(1.0),
            "p2_low_an_flag": _obs(0.0),
        }

    candidates = []
    for key in ("genome", "exome"):
        if key not in data:
            raise AnnotationError("gnomad_incomplete_frequency_response")
        part = data[key]
        if part is None:
            continue
        if not isinstance(part, dict) or not isinstance(part.get("filters"), list):
            raise AnnotationError("gnomad_invalid_filters")
        if any(not isinstance(f, str) for f in part["filters"]):
            raise AnnotationError("gnomad_invalid_filters")
        if part["filters"]:
            continue
        af, an = part.get("af"), part.get("an")
        if (isinstance(af, bool) or not isinstance(af, (int, float)) or not math.isfinite(af)
                or not 0 <= af <= 1 or type(an) is not int or an <= 0):
            raise AnnotationError("gnomad_invalid_frequency")
        candidates.append((af, an))
    if not candidates:
        # Present but filtered/non-PASS everywhere this query looked.
        return {
            "p4_freq_log10": {"value": None, "status": "filtered_or_unavailable"},
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
    payload = _request_json("post", "gnomad_constraint", GNOMAD_API,
                            json={"query": query, "variables": {"gene": gene_symbol}}, timeout=timeout)
    gene = _graphql(payload, "gene")
    if gene is not None and "gnomad_constraint" not in gene:
        raise AnnotationError("gnomad_constraint_invalid_response")
    constraint = gene["gnomad_constraint"] if gene is not None else None
    if constraint is not None and not isinstance(constraint, dict):
        raise AnnotationError("gnomad_constraint_invalid_response")
    out = {}
    for key in ("oe_lof_upper", "oe_mis_upper"):
        if constraint is not None and key not in constraint:
            raise AnnotationError("gnomad_constraint_invalid_response")
        value = constraint[key] if constraint is not None else None
        name = "p3_gene_constraint_" + key
        if value is None:
            out[name] = {"value": None, "status": "missing"}
        else:
            out[name] = _obs(value)
            if value < 0:
                raise AnnotationError("gnomad_constraint_invalid_numeric")
    return out


def fetch_conservation(chrom: str, pos: int, timeout: float = 15.0) -> dict:
    """Fetch a reference-only phyloP track; never substitute it into a model."""
    payload = _request_json("get", "ucsc", UCSC_API,
        params={"genome": "hg38", "track": "phyloP100way", "chrom": f"chr{chrom}", "start": pos - 1, "end": pos},
        timeout=timeout)
    if not isinstance(payload, dict) or not isinstance(payload.get("phyloP100way"), list):
        raise AnnotationError("ucsc_invalid_response")
    rows = payload["phyloP100way"]
    out = {"ensembl_conservation": {"value": None, "status": "unverified_source"},
           "ucsc_phyloP100way_reference": {"value": None, "status": "missing"}}
    if rows:
        if len(rows) != 1 or not isinstance(rows[0], dict) or "value" not in rows[0]:
            raise AnnotationError("ucsc_invalid_response")
        out["ucsc_phyloP100way_reference"] = _obs(rows[0]["value"])
    return out


def fetch_gpn_star(chrom: str, pos: int, ref: str, alt: str) -> dict:
    """The three GPN-Star LLR scores, via a remote Parquet row lookup (no local model, no GPU).

    Requires `polars`. Uses `songlab/gpn-star-scores` on Hugging Face — the
    exact same public precomputed dataset the frozen cache's gpn_v100/m447/p243
    values came from.
    """
    checkpoints = {"gpn_v100": "gpn-star-hg38-v100-200m", "gpn_m447": "gpn-star-hg38-m447-200m", "gpn_p243": "gpn-star-hg38-p243-200m"}
    out = {}
    for feature_name, dataset in checkpoints.items():
        path = f"{GPN_STAR_ROOT}/{dataset}/llr/llr_chr{chrom}.parquet"
        try:
            import polars as pl

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
            try:
                if row.height != 1:
                    raise AnnotationError("gpn_ambiguous_match")
                out[feature_name] = _obs(float(row["llr_calibrated"][0]))
            except Exception as exc:  # per-feature boundary includes Polars schema errors
                out[feature_name] = {"value": None, "status": "missing", "error": str(exc)[:200]}
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
