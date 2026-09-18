"""Render a human-readable Markdown report per variant from predict_variants.py output.

This is the piece the original engineering pass explicitly deferred ("LLM
deferred" in the project's working notes): predict_variants.py emits structured
JSON, but nothing turned it into something a person can read at a glance. This
script does that deterministically, with no model calls of its own — it only
narrates numbers that are already in predictions.jsonl.

Usage:
    python render_variant_report.py --predictions out/predictions.jsonl [--variant-id 48] [--out report.md]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

MODEL_ORDER = ["ORIGINAL18", "GPN3", "ORIGINAL18_GPN", "FULL_UNIFIED", "FULL_STRATIFIED"]

MODEL_INFO = {
    "ORIGINAL18": {
        "label": "18 features only",
        "short": "18 hand-engineered features",
        "description": (
            "The original 18 hand-engineered features: consequence flags, "
            "allele frequency, gene constraint, conservation, protein-domain "
            "overlap. No GPN or AlphaGenome-family input."
        ),
    },
    "GPN3": {
        "label": "GPN only (3 scores)",
        "short": "3 GPN checkpoint scores, no hand-engineered features",
        "description": (
            "Only the three GPN (Genomic Pre-trained Network) scores below, "
            "from three different GPN model checkpoints. No hand-engineered "
            "features at all."
        ),
    },
    "ORIGINAL18_GPN": {
        "label": "18 features + GPN",
        "short": "18 features + 3 GPN scores",
        "description": "The 18 original features plus the three GPN scores.",
    },
    "FULL_UNIFIED": {
        "label": "18 features + GPN + AVI (pooled)",
        "short": "18 features + GPN + AVI, one pooled model",
        "description": (
            "The 18 original features + GPN scores + AVI, one pooled model "
            "across every supported consequence class. AVI is a composite "
            "score that itself blends AlphaGenome, AlphaMissense and related "
            "outputs — it is NOT a pure AlphaGenome score. There is no "
            "dedicated AlphaGenome-only model in this release."
        ),
    },
    "FULL_STRATIFIED": {
        "label": "18 features + GPN + AVI (stratified)",
        "short": "18 features + GPN + AVI, protein/RNA-track sub-models",
        "description": (
            "Same feature set as FULL_UNIFIED, but routed to a protein-track "
            "sub-model for missense variants and an RNA-track sub-model for "
            "the other three supported classes, instead of one pooled model."
        ),
    },
}

COMPONENT_LABELS = {
    "gpn_v100": "GPN score (vertebrate-100 checkpoint)",
    "gpn_m447": "GPN score (mammal-447 checkpoint)",
    "gpn_p243": "GPN score (primate-243 checkpoint)",
    "avi": "AVI (composite AlphaGenome/AlphaMissense-family score)",
    "splice_sites": "Splice-site mechanism score (annotation only, not a model input)",
}

STATUS_EXPLANATIONS = {
    "scored_research": "Scored.",
    "missing_annotation": "This variant is not in the frozen verified cache, so nothing is scored. "
    "This tool only replays known, pre-annotated ClinVar hearing-loss variants — "
    "see docs/predict.md for why it cannot score an arbitrary novel variant.",
    "unverified_external_annotation": "An externally supplied annotation was provided but could not be "
    "verified against the frozen cache, so it was rejected rather than trusted.",
    "unsupported_non_snv": "Only single-nucleotide substitutions are supported; this input is an "
    "indel or a multi-base change.",
    "unsupported_variant_type": "The variant's consequence class is outside the four supported types "
    "(missense, canonical_splice, splice_region, synonymous).",
    "unsupported_assembly": "Only GRCh38 coordinates are supported.",
    "unsupported_chromosome": "The chromosome name was not recognised.",
    "invalid_position": "The position is not a valid positive integer.",
    "invalid_allele": "REF/ALT must be composed of A/C/G/T only.",
    "identical_alleles": "REF and ALT are identical; this is not a variant.",
    "duplicate_variant_id": "Another input row reused this variant_id.",
    "duplicate_variant_key": "Another input row normalizes to the same chrom-pos-ref-alt.",
    "reference_mismatch": "The supplied REF does not match the frozen/reference-verified base at this position.",
    "unresolved_gene_identity": "The cached gene symbol(s) for this variant could not be resolved to a "
    "single gene-group.",
    "unsupported_group_mapping": "This variant's genes span more than one historical gene-group, or mix "
    "known and unrepresented genes.",
}

DECISION_EXPLANATIONS = {
    "higher_risk_support": "Research support for the higher-risk (pathogenic-leaning) class at this model's frozen 5% threshold.",
    "lower_risk_support": "Research support for the lower-risk (benign-leaning) class at this model's frozen 5% threshold.",
    "uncertain_both": "The score falls in both classes' threshold bands — treated as uncertain, not a tie-break call.",
    "uncertain_empty": "The score falls in neither class's threshold band — treated as uncertain.",
    "uncertain_subtype_validation": "This variant's consequence class is not missense. Directional research "
    "calls are only reported for missense; other supported classes always report uncertain here because "
    "their directional calls have not been validated.",
    "uncertain_unavailable": "Not scored — see the status/reason above.",
}


def fmt_score(value):
    return "—" if value is None else f"{value:.4f}"


def fmt_component(name, comp):
    label = COMPONENT_LABELS.get(name, name)
    if not comp or comp.get("status") != "observed":
        return f"- **{label}**: not observed"
    return f"- **{label}**: {comp['value']:.4g}"


def render_variant(record: dict) -> str:
    lines = []
    vid = record.get("variant_id", "?")
    key = record.get("variant_key") or f"{record.get('chrom')}-{record.get('pos')}-{record.get('ref')}-{record.get('alt')}"
    lines.append(f"## Variant `{vid}` ({key}, {record.get('assembly', 'GRCh38')})")
    lines.append("")

    genes = record.get("gene_symbols") or []
    variant_type = record.get("variant_type")
    status = record.get("annotation_status")

    lines.append("### Identity and annotation")
    lines.append(f"- **Gene(s)**: {', '.join(genes) if genes else 'unknown'}")
    lines.append(f"- **Consequence class**: {variant_type or 'unresolved'}")
    lines.append(f"- **Annotation status**: `{status}`")
    if status and status != "verified_cache":
        note = STATUS_EXPLANATIONS.get(status, "See run_manifest.json `limits` for what this status means.")
        lines.append(f"  - {note}")
    lines.append(
        f"- **Cache membership**: "
        f"{'in the 7,125-variant development cohort' if record.get('in_development_7125') else 'not in the development cohort'}; "
        f"{'in the 4,050-variant common-feature cohort' if record.get('in_common_4050') else 'not in the common-feature cohort'}"
    )
    context = record.get("group_context")
    if context:
        context_note = {
            "known_gene_held_out": "this variant's gene is known; the fold used never trained, "
            "selected, or calibrated on that gene's group.",
            "novel_gene_context": "this variant's gene was not in the training data, so there is no "
            "held-out guarantee for it — treat any score for it with extra caution.",
        }.get(context, context)
        lines.append(f"- **Gene-fold context**: `{context}` — {context_note}")
    lines.append("")

    components = record.get("components") or {}
    if components:
        lines.append("### Raw component scores (not model outputs)")
        for name in ("gpn_v100", "gpn_m447", "gpn_p243", "avi", "splice_sites"):
            if name in components:
                lines.append(fmt_component(name, components[name]))
        lines.append("")

    lines.append("### Model scores")
    lines.append(
        "All scores are **uncalibrated classifier outputs, not probabilities of "
        "pathogenicity** — a higher number means more research support for the "
        "higher-risk class at this model's own frozen threshold, not \"X% chance "
        "of being pathogenic\"."
    )
    lines.append("")
    lines.append("| Model | What it uses | Score | Research decision |")
    lines.append("|---|---|---:|---|")
    models = record.get("models") or {}
    for name in MODEL_ORDER:
        m = models.get(name, {})
        info = MODEL_INFO[name]
        decision = m.get("research_decision", "uncertain_unavailable")
        lines.append(f"| **{info['label']}** | {info['short']} | {fmt_score(m.get('score'))} | `{decision}` |")
    lines.append("")

    lines.append("### What each model actually is")
    for name in MODEL_ORDER:
        lines.append(f"- **{MODEL_INFO[name]['label']}**: {MODEL_INFO[name]['description']}")
    lines.append("")

    lines.append("### What each decision means for this variant")
    for name in MODEL_ORDER:
        m = models.get(name, {})
        if m.get("status") != "scored_research":
            continue
        decision = m.get("research_decision")
        explanation = DECISION_EXPLANATIONS.get(decision, decision)
        imputed = m.get("imputed_features") or []
        extra = f" (imputed features: {', '.join(imputed)})" if imputed else ""
        lines.append(f"- **{MODEL_INFO[name]['label']}**: {explanation}{extra}")
    lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predictions", required=True, help="Path to a predictions.jsonl from predict_variants.py")
    ap.add_argument("--variant-id", help="Render only this variant_id (default: all rows)")
    ap.add_argument("--out", help="Write to this file instead of stdout")
    args = ap.parse_args()

    records = [json.loads(line) for line in Path(args.predictions).read_text().splitlines() if line.strip()]
    if args.variant_id:
        records = [r for r in records if r.get("variant_id") == args.variant_id]
        if not records:
            print(f"No variant_id={args.variant_id!r} found in {args.predictions}", file=sys.stderr)
            sys.exit(1)

    header = (
        "# Variant prediction report\n\n"
        "> Research tool. Scores are uncalibrated and not probabilities of "
        "pathogenicity. Only variants already in the frozen verified cache can "
        "be scored — see `docs/predict.md` and `MODEL_CARD.md`.\n\n"
    )
    body = "\n---\n\n".join(render_variant(r) for r in records)
    text = header + body + "\n"

    if args.out:
        Path(args.out).write_text(text)
    else:
        print(text)


if __name__ == "__main__":
    main()
