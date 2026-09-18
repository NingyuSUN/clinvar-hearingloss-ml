# Predicting a variant

`predict/` packages five trained models behind one CLI so you can look up a
variant and see all five scores plus a plain-language report. Read this whole
page before using it — the single most important fact about this tool is in
the next paragraph.

## The one thing to understand before using this

**This only works for variants that are already in the frozen, hash-verified
annotation cache shipped in `predict/bundle/`.** That cache covers roughly
4,050 of the 7,125 variants used to develop this project (the rest are
outside the four supported consequence classes, or are indels/multi-base
changes rather than single-nucleotide substitutions). `predict/examples/known_variants_7125.csv`
lists every variant this tool was built around, so you can pick a real
`variant_id` to try.

A genuinely new variant — one not already in that cache — will come back with
**null scores on all five models**, tagged `missing_annotation`. This is not
a bug or an oversight: an earlier version of this tool accepted
caller-supplied external annotations for arbitrary variants, and a review
found that a caller-supplied source name and file hash were being treated as
proof that the annotation was genuine, when they only proved internal
consistency. That let fabricated features get scored as if they were real. The
fix was to lock inference to the verified cache and make any unverified input
return null rather than a plausible-looking number. Building a genuine
live-annotation pipeline (real-time VEP, gnomAD, GPN inference, and an
AlphaGenome call for a variant nobody has annotated before) is real,
substantial future work, not something bolted on here — see
[`PROJECT_STATUS.md`](../PROJECT_STATUS.md).

## The five models

| Model | Uses |
|---|---|
| `ORIGINAL18` | The original 18 hand-engineered features only (consequence, allele frequency, gene constraint, conservation, protein-domain overlap). Same model family as the headline result in `README.md`. |
| `GPN3` | Only three GPN (Genomic Pre-trained Network) scores, from three different GPN model checkpoints. No hand-engineered features. |
| `ORIGINAL18_GPN` | The 18 original features plus the three GPN scores. |
| `FULL_UNIFIED` | 18 features + GPN + **AVI**, one pooled model across all four supported consequence classes. |
| `FULL_STRATIFIED` | Same feature set as `FULL_UNIFIED`, but routed to a protein-track sub-model for missense variants and an RNA-track sub-model for the other three classes, instead of one pooled model. |

**AVI is a composite score that blends AlphaGenome, AlphaMissense and related
outputs — it is not a pure AlphaGenome score.** There is no dedicated
AlphaGenome-only model in this release; `FULL_UNIFIED` and `FULL_STRATIFIED`
are the closest things to it, and both also include GPN and the 18 original
features.

## What a score means (and doesn't)

Every score is an **uncalibrated XGBoost classifier output, not a probability
of pathogenicity**. `predictions.jsonl` always carries `"probability": null`
and `"probability_status": "not_calibrated"` — that's deliberate, not a gap
in the output format. No one has evaluated whether 0.8 means "80% likely
pathogenic" for these models, so the tool never implies that.

Each model also reports a `research_decision` against its own frozen 5%
threshold: `higher_risk_support`, `lower_risk_support`, or `uncertain_*`. A
directional call like that is only ever reported for **missense** variants.
For the other three supported classes (`canonical_splice`, `splice_region`,
`synonymous`) the decision is always `uncertain_subtype_validation`,
regardless of the raw score — those subtypes' directional calls have not been
independently validated, so the tool does not pretend otherwise.

## Install and run

```bash
pip install -r predict/requirements-inference.txt   # numpy + xgboost only

python predict/predict_variants.py \
  --input predict/examples/input.csv \
  --bundle predict/bundle \
  --output my_output/

# Turn the structured output into a readable report:
python predict/render_variant_report.py \
  --predictions my_output/predictions.jsonl \
  --variant-id 48
```

Input CSV columns: `variant_id,assembly,chrom,pos,ref,alt`. Only GRCh38,
1-based coordinates. `--output` must be a fresh (empty or non-existent)
directory. Extra input columns are ignored (and logged in
`run_manifest.json`); clinical labels are never read as features even if a
column happens to be named `label` or similar — `predict/prediction_core.py`
actively rejects any annotation record containing one.

Outputs:
- `predictions_wide.csv` — one row per variant, all five scores/statuses/decisions side by side.
- `predictions.csv` — one row per variant × model (long format).
- `predictions.jsonl` — full structured detail per variant, including the raw GPN/AVI component scores and gene-fold routing. This is what `render_variant_report.py` reads.
- `run_manifest.json` — versions, input/bundle/code hashes, per-model status counts, and the same `limits` list summarized above.

No network calls happen during inference, and nothing is retrained. Scoring
the full 7,125-variant development list takes under 15 seconds on a modern
CPU.

## Trying it out

```bash
make predict-demo
```

runs the CLI on `predict/examples/input.csv` and then renders a report for
every variant in it, so you can see the whole pipeline end to end without
picking a variant yourself first.

## Gene-fold routing, briefly

Each model was trained with 5-fold gene-held-out cross-validation. When you
ask about a variant in a **known** gene, the tool picks the one fold that
never trained, selected, or calibrated on that gene's group — so the number
you get is not the model grading its own homework. A variant in a gene the
models have never seen gets `novel_gene_context` (fold 0, no held-out
guarantee) instead of a fabricated default.

## Limitations (see also `MODEL_CARD.md`)

- Cache-only: no scoring for variants outside the frozen cache (see above).
- Not calibrated: scores are not probabilities.
- Directional calls only for missense; other classes report uncertain by design.
- This is a fixed `seed101` research profile, not a new performance estimate — reusing it across many variants does not add new external validation evidence.
- `predict/bundle/` is hash-verified end to end (`predict.Bundle` re-checks every file's sha256 against `bundle_manifest.json` on load); a corrupted or hand-edited bundle file fails loudly rather than silently scoring incorrectly.
