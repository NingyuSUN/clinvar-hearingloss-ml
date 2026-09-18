# Project Status

**Status: ongoing research.** The frozen headline protocol (`P4-LABEL-SPLIT-THRESHOLD-2026-09-10`)
has been run once and is reported in `README.md` / `docs/results.md`; the
project is not closed, and open questions below are tracked as follow-on work
rather than known bugs.

## Completion checklist

- [x] Label rule, cohort gate, gene-group split, feature list and operating
      point frozen before training (`docs/methods.md`).
- [x] Headline gene-held-out evaluation run (10 seeds × 5 folds), reported
      against the random-split baseline to make the leakage gap explicit.
- [x] Feature ablation (paired per-fold) and TreeSHAP attribution, cross-checked
      against each other (`docs/shap.md`).
- [x] Sensitivity analyses: all-tier and ≥2-star review-status cohorts
      (`results/review_status_sensitivity.csv`).
- [x] Repository hardening: MIT license, CI (tests + public-artifact
      validation on every push), line-ending normalization, model card.
- [x] Downloadable five-model prediction CLI (`predict/`) with a plain-language
      report generator, covering the ~4,050-variant frozen cache; verified
      byte-identical to the original engineering pass's output before being
      folded into this repo (see `CHANGELOG.md`).
- [x] Live-annotation pipeline for genuinely novel variants (`predict_novel.py`):
      real-time Ensembl VEP + gnomAD v4.1.1/v2.1.1 + a public GPN-Star lookup
      (no local model or GPU needed) + an optional AlphaGenome API call with
      the user's own key. Every reproducible feature checked against the
      frozen cache matched exactly except conservation, which uses UCSC
      phyloP as a documented substitute for the original, never-fully-audited
      Ensembl GERP-style score — see `docs/predict_novel.md`. This path has
      had **no independent external validation** and is explicitly weaker
      than the cache-only path, not a replacement for it.
- [ ] Probability calibration and independent external-label validation for
      the `predict/` models (both cache-only and live-annotated); scores are
      currently uncalibrated by design.
- [ ] Reconstruct the original `ensembl_conservation` provenance well enough
      to either confirm or replace the UCSC phyloP substitute used by
      `predict_novel.py` with something closer to the training data's actual
      source.
- [ ] `predict_novel.py`'s VEP transcript selection skips the frozen
      pipeline's "target-gene" preference tier (the internal gene-panel list
      it depended on wasn't preserved); it currently uses MANE Select >
      canonical > protein-coding > any instead.
- [ ] A dedicated AlphaGenome-only model. The closest things that exist today
      (`FULL_UNIFIED`, `FULL_STRATIFIED`) use a composite AVI score that also
      includes AlphaMissense and related signals.
- [ ] Strict-missense-only re-evaluation of the coding non-truncating subset
      (flagged as separate follow-on work in the 2026-09-11 stratification
      commit; that commit's documentation changes were reverted 2026-09-18
      pending author review — see `CHANGELOG.md`).
- [ ] Audit the conservation feature's exact provenance and the ClinVar
      phenotype-condition mapping (open item in `docs/limitations.md`).

## Decision

Kept as a public research/portfolio repository. Not released as, and not
intended to be used as, a clinical decision-support tool.

## Evidence snapshot

- Headline: 7,125 variants, 167 gene groups, gene-held-out AUC 0.930
  (0.924–0.936) vs. random-split AUC 0.996 (0.995–0.997).
- Coding non-truncating subset (N=3,121): gene-held-out AUC 0.814
  (0.800–0.829).
- A four-consequence-flag-only model already reaches AUC 0.85 on the full
  cohort — most of the 0.93 headline is loss-of-function identification.

## Optional future research

- Strict missense-only re-training and evaluation (the raw counts for this
  split already exist from the 2026-09-11 stratification work; only the
  documentation change was reverted, not the underlying analysis).
- Per-consequence-class operating points instead of a single global R90
  threshold.
- Re-evaluation against a newer ClinVar release to check drift.
- Validating `predict_novel.py`'s five models against external labels that
  never touched training/selection/calibration, now that the live-annotation
  pipeline exists.
- A true AlphaGenome-only model, if the AlphaGenome-specific component can be
  cleanly separated out of the current composite AVI score.
