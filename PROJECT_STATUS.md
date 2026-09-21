# Project Status

**Status: ongoing research.** The frozen headline protocol (`P4-LABEL-SPLIT-THRESHOLD-2026-09-10`)
has been run once and is reported in `README.md` / `docs/results.md`; the
project is not closed, and open questions below are tracked as follow-on work
separately from the published run. The [documentation review](docs/validation_status.md)
clarifies cohort terminology, metric aggregation, and unresolved validation issues.

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
- [x] Experimental live annotation with offline failure regression tests:
      HTTP/GraphQL errors, variant identity, frozen AF transform, held-out-group
      integrity and honest cache membership. Only GPN3 can score compatible
      live inputs; the four conservation-dependent models are blocked.
      UCSC phyloP is reference-only. No independent external validation or
      new live service verification is claimed; see `docs/predict_novel.md`.
- [ ] Probability calibration and independent external-label validation for
      the `predict/` models (both cache-only and live-annotated); scores are
      currently uncalibrated by design.
- [ ] Reconstruct and validate the original `ensembl_conservation` source
      before enabling conservation-dependent models for live variants.
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
  vs. random-split AUC 0.996. The revised conditional headline interval is
  0.892–0.974; the coding non-truncating interval is 0.702–0.920.
- Coding non-truncating subset (N=3,121): gene-held-out AUC 0.815
  (empirical point estimate; historical row-bootstrap mean rounded to 0.814).
- A four-consequence-flag-only model already reaches AUC 0.85 on the full
  cohort, indicating strong consequence-category separation without quantifying
  its share of the full-model AUC.

## Validation follow-up

- [x] Implement conditional gene-group intervals and descriptive paired comparisons.
- [x] Execute and verify the frozen-prediction reanalysis (5 cohorts, 2,000
      gene-group draws each); see [statistical revision](docs/statistical_revision.md).
- [ ] Add regression checks for live-source errors, REF validation, and
      training/inference feature parity.

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
