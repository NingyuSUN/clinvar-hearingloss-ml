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
- [x] Scored-missense-only training and matched evaluation in the exploratory
      [four-arm follow-up](experiments/esm_missense_2x2/README.md), completed
      2026-10-21. Coverage is 1,328 scored missense records, not all possible
      missense variants; the earlier documentation history remains in `CHANGELOG.md`.
- [ ] Audit the conservation feature's exact provenance and the ClinVar
      phenotype-condition mapping (open item in `docs/limitations.md`).
- [x] Prespecified protein language model extension (`PLM-ESM-2026-10`,
      [protocol](docs/plm_protocol.md), [results](docs/plm_results.md)).
      Most of the primary mixed-cohort gain is reproduced by the scored-status
      control; within scored missense variants (1,328), ESM-1v adds +0.017 AUC
      [+0.005, +0.033]. The separate 2026-10-21 follow-up now also evaluates
      missense-only training; the two runs use different validation targets.

## Current experiment index — 2026-10-21 (Asia/Tokyo)

- **Objective and phase:** completed internal exploratory analysis and public
  research package for the ESM missense-training question.
- **Inputs and code:** 3,121 mixed coding non-truncating records; evaluation on
  1,328 scored missense records. Scientific base `d27ff142`; the
  [protocol](experiments/esm_missense_2x2/PROTOCOL.md), pinned input hashes,
  exact environment, release hashes and portable reproduction commands are in
  the [experiment directory](experiments/esm_missense_2x2/README.md).
- **Verified result:** missense-only Original18 + ESM versus Original18 has
  AUC 0.986222 versus 0.976574, FN 128 versus 91 of 925 P/LP, and FP 8 versus
  20 of 403 B/LB. Ranking improved while sensitivity fell at the chosen
  operating rule. Training-scope comparisons are inconclusive.
- **Evidence and limits:** aggregate results and independent arithmetic/model
  replay audits are linked in the experiment directory. Narrow validation
  gene diversity (40/50 subsets have two groups), conditional intervals and
  same-cohort reuse limit inference; external and clinical validation are absent.
- **Next step:** study score/threshold/voting behavior and threshold robustness
  using an independently specified evaluation design; do not tune on these test labels.
- **Git scope:** this package adds code, aggregate results and documentation
  to the public repository. Raw per-variant OOF predictions, fitted models and
  operational logs are not distributed. A remote configuration alone does not
  establish push status; use the branch/PR commit shown by GitHub for the
  published revision. This entry does not assert that a working checkout is
  synchronized with a remote branch.

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
- [ ] Expand the existing live-source and REF regression checks, and verify
      broader training/inference feature parity.

## Optional future research

- Threshold robustness and score/threshold/voting decomposition after the
  completed scored-missense-only training experiment.
- Per-consequence-class operating points instead of a single global R90
  threshold.
- Re-evaluation against a newer ClinVar release to check drift.
- Validating `predict_novel.py`'s five models against external labels that
  never touched training/selection/calibration, now that the live-annotation
  pipeline exists.
- A true AlphaGenome-only model, if the AlphaGenome-specific component can be
  cleanly separated out of the current composite AVI score.
