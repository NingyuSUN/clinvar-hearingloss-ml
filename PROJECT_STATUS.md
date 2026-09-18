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
