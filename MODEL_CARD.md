# Model Card: ClinVar Hearing-Loss Variant Pathogenicity Model

## Summary

An XGBoost classifier that predicts ClinVar pathogenic/benign status for
hearing-loss-gene-panel variants, evaluated under a **gene-held-out** grouped
cross-validation so that the reported performance cannot come from the model
having memorised a gene it saw at training time. The headline gene-held-out
AUC (0.930) is markedly lower than the same model's AUC under an ordinary
random split (0.996) — that ~0.07 gap is the study's central finding, not a
footnote: most of a naive evaluation's apparent accuracy is the model
recognising the gene, not learning variant-level pathogenicity signal.

This is a research/portfolio project, not a validated clinical tool. See
`docs/limitations.md` for the full limitations list.

## Intended use

- Demonstrating a leakage-aware evaluation protocol (frozen label rule,
  frozen cohort, frozen gene-group split, frozen feature list, frozen
  operating point — all fixed in `docs/methods.md` before training) for
  variant-pathogenicity modelling.
- A methods reference for anyone building a similar gene-panel pathogenicity
  classifier who wants to see the gap between random-split and
  gene-held-out AUC made explicit.

## Out-of-scope use

- **Not a clinical diagnostic tool.** No claim is made about performance on
  variants a clinician would actually need help with (VUS / Conflicting are
  excluded from training and evaluation by construction — see
  `docs/limitations.md`).
- Not validated on any cohort outside the ClinVar hearing-loss gene panel used
  here, and not re-evaluated against newer ClinVar releases.
- The R90 (90% recall) operating point is a single global threshold; per
  `README.md` it hides a large per-consequence-class gap (missense recall
  ~53% vs. truncating recall ~99.7%) and must not be read as a uniform
  90%-recall guarantee across variant classes.

## Data

- **Source**: NCBI ClinVar `variant_summary.txt.gz` (public domain), Ensembl
  VEP consequence terms, gnomAD v4.1.1 allele frequency, gnomAD v2.1.1 gene
  constraint, and an Ensembl comparative-genomics conservation score. Full
  provenance in `docs/data.md`.
- **Label**: ClinVar aggregate `ClinicalSignificance` (P/LP → 1, B/LB → 0);
  VUS and Conflicting are excluded, not imputed.
- **Cohort**: 7,125 variants (feature-complete gate ∩ ClinVar review status
  ≥ 1 star), 167 gene groups. All-tier and ≥2-star cohorts are reported as
  sensitivity analyses in `results/review_status_sensitivity.csv`.

## Performance

See `README.md` and `docs/results.md` for full tables. Headline numbers:

| Evaluation | AUC | 95% CI |
|---|---:|---:|
| Gene-held-out (the number that matters) | **0.930** | 0.924 – 0.936 |
| Ordinary random split (shown to illustrate leakage) | 0.996 | 0.995 – 0.997 |
| Coding non-truncating subset, gene-held-out | 0.814 | 0.800 – 0.829 |

## Limitations and risks

- A model using only the four consequence-class flags already reaches
  AUC 0.85 on the full cohort; most of the full-cohort 0.93 is loss-of-function
  identification, not fine-grained pathogenicity discrimination.
- Allele frequency is the dominant feature on the harder (coding
  non-truncating) subset, and ClinVar's own benign calls use frequency-based
  ACMG rules (`BA1`/`BS1`) — part of that signal is the model re-deriving the
  labelling rule rather than learning independent biology.
- Gene-level features (gene constraint) are constant within a gene; under
  gene-held-out CV they generalise to unseen genes but cannot capture
  within-gene variant-level differences.
- The conservation feature's exact provenance and the ClinVar
  phenotype-condition mapping are not fully audited (`docs/limitations.md`).
- **A documentation correction is on record**: a 2026-09-11 commit renamed
  the "coding non-truncating" cohort's occurrences of "missense" throughout
  the docs, since that 3,121-row subset is actually 1,329 missense / 1,733
  synonymous / 59 other coding variants, not pure missense. That commit was
  reverted on 2026-09-18 pending author review (see `CHANGELOG.md`); readers
  should not assume the current README wording is the final word on this
  point.

## Release decision

Public research/portfolio repository. No clinical or production use is
authorized. Any downstream use for variant triage would require, at minimum,
re-evaluation on VUS/Conflicting variants and a per-consequence-class
operating point rather than the single global R90 threshold used here.
