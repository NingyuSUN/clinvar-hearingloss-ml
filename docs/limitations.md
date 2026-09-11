# Limitations

> **Stratification correction (2026-09-11):** The historical 3,121-row `coding_nontruncating` subset includes 1,329 missense, 1,733 synonymous and 59 other coding variants. Its AUC 0.814, R90, ablation and SHAP results are **not strict-missense results**. Historical filenames containing `missense` are retained for reproducibility. See [stratification details](stratification.md).

This benchmark measures something narrower than "predicting hearing-loss variant
pathogenicity", and the headline number should not be read as that.

## The task is easier than it looks

- **VUS and Conflicting variants are excluded** by the label rule. What remains is
  a selected set of variants that ClinVar submitters could already classify
  confidently. Performance on this selected set does not establish performance on
  the broader uncertain-variant population.
- **Consequence type near-determines the label for a third of the cohort.**
  Truncating variants are 99.9 % pathogenic and canonical-splice 99.6 %
  (ACMG `PVS1`); non-coding variants are 95 % benign. A four-flag consequence-only
  model already scores AUC 0.85. The full-cohort 0.93 is mostly loss-of-function
  identification. The coding non-truncating AUC (0.81) includes synonymous variants; strict missense performance has not been evaluated.

## Circularity with the labels

- **Allele frequency co-determines the labels it is used to predict.** ClinVar
  benign calls use ACMG `BA1` (allele frequency > 5 % → stand-alone benign) and
  `BS1` (elevated frequency → benign evidence). `p4_freq_log10` is the dominant
  coding non-truncating feature; part of "common → benign" is the model reproducing the rule
  the submitters applied.
- Gene-level features (gene constraint) are constant within a gene. Under
  gene-held-out CV they still generalise to unseen genes, but they cannot
  distinguish variants within a gene.
- The two gene-constraint features (`oe_lof_upper`, `oe_mis_upper`) are
  correlated with each other (r = 0.62); their individual SHAP attribution
  direction is not reliable (see `docs/shap.md`). Their *combined* importance is
  well supported (agrees with the ablation), but no directional claim should be
  made for either one alone without a dedicated analysis.

## The operating point

- The threshold targets 90 % recall on inner validation; observed outer-test
  recall is ~84 % overall and ~53 % within the coding non-truncating subgroup of that headline model. The separately trained coding-subset model has different pooled recall/precision (~91 %/~53 %); see [stratification details](stratification.md). Any deployment claim needs a per-consequence
  operating point, not one global threshold.
- The R90 threshold is chosen on a small inner-validation set and does not
  transfer exactly: inner-validation recall ~0.90, outer-test recall ~0.84.

## Unaudited inputs

- **`ensembl_conservation`** — the score's release, coordinate mapping and exact
  definition are not re-verified in this cycle. Variants lacking it are dropped
  from the cohort rather than imputed.
- **Phenotype mapping** — the hearing-loss selection is by MedGen ID / phenotype
  string on the ClinVar row. Whether every retained P/LP assertion is *for the
  hearing-loss phenotype* (rather than another condition of the same gene) is not
  individually checked.
- **21 multi-gene coding overlaps** are excluded by the single-gene gate rather
  than resolved.

## Scope

- Not a clinical diagnostic tool. Research / methods work.
- Gene constraint, conservation and domain membership are informative features,
  not causal evidence of pathogenicity.
- SHAP results explain the frozen models. Strict-missense re-training and
  comparisons to published variant-effect predictors remain separate work.
