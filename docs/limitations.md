# Limitations

This benchmark measures something narrower than "predicting hearing-loss variant
pathogenicity", and the headline number should not be read as that.

## The task is easier than it looks

- **VUS and Conflicting variants are excluded** by the label rule. What remains is
  a selected set of variants that ClinVar submitters could already classify
  confidently. Performance on this set is an upper bound on performance over all
  variants a clinician actually sees.
- **Consequence type near-determines the label for a third of the cohort.**
  Truncating variants are 99.9 % pathogenic and canonical-splice 99.6 %
  (ACMG `PVS1`); non-coding variants are 95 % benign. A four-flag consequence-only
  model already scores AUC 0.85. The full-cohort 0.93 is mostly loss-of-function
  identification. The scientifically meaningful number is missense (AUC 0.81).

## Circularity with the labels

- **Allele frequency co-determines the labels it is used to predict.** ClinVar
  benign calls use ACMG `BA1` (allele frequency > 5 % → stand-alone benign) and
  `BS1` (elevated frequency → benign evidence). `p4_freq_log10` is the dominant
  missense feature; part of "common → benign" is the model reproducing the rule
  the submitters applied.
- Gene-level features (gene constraint) are constant within a gene. Under
  gene-held-out CV they still generalise to unseen genes, but they cannot
  distinguish variants within a gene.

## The operating point

- A global R90 threshold reaches 90 % recall overall while missing ~50 % of
  missense pathogenic variants. Any deployment claim needs a per-consequence
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
- Interpretability (SHAP, per-gene error analysis) and any comparison to published
  variant-effect predictors are out of scope for this evaluation.
