# Limitations

This benchmark measures something narrower than "predicting hearing-loss variant
pathogenicity", and the headline number should not be read as that.

## The task is easier than it looks

- **VUS and Conflicting variants are excluded** by the label rule. What remains is
  a selected set of variants that ClinVar submitters could already classify
  confidently. Performance on this selected set does not establish performance over all
  variants encountered in practice.
- **Consequence type near-determines the label for a third of the cohort.**
  Truncating variants are 99.9 % pathogenic and canonical-splice 99.6 %
  (ACMG `PVS1`); non-coding variants are 95 % benign. A four-flag consequence-only
  model already scores AUC 0.85. This demonstrates strong separation by consequence category; it does not
  quantify that category's contribution to the full model's AUC. The coding non-truncating subset has pooled AUC about 0.81; this is not
  a strict-missense estimate.

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

- The validation-selected global R90 threshold gives outer-test pooled recall
  about 84.5 %, while missing about 47 % of coding non-truncating pathogenic
  variants. Any deployment claim needs a per-consequence
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
- SHAP attribution is reported in `docs/shap.md`; it is not causal evidence.
  Per-gene error analysis and comparisons to published variant-effect predictors
  are not established by the headline evaluation.


## Statistical reporting

Revised intervals resample whole gene groups conditional on saved predictions;
they do not capture refitting, model selection or all CV prediction dependence.
Paired-fold comparisons are descriptive, without significance or equivalence
claims. See [statistical revision](statistical_revision.md) for execution status.
