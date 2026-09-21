# Methods — the frozen protocol

Every choice below is fixed before a model is trained. `src/hlpath/protocol.py` is
the executable form.

## Label

Source column: ClinVar **aggregate** `ClinicalSignificance`, not `ClinSigSimple`
(which flags a variant as pathogenic if any single submission is P/LP and would
mislabel aggregate-benign variants).

| Aggregate value | Label |
|---|---|
| `Pathogenic`, `Likely pathogenic`, `Pathogenic/Likely pathogenic` | 1 |
| `Benign`, `Likely benign`, `Benign/Likely benign` | 0 |
| Uncertain significance, Conflicting, risk factor, drug response, … | excluded |

"Likely" is kept (standard practice). Multi-term strings (`Pathogenic; association`)
take the first term.

## Cohort

Feature-complete gate — a variant is eligible only if all of:

1. it has a label;
2. its consequence annotation is reliable (a single or consensus transcript);
3. gnomAD gene constraint is present;
4. its gnomAD frequency has a definitive status (observed-PASS, or confirmed
   absent — not "query failed");
5. the coding / non-coding call was made;
6. `ensembl_conservation` is present — **rows without it are dropped, not imputed**;
7. `GeneSymbol` names a single gene.

**Headline cohort** = feature-complete ∩ review status ≥ 1 star (drops
"no assertion criteria provided"). N = 7,125.

Sensitivity cohorts, re-run with the same protocol:

| Cohort | Definition | N |
|---|---|---:|
| all-tier | feature-complete (0-star included) | 7,736 |
| ≥ 2-star | feature-complete ∩ review star ≥ 2 | 3,554 |
| expert panel | feature-complete ∩ review star ≥ 3 | 116 (spot check only) |

ClinVar gold stars: `no assertion criteria provided` = 0; `single submitter` /
`conflicting classifications` = 1; `multiple submitters, no conflicts` = 2;
`reviewed by expert panel` = 3; `practice guideline` = 4.

## Gene split

The split unit is a **gene group**: a connected component of `GeneSymbol`
co-occurrence. If two genes ever appear together in a multi-gene annotation, the
whole component is held out together, so a gene is never partly in train and
partly in test.

- **Outer**: grouped 5-fold cross-validation on gene groups. Folds are balanced on
  both total row count and pooled label prevalence (a longest-processing-time
  assignment plus one refinement pass; deterministic given the seed).
- **Inner**: inside each training fold, a balanced gene-group hold-out
  (~20 % of the training groups) is used for iteration selection and the R90
  threshold. The outer-test genes are never used for anything.
- An ordinary stratified random row split is run only as a leakage diagnostic.

## Model

XGBoost binary classifier: `eta` 0.03, `max_depth` 5, `subsample` /
`colsample_bytree` 0.9, `lambda` 1.0.

- Trained to a 700-round cap with **no in-loop early stopping**.
- The best iteration is the **argmin of the inner-validation log-loss curve**
  (minimum 10 iterations). Log-loss, not AUC — with gene-grouped inner validation a
  held-out group is often single-consequence, so validation AUC saturates while
  the model is still uncalibrated.
- Missing values are imputed inside every split: continuous features → the
  inner-train median; binary indicator flags → 0.

## Operating point

**R90**: the smallest-false-positive threshold with recall ≥ 0.90 on the inner
validation. The threshold is frozen there and applied to the outer-test
predictions **of the same fitted model** (refitting on more data shifts the
probability calibration, so a transferred threshold would not be valid).

`metrics.r90_threshold` returns a failure flag — never a fabricated recall = 1 —
when the validation set is single-class or the scores are not finite.

## Features (18, six groups)

| Group | Columns |
|---|---|
| consequence | `p1_has_frameshift`, `p1_has_stop_gained`, `p1_has_canonical_splice`, `p1_has_splice_region` |
| allele shape | `variant_length`, `has_del`, `has_ins`, `has_dup` |
| conservation | `ensembl_conservation` |
| gene constraint | `p3_gene_constraint_oe_lof_upper`, `p3_gene_constraint_oe_mis_upper` |
| frequency | `p4_freq_log10`, `p2_af_known_missing`, `p2_not_observed_in_gnomad_r4`, `p2_low_an_flag` |
| protein domain | `p3_in_specific_domain`, `p3_specific_domain_count`, `p3_variant_is_coding` |

`GeneSymbol` is never a feature. The ablation adds these groups one at a time and
also removes each group from the full model (`FEATURE_SETS` in `protocol.py`).

## Statistics — revised reporting

The original model-fitting protocol above is unchanged. The 2026-09-21
post hoc statistical revision changes uncertainty estimation and reporting:

- Average OOF scores across seeds once per variant. Apply each original fold's
  validation-selected threshold before voting; a vote of exactly 0.5 is positive,
  preserving the historical rule. This is a repeated-OOF summary, not evaluation
  of a single deployed model.
- Report the empirical pooled AUC, recall and precision as point estimates.
  The bootstrap mean is a separate diagnostic, not the point estimate.
- Resample biological gene-associated groups as whole clusters, with replacement,
  for 2,000 draws. Each draw samples G groups and retains all their variants,
  including repeated copies when a group is sampled repeatedly. Metrics remain
  variant-weighted; cluster sizes and total variant counts can vary by draw.
- Report percentile 95% intervals conditional on fixed models and predictions.
  This preserves within-group dependence but does not capture retraining,
  model-selection uncertainty or all dependence induced by overlapping CV fits.
  It is not a confidence interval for clinical performance on new cohorts.
- Count undefined draws per metric. AUC requires both classes, recall requires
  positives and precision requires positive calls. An interval is withheld if
  fewer than 100 or fewer than 90% of requested draws are valid.
- Compare feature sets using matched-fold mean differences and descriptive SDs.
  Seed-mean ranges are descriptive too. No independent-fold standard error,
  t statistic, p-value or significance flag is produced. No equivalence claim
  follows from a small or inconsistent difference.
- Validate one prediction per variant per seed, complete seed/fold coverage,
  stable labels/groups, finite successful thresholds and recorded test counts.
  The frozen-data adapter additionally checks the original 10 × 5 protocol,
  gene isolation and alignment between compared models.

Historical row-bootstrap intervals and uncorrected t statistics remain archived
under `results/`; they are superseded for inference. New outputs are versioned
separately. See [statistical revision](statistical_revision.md) for execution status.

Methodological context: [cluster bootstrap for correlated ROC data](https://pmc.ncbi.nlm.nih.gov/articles/PMC8586066/)
and [dependence in repeated-CV comparisons](https://scikit-learn.org/1.5/auto_examples/model_selection/plot_grid_search_stats.html).
