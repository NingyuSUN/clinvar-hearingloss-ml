# Validation status and reporting boundaries

## Cohorts and metrics

The 3,121-row subset is **coding non-truncating**, not strict missense.
Historical filenames and JSON keys containing `missense` remain unchanged.
Prediction-bundle missense routing uses a separate definition.

Pooled seed-averaged AUC, per-fold mean AUC and the mean of per-seed pooled AUCs
answer different questions. Compare feature sets using the same aggregation.
Frequency-only 0.772 and full-model 0.856 are matching per-fold means; neither
should be subtracted from the full-model pooled empirical AUC 0.81476.

R90 is selected on inner validation. Headline pooled recall is about 0.845;
at those thresholds, coding non-truncating recall is about 0.530. The separately
fitted subset model answers a different question.

## Statistical revision

[Completed revision and reproduction](statistical_revision.md) documents the replacement
of variant-row bootstrap with conditional gene-group bootstrap. Original point
estimates are kept distinct from bootstrap means. Paired-fold comparisons now
report descriptive effects without independent-fold t tests or significance flags.
Historical results are preserved; old intervals are not the revised inference.

The revised intervals condition on saved predictions. They do not estimate
retraining, model-selection or full new-cohort uncertainty, and they do not
remove every dependence induced by overlapping CV training sets.

## Exploratory ESM missense-training follow-up

The [2026-10-21 four-arm experiment](../experiments/esm_missense_2x2/README.md)
compares all arms on the same 1,328 scored missense variants, with newly fitted
mixed-training controls and a common missense-only inner validation target.
After missense-only training, adding ESM improves AUC by 0.009648 but increases
FN from 91 to 128 of 925 P/LP, while FP decreases from 20 to 8 of 403 B/LB.
Ranking and the validation-threshold voting decision are different summaries.
These results do not establish poorer probability calibration or a unique
threshold mechanism. Training-scope contrasts are inconclusive, not equivalent.

The experiment's 32 synthetic tests and independent saved-model/arithmetic
audits check implementation and result consistency. They do not supply external
validation. Forty of fifty inner validation subsets contain only two missense
gene groups. Conditional gene-group intervals omit retraining and threshold
selection, and the post hoc study reuses the original ClinVar cohort.

## Research boundaries

- Probability calibration and independent external validation are outside the
  completed evaluation. Scores must not be interpreted as clinical probabilities.
- Frozen-cache replay is a research demonstration, not external validation.
- Live annotation blocks conservation-dependent models; UCSC phyloP remains
  reference-only because source/scale compatibility has not been established.
- Offline regression checks cover source failures, reference identity and
  frequency parity with the frozen cache. They do not verify current service
  availability or broad training/inference feature parity.

[Project status](../PROJECT_STATUS.md) · [Methods](methods.md) ·
[Results](results.md) · [Live prediction](predict_novel.md)
