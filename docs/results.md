# Results

The historical subset is coding non-truncating, not strict missense. Historical numeric
artifacts are preserved. Revised aggregate outputs are in [statistics_20260921](../results/statistics_20260921/). See [validation status](validation_status.md) for
metric aggregation and unresolved statistical-inference limitations.

Gene-held-out grouped 5-fold cross-validation, 10 frozen seeds. XGBoost.
Raw tables are in `results/`; `run_manifest.json` records the XGBoost version and
the matrix hash.

Three ways to summarise AUC:

| | what it measures |
|---|---|
| per-fold mean ± SD | average fold-model AUC on its own held-out genes — SD is **not** a CI |
| seed-averaged pooled AUC | ranking of each variant's average OOF score over the whole cohort |
| gene-group-macro | unweighted mean of within-group AUC over groups with both label classes |

## Headline cohort (N = 7,125)

| Full model (18 features) | value |
|---|---|
| AUC, per-fold mean ± SD | 0.950 ± 0.016 |
| Pooled AUC [conditional 95% interval] | **0.930 [0.892, 0.974]** |
| AUC, gene-group-macro | 0.978 |
| R90 recall (pooled) | 0.845 [0.800, 0.888] |
| R90 precision (pooled) | 0.855 [0.699, 0.972] |

Intervals use 2,000 whole-gene-group bootstrap draws with fixed OOF scores and
votes. They preserve within-group dependence but omit refitting, selection and
full CV prediction dependence. All 2,000 draws were valid for every reported metric.

The inner-validation recall is ~0.90 by construction; the outer-test recall is
~0.84 — the R90 threshold is fitted on a small inner-validation set and does not
transfer exactly.

### Additive feature ladder — per-fold AUC

| Feature set | n feat | Headline AUC | coding non-truncating AUC |
|---|---:|---:|---:|
| conservation only | 1 | 0.700 | 0.666 |
| consequence only | 4 | 0.852 | 0.531 |
| frequency only | 4 | 0.819 | 0.772 |
| allele + conservation (Base) | 5 | 0.797 | 0.676 |
| + consequence | 9 | 0.925 | 0.682 |
| + gene constraint | 11 | 0.926 | 0.707 |
| + frequency | 15 | 0.946 | 0.853 |
| **+ protein domain (Full)** | 18 | **0.950** | **0.856** |

### Ladder steps — paired per-fold difference

| Step | Headline mean Δ | Coding non-truncating mean Δ |
|---|---|---|
| Base → + consequence | +0.128 | +0.006 |
| + consequence → + constraint | +0.001 | +0.025 |
| + constraint → + frequency | +0.020 | +0.146 |
| + frequency → Full | +0.004 | +0.002 |

These are descriptive mean paired-fold differences. Repeated folds share training
data, so they are not independent experiments. The revised code does not output
ordinary t tests or significance flags; no superiority or equivalence claim is made.

### Remove one feature group from the full model

| Remove | Headline mean Δ | Coding non-truncating mean Δ |
|---|---|---|
| consequence flags | −0.085 | +0.001 |
| gene constraint | **+0.011** | **−0.054** |
| conservation | −0.012 | −0.040 |
| allele frequency | −0.010 | −0.136 |
| protein domain | −0.004 | −0.004 |
| allele shape | +0.001 | +0.003 |

Removing gene constraint increases mean fold AUC in the headline cohort and
decreases it in the coding non-truncating subset. These descriptive effects
do not establish a causal mechanism such as memorization of gene identity.

## Coding non-truncating subset ( N = 3,121, prevalence 0.32)

| Full model | value |
|---|---|
| AUC, per-fold mean ± SD | 0.856 ± 0.038 |
| Pooled AUC [conditional 95% interval] | **0.815 [0.702, 0.920]** |
| AUC, gene-group-macro | 0.941 |
| R90 recall (pooled) | 0.910 [0.858, 0.951] |
| R90 precision (pooled) | 0.525 [0.330, 0.771] |

Frequency-only has per-fold mean AUC 0.772; the matching full-model per-fold
mean is 0.856: matched-fold mean difference +0.08385. Under seed-averaged pooled
aggregation, frequency-only AUC is 0.74584 and full-model AUC is 0.81476, a
descriptive difference of +0.06893. These aggregations are not interchangeable.
No independent-fold significance is claimed. The previously displayed 0.814
was the historical row-bootstrap mean; 0.815 is the rounded empirical point
estimate from the same predictions, not a model improvement.

## Consequence-stratified recall (headline model, global R90 threshold)

| Consequence class | N | Prevalence | AUC | R90 recall | R90 precision |
|---|---:|---:|---:|---:|---:|
| truncating | 1,952 | 0.999 | — | 0.997 | 0.999 |
| canonical splice | 563 | 0.996 | — | 0.966 | 0.998 |
| coding non-truncating | 3,121 | 0.317 | 0.778 | **0.530** | 0.567 |
| non-coding / other | 1,489 | 0.051 | 0.833 | 0.158 | 0.098 |

(`AUC` is omitted here for the nearly single-class LoF strata. Recall is the fraction of pathogenic
variants in each stratum recovered at the one global R90 threshold.)

## Diagnostics

### Gene-held-out vs random row split (full model, pooled AUC)

| Split | AUC |
|---|---:|
| gene-held-out | 0.9301 |
| random row split | 0.9963 |
| gap | **0.0662** |

The comparison above uses `results/split_gap.csv`. The 0.9267 value in
`pooled_oof_by_cohort.csv` is the mean of per-seed pooled AUCs, not the AUC
of seed-averaged predictions. Neither comparison identifies a biological mechanism.

### Review-status sensitivity (conditional gene-group intervals)

| Cohort | N | Pooled AUC [conditional 95% interval] |
|---|---:|---|
| all-tier (0-star included) | 7,736 | 0.930 [0.882, 0.975] |
| ≥ 1 star (headline) | 7,125 | 0.930 [0.892, 0.974] |
| ≥ 2 star | 3,554 | 0.967 [0.942, 0.987] |

The cohorts differ in composition as well as review status. These overlapping
analyses are not independent replications and do not establish that higher review
status causes better model performance. The revised conditional intervals
are wider than the historical variant-row intervals. See the
[revision report](statistical_revision.md) for the comparison and provenance.

## Feature attribution

`docs/shap.md` has the full write-up. Headline result: TreeSHAP attribution on
the fitted models ranks the coding non-truncating feature groups in the same order as the
ablation above (frequency > gene constraint > conservation > domain, Spearman
0.94) — complementary descriptions of the same modelling exercise, not independent
replications or evidence of causal feature effects.
