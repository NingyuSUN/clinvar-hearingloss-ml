# Results archive

Root-level numeric files preserve the historical frozen run. Do not overwrite
them when reproducing or revising statistics; use a new run directory.

- `ladder_headline.csv` / `ladder_missense.csv`: descriptive per-fold AUC.
- `paired_ladder_*.csv`: historical differences and **uncorrected t statistics**;
  their significance flags are superseded and must not support current claims.
- `bootstrap_ci.csv` / `summary.json`: historical **variant-row** bootstrap,
  not patient bootstrap; their intervals omit within-gene dependence.
- `stratified_headline.csv`: consequence-stratified recall and precision.
- `split_gap.csv`: seed-averaged pooled AUC comparison.
- `pooled_oof_by_cohort.csv`: means of per-seed pooled metrics, a different aggregation.
- `run_manifest.json`: original training provenance.

See [statistical revision](../docs/statistical_revision.md) for the replacement
method, versioned outputs and execution status. Historical `missense` filenames
refer to the coding non-truncating subset.

## Revised statistics

`statistics_20260921/` contains the completed post hoc reanalysis of the same
frozen OOF predictions: conditional whole-gene-group intervals, descriptive
paired effects and a provenance manifest. No models were retrained. See the
[revision report](../docs/statistical_revision.md) before interpreting the intervals.
