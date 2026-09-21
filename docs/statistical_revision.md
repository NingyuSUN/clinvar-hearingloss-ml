# Statistical revision — 2026-09-21

## Status

Completed locally with explicit user authorization on 2026-09-21. No model was
retrained. Five cohorts used 2,000 gene-group draws each; all draws were valid for
AUC, recall and precision. Cohort counts, headline/subset confusion matrices and
headline/random pooled AUC matched the historical reference. The implementation
passed independent design/code review, independent numerical/provenance checks,
and 29 targeted regression/artifact tests (9.91 seconds). The reviewer independently
recomputed AUC and confusion-derived metrics for all five cohorts and verified
12 input, 4 code and 7 output hashes.

Versioned outputs: [aggregate tables and manifest](../results/statistics_20260921/).

| Cohort | Variants | Gene groups | Empirical AUC [conditional 95% interval] |
|---|---:|---:|---:|
| headline | 7,125 | 167 | 0.930 [0.892, 0.974] |
| coding_nontruncating | 3,121 | 148 | 0.815 [0.702, 0.920] |
| random_split | 7,125 | 167 | 0.996 [0.994, 0.998] |
| all_tier | 7,736 | 183 | 0.930 [0.882, 0.975] |
| ge2star | 3,554 | 141 | 0.967 [0.942, 0.987] |

The headline interval changed from the historical row-bootstrap [0.924, 0.936]
to [0.892, 0.974]; the coding non-truncating interval changed from [0.800, 0.829]
to [0.702, 0.920]. Uncertainty is substantially wider when whole biological groups
are resampled. Point predictions and confusion matrices are unchanged.

The subset empirical AUC is 0.8147647, which rounds to **0.815**. The old displayed
0.814 came from the row-bootstrap mean (0.8144374), not the original-sample point
estimate. This is a reporting correction, not an improvement in the model.

## What changes

`bootstrap_ci` now requires biological gene groups and resamples complete groups.
It returns the empirical estimate, a separate bootstrap mean, conditional
percentile limits and valid/undefined draw counts. `paired_delta` now returns
matched-fold descriptive differences without ordinary t tests or significance
flags. See [methods](methods.md) for the estimand and limitations.

OOF aggregation rejects duplicate rows, inconsistent identities/labels/groups,
missing thresholds and incomplete coverage. Random-split output now keeps
biological `gene_group` separately from the `split_group` used for assignment.

The new frozen-data adapter reads archived analysis inputs only. Its code comes
from this repository. It checks the 101–110 seeds × five folds, gene isolation,
matched evaluation rows and reference cohort/confusion/AUC values. It never
runs the archived training or analysis implementation.

## Reproduce without retraining

Install this repository (`pip install -e .`) and pytest, then run with the archived
`frozen_rerun_2026-09-10` directory as `--source`:

```bash
python scripts/revise_frozen_statistics.py --source /path/to/frozen_rerun_2026-09-10 --output runs/statistics_20260921 --boot 2000
python -m pytest tests/test_statistics.py tests/test_statistics_artifacts.py tests/test_public_artifacts.py -q
```

The output directory must be new or empty. Outputs include aggregate cluster
intervals, descriptive paired-fold deltas, same-aggregation pooled comparisons,
and a manifest hashing input files, source code and generated outputs. Variant
level prediction files are not copied into the public result directory.

## Interpretation boundaries

This is a post hoc reanalysis of fixed predictions, not a new training run or an
independent replication. Group resampling does not quantify all CV dependence,
retraining variation, source drift or clinical deployment uncertainty. Repeated
seeds reuse the same variants. No formal model superiority or equivalence test
is claimed. Probability calibration and independent external validation remain
outside the completed research scope.

Historical `results/summary.json`, `bootstrap_ci.csv` and `paired_ladder_*.csv`
are preserved byte-for-byte. Their row-bootstrap intervals and uncorrected
significance flags must not be used as the revised inference.
