# ESM in a matched missense-only evaluation

Adding a frozen ESM-1v score improved missense ranking after missense-only training, but increased false negatives under the prespecified R90 and voting rule. This exploratory follow-up tests whether restricting the downstream classifier's training population to scored missense variants resolves the ranking–detection trade-off observed in the earlier coding non-truncating analysis. It does not replace that analysis.

## Design

The four arms share the same **1,328 scored missense records: 925 pathogenic/likely pathogenic (P/LP), 403 benign/likely benign (B/LB), and 135 gene groups**. These are variant records, not patients.

| Arm | Inner training population | Predictors |
| --- | --- | --- |
| A | Mixed coding non-truncating | Original 18 features |
| B | Mixed coding non-truncating | Original 18 features + `esm1v_llr` |
| C | Scored missense only | Original 18 features |
| D | Scored missense only | Original 18 features + `esm1v_llr` |

The mixed source cohort contains 3,121 records (989 P/LP, 2,132 B/LB; 148 gene groups). The 1,793 records outside the scored subset comprise 1,733 synonymous, 59 other coding records, and one historical missense record affecting two residues. Results do not automatically extend to unscored records or all missense variants.

All arms inherit the same mixed-cohort gene-group assignments: 10 seeds × 5 outer folds. C/D restrict the inner training rows only after splitting. Every arm selects its boosting iteration and R90 threshold on the **same missense-only inner validation records**, with all validation/test gene groups excluded from fitting. The experiment contains **200 new downstream XGBoost fits**; ESM-1v itself is neither retrained nor fine-tuned.

Because the common missense validation target differs from the earlier mixed validation target, A/B are newly fitted controls. Their estimates must not be directly subtracted from the earlier protocol's estimates to attribute a change to ESM or training scope. See [PROTOCOL.md](PROTOCOL.md) for the frozen design and parameters.

## Results

| Arm | AUC | TP | FN / 925 P/LP | FP / 403 B/LB | TN | Recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A | 0.972353 | 845 | 80 | 24 | 379 | 91.351% |
| B | 0.984591 | 795 | 130 | 9 | 394 | 85.946% |
| C | 0.976574 | 834 | 91 | 20 | 383 | 90.162% |
| D | 0.986222 | 797 | 128 | 8 | 395 | 86.162% |

For the prespecified primary comparison **D − C**, adding ESM after missense-only training gives:

| Metric | Difference | Conditional 95% interval |
| --- | ---: | ---: |
| AUC | +0.009648 | +0.003683 to +0.016476 |
| Recall | −4.000 percentage points | −7.887 to −0.812 percentage points |
| False-positive rate | −2.978 percentage points | −5.656 to −0.983 percentage points |

The net increase of 37 false negatives comprises **30 rescued P/LP records and 67 newly missed P/LP records**; 61 were missed by both C and D. Twelve B/LB false positives were removed and none were newly introduced; eight remained false positives in both arms.

The AUC, recall, false-positive-rate and specificity intervals for C − A and D − B all cross zero. These comparisons provide no clear evidence that missense-only training improves performance over mixed training; they do **not** establish equivalence. All four contrasts, including B − A, are in [results/contrasts.csv](results/contrasts.csv). Per-seed summaries show the same D − C directions for AUC, recall and false-positive rate in all ten seeds, but repeated partitions of one cohort are not ten independent studies.

## Why AUC and missed variants move differently

AUC evaluates the ranking of each variant's **mean score across ten seeds**, with half credit for ties. Classification uses each seed's own validation-derived threshold and calls a variant positive when **at least 5 of 10** out-of-fold decisions are positive; a 5–5 tie is positive. Thus ranking and error counts summarize different aspects of the predictions.

R90 chooses the largest validation threshold retaining at least 90% positive recall, minimizing false positives subject to that constraint. It is an **inner-validation target**, not a guarantee of 90% recall on held-out genes. At fixed validation scores, this threshold is determined by the positive-score order statistic. Negative validation records can influence model/iteration selection upstream, but do not directly determine that order statistic. These results do not establish that a numerically higher threshold, poorer calibration, or any single mechanism caused the extra misses; a dedicated score/threshold/voting decomposition was not performed.

The supported conclusion is a conditional ranking gain alongside a worse sensitivity operating point and fewer false positives. Removing synonymous records from C/D training does not remove that trade-off.

## Uncertainty and validation limits

- The 95% intervals use 2,000 paired percentile bootstrap draws of all 135 gene groups, seed 0, with common group multiplicities across arms. They condition on the stored fitted predictions and do not include repeated model fitting, iteration selection or threshold estimation. Secondary contrasts are exploratory and not multiplicity-adjusted.
- **40 of 50 inner missense-validation subsets contain only two gene groups**, with 170 records (149 P/LP, 21 B/LB). The other ten subsets contain 17–28 groups. Missense inner training contains 871–916 records; outer missense tests contain 19–39 groups. Matching counts do not establish identical group identities across those 40 subsets. This inherited design was retained after the pre-fit audit; no outcome-dependent reseeding was used.
- All group-isolation and numerical checks passed, but narrow validation gene diversity limits threshold transportability. Gene grouping does not eliminate every cross-group homology dependence.
- This is a post hoc exploratory study on the same ClinVar-derived cohort, with a protocol fixed before these new fits. It is not newly blinded, external, functional or clinical validation. Different training population sizes also prevent interpreting the training-scope contrast as an isolated biological intervention.

## Reproduce

Run from the repository root in a non-synchronized Linux working directory. The runner uses POSIX file locking. Use a fresh output directory; existing frozen outputs must not be overwritten. The original input data must match the SHA256 values in [verification/provenance.json](verification/provenance.json). Input preparation and data scope are documented in [the repository data guide](../../docs/data.md).

The recorded runtime is Python **3.12.14**, NumPy **2.5.3**, pandas **3.0.1**, XGBoost **3.4.1**, scikit-learn **1.7.2**, and SciPy **1.18.1**. [environment_freeze.txt](environment_freeze.txt) records the complete installed package versions; exact reproduction requires the same environment and unchanged base scientific modules from commit `d27ff142fed63a110114f8efa928cdd26bedd2df`.

```bash
python -m venv .venv-esm-2x2
source .venv-esm-2x2/bin/activate
python -m pip install -r experiments/esm_missense_2x2/environment_freeze.txt

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=10
export PYTHONPATH=src:experiments/esm_missense_2x2

python -m unittest discover -s experiments/esm_missense_2x2 -p 'test_*.py' -v

python experiments/esm_missense_2x2/run_missense_2x2.py prepare \
  --matrix data/modeling_matrix.csv.gz --plm-features data/plm/plm_features.csv \
  --out artifacts/esm_missense_2x2_reproduction --threads 10 --repo-root .

python experiments/esm_missense_2x2/audit_fit_artifacts.py \
  --repo-root . --runner-dir experiments/esm_missense_2x2 \
  --matrix data/modeling_matrix.csv.gz --plm-features data/plm/plm_features.csv \
  --run-dir artifacts/esm_missense_2x2_reproduction --stage splits \
  --out artifacts/esm_missense_2x2_reproduction/split_audit.json

nice -n 10 python experiments/esm_missense_2x2/run_missense_2x2.py fit \
  --matrix data/modeling_matrix.csv.gz --plm-features data/plm/plm_features.csv \
  --out artifacts/esm_missense_2x2_reproduction --threads 10 --repo-root .

python experiments/esm_missense_2x2/analysis.py \
  --input artifacts/esm_missense_2x2_reproduction/oof_predictions.csv \
  --out artifacts/esm_missense_2x2_reproduction/statistics \
  --bootstrap 2000 --bootstrap-seed 0

python experiments/esm_missense_2x2/audit_fit_artifacts.py \
  --repo-root . --runner-dir experiments/esm_missense_2x2 \
  --matrix data/modeling_matrix.csv.gz --plm-features data/plm/plm_features.csv \
  --run-dir artifacts/esm_missense_2x2_reproduction --stage complete \
  --replay-threads 10 \
  --out artifacts/esm_missense_2x2_reproduction/fit_audit.json

python experiments/esm_missense_2x2/audit_results.py \
  --predictions artifacts/esm_missense_2x2_reproduction/oof_predictions.csv \
  --statistics artifacts/esm_missense_2x2_reproduction/statistics --repo-root . \
  --out artifacts/esm_missense_2x2_reproduction/result_audit.json

python experiments/esm_missense_2x2/summarize_seed_stability.py \
  --statistics artifacts/esm_missense_2x2_reproduction/statistics
```

The original run passed 32 synthetic tests and a separate 54-case R90 order-statistic check. Independent audits reconstructed all 200 arm/fold roles and training-only medians, replayed 200 saved models on 400 validation/test panels (maximum score discrepancy 0), and recomputed the four confusion matrices, AUCs, primary intervals and paired error changes. These establish engineering consistency and arithmetic, not external scientific validity.

The public source differs from the frozen execution source only in six module docstrings that describe the tools for a general research audience. No executable statements, parameters or test logic changed; `fit_core.py` is byte-identical. Original result/audit hashes continue to identify the original execution files. [PACKAGE_MANIFEST.json](PACKAGE_MANIFEST.json) and [verification/provenance.json](verification/provenance.json) distinguish original and released hashes. A fresh `prepare` binds a reproduction to the released files; the released runner must not be substituted into the old run and claimed to satisfy its original byte-hash binding.

## Files and release scope

- `results/`: full-precision aggregate tables, paired intervals, seed summaries and exact C/D discordance counts. The JSON files are byte-identical to the original results. The four CSV files change only line endings from CRLF to LF to match repository normalization; every header and data field is unchanged. The manifest records both original and published hashes.
- `verification/`: original aggregate fit/result audit evidence, original prediction hashes, fixed scientific configuration and an aggregate validation-diversity check. Historical artifact hashes identify files that are not distributed in this package.
- Python source and synthetic tests: fitting, aggregation, bootstrap, independent replay and arithmetic checks.
- `PROTOCOL.md`: a public rendering of the pre-fit protocol, preserving its scientific design and pre-fit split limitation.

Per-variant prediction tables, model binaries, checkpoint directories and full run logs are not included in this release. The commands regenerate them from the pinned input data; the aggregate evidence does not by itself permit replaying historical models without those artifacts. Original analysis outputs remain preserved.
