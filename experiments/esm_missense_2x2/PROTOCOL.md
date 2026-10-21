# ESM-MISSENSE-2X2-20261021

This is the public rendering of the specification frozen before the new model fits on 2026-10-21. The follow-up was motivated by previously observed results and remains exploratory; it is not a retroactive amendment to the original PLM-ESM-2026-10 analysis. Publication editing removes execution-specific administration without changing the scientific design. Original protocol and source hashes are recorded in `PACKAGE_MANIFEST.json`.

## Question and fixed inputs

Does adding the frozen ESM-1v score to the original 18 features improve discrimination and pathogenic-class detection when the downstream classifier is trained only on eligible missense variants? Compare training scopes on identical gene-held-out missense test records.

Scientific base: commit `d27ff142fed63a110114f8efa928cdd26bedd2df`. Inputs: `data/modeling_matrix.csv.gz` and `data/plm/plm_features.csv`. Before fitting, verify the input hashes against the original seal, join one-to-one by `VariationID`, and verify scored flags and finite ESM-1v values. Do not change transcripts, labels, gene groups or features. Reuse the frozen ESM scores; do not rerun or fine-tune the protein model.

The source `coding_nontruncating` cohort has 3,121 records (989 P/LP, 2,132 B/LB; 148 gene groups). The scored missense subset has 1,328 records (925 P/LP, 403 B/LB; 135 gene groups). The excluded historical missense record changes two residues. The evaluation does not automatically generalize to unscored variants or all possible missense variants.

| Arm | Inner training data | Features |
| --- | --- | --- |
| A | Mixed coding non-truncating | Original18 |
| B | Mixed coding non-truncating | Original18 + `esm1v_llr` |
| C | Scored missense only | Original18 |
| D | Scored missense only | Original18 + `esm1v_llr` |

## Splitting, fitting and operating rule

1. Use the original `balanced_group_folds` algorithm on the full mixed cohort with seeds 101–110 and five outer folds. Apply the original inner partition rule to the outer training groups: five-way group allocation, fold 0 held out, seed `seed*97+outer_fold`. Freeze assignments before fitting. Filter C/D training rows only after assignment; do not generate separate missense splits.
2. Outer test gene groups, including their synonymous records, are absent from every arm's training and validation. Inner validation groups are absent from fitting. Require exact shared outer/inner assignments and evaluated missense IDs across arms, with both classes in every arm's inner training, common missense inner validation and missense outer test. Fail visibly on an invalid split; do not reseed after observing results.
3. A/B fit on mixed inner training; C/D fit on scored-missense inner training. All four select iterations and thresholds on the same missense-only inner validation records. Fit preprocessing separately on each arm's own inner training: continuous features use their training median; binary features use zero. Missing ESM in mixed B training follows the continuous-feature median rule. No validation/test-derived medians.
4. Retain original XGBoost parameters: binary logistic, learning rate 0.03, depth 5, subsample 0.9, feature subsample 0.9, L2 penalty 1, minimum child weight 1, histogram trees, and the original seed. Fit 700 rounds; choose the first minimum missense-validation log loss from zero-based index 10 onward. No additional hyperparameter search and no refit on inner validation. Use ten CPU threads and one fit at a time.
5. Choose the inner-validation R90 threshold with the fewest false positives subject to recall ≥0.90, breaking ties at the highest threshold. Equivalently, at fixed scores this is the largest threshold retaining at least 90% of positive records. Freeze it; classify the outer missense test score as positive at `score >= threshold`. Negative validation scores affect upstream model/iteration selection but not the fixed-score positive-order-statistic threshold directly. The inner target is not a held-out recall guarantee; do not tune against test outcomes.
6. The common missense validation target differs from the original mixed validation target. A/B are newly fitted controls, not copied historical estimates. There are 4 × 10 × 5 = 200 fits, with each eligible missense record evaluated once per seed per arm.

## Prediction sealing and analysis

Save full-precision out-of-fold (OOF) scores, thresholds, decisions, labels, groups, record IDs, seeds and fold indices. Retain per-fold sample/class counts, iteration selection, training medians, model files and validation curves. Bind them to the frozen input, source, environment and split hashes. Do not inspect performance summaries during fitting. Seal the completed prediction artifacts with SHA256 before statistical analysis; preserve prior outputs and do not alter run semantics to obtain a preferred result.

The primary comparison is **D − C**. Additional exploratory comparisons are B − A, C − A and D − B. Report all four arms.

For each variant, AUC uses its mean score over ten seeds; tied scores receive half credit. Decisions use the ten fold-specific threshold decisions, declaring positive at **at least five positive votes**. Report unique-record denominators, AUC, TP/FN/FP/TN, recall, precision, false-positive rate and specificity. Describe per-seed metrics without treating the seeds as independent data sets. Tabulate exact C/D discordant calls separately for P/LP and B/LB: rescued and newly missed records, and removed and newly introduced false positives.

Use 2,000 paired bootstrap draws (seed 0), sampling all 135 complete gene groups with replacement. Use identical group multiplicities for all arms. Calculate paired metric differences and percentile 95% intervals. Report undefined draws; withhold intervals if fewer than 90% of draws are valid, with a minimum of 100. These intervals condition on the fixed fitted predictions and omit full refitting, iteration-selection and threshold-estimation uncertainty. Secondary comparisons are exploratory and not multiplicity-adjusted.

A ΔAUC interval above zero supports conditional ranking improvement. A zero-crossing interval is inconclusive, not equivalence. Interpret recall differences with actual FN counts and intervals. Fewer FP together with more FN is a trade-off, not an overall decision improvement. Do not select a deployment threshold from test labels or introduce a new non-inferiority margin.

## Checks and limitations

Engineering checks require a unique merge, fixed hashes/versions, expected cohorts, group disjointness, matched targets, training-only imputation, 200 successful fits and 53,120 OOF rows, ten repeated predictions per variant/arm, finite scores/thresholds, exact confusion counts and independent recomputation. Synthetic tests cover ties, paired group resampling, voting, discordance, invalid predictions and split isolation. Independent audits reconstruct roles and medians from source data, replay saved model predictions and check result arithmetic.

The pre-fit split audit found **40/50 inner missense-validation sets with only two gene groups**, 170 records (149 P/LP and 21 B/LB). The remaining ten sets contain 17–28 groups. Missense inner training contains 871–916 records; outer missense test sets contain 19–39 groups. Aggregate equality of counts does not establish equality of group identities across those 40 sets. Retain the inherited assignments and disclose their limited validation gene diversity; do not claim R90 transportability to new genes.

This uses the same observed ClinVar-derived cohort, not a new blinded or external dataset. No functional or clinical validation is established. Smaller missense-only training may help or hurt performance; training-scope contrasts also differ in sample size and are not isolated biological causal interventions. Fixed-prediction bootstrap intervals omit training uncertainty, rare groups may be unstable, and gene-group resampling does not represent every cross-group homology dependence. Engineering correctness does not resolve these scientific limitations.
