# Detailed consequence stratification (2026-09-11)

The original `coding_nontruncating` group is a coarse coding category, not a missense-only definition. Historical results and filenames are retained, with corrected reader-facing names. No model was retrained in this repair.

| Within the historical 3,121-row subset | P/LP | B/LB | Total |
|---|---:|---:|---:|
| Missense | 926 | 403 | 1,329 |
| Synonymous | 10 | 1,723 | 1,733 |
| Other coding classes combined (inframe/start-stop/other) | 53 | 6 | 59 |

The strict missense SNV subset contains **1,325 variants (924 P/LP, 401 B/LB)**. No AUC, R90, ablation or SHAP result for a model trained specifically on this strict subset is claimed. In particular, the old 0.814 AUC cannot be used as that baseline.

## Definition and provenance

Detailed strata use `p1_consequence_terms`, the JSON consequence mapping for the already-selected P1 target-gene/transcript policy. They do not use the unrestricted all-transcript union. The project-specific hierarchy puts truncating before canonical splice, then other classes. This preserves the legacy grouping and is not the Ensembl severity ordering (which places canonical splice before stop-gained/frameshift); overlapping canonical-splice and truncating flags are retained. Within the remaining coding classes, missense, synonymous, inframe-indel and start/stop-altering conflicts are explicitly marked `transcript_mixed`. Missense with a splice-region term retains a missense primary class and an additional splice-region flag. All selected terms and transcript provenance are retained.

Intronic variants are not called deep-intronic without a splice-distance annotation. Unknown or unrecognized terms remain explicit. The strict missense SNV export requires the existing headline flag, a missense primary class and a syntactically valid single-base REF/ALT key from the existing annotation pipeline (reference concordance is not revalidated here). It preserves the frozen label and gene group.

The additive data and reproducible script are in the Drive project under `runs/stratification_repair_2026-09-11/`: `stratify.py`, `stratified_modeling_data.csv`, `strict_missense_snv.csv`, `cohort_counts.csv`, `summary.json`, and `provenance.json`. Frozen P1–P5 matrices, labels, splits, feature lists, predictions and trained-result files are unchanged. New `p6_` fields are stratification metadata, not model features.

The public repository retains the coarse historical `consequence_class` for replay. A new strict-subset training experiment requires its own manifest and results; it must not overwrite the historical experiment.

## Validation status

Nineteen focused semantic tests passed, including class nesting and overlapping consequences. All 8,056 master rows match the annotation table one-to-one; original columns and order are preserved. There are 7,125 headline rows. Claude Opus and Gemini independently reviewed the supplied code and aggregate counts. Their actionable findings were addressed locally; see `review_adjudication.md`. This is a code/definition review, not clinical validation or a claim of model accuracy.

## Overlap and composition diagnostics

Of the 1,325 missense SNVs, **80 also have selected-transcript splice-region terms**. Thus “strict” means the selected consequence class and SNV representation, not proof of a pure protein mechanism. The separate `strict_missense_snv_no_spliceregion.csv` contains **1,245** variants for a future sensitivity analysis. Headline rows have no `transcript_mixed` or `unknown_other` primary classes; the full 8,056-row table retains 140 unknown/other rows and explicit missing/invalid keys.

In the 3,121-row mixed subset, 10/1,733 synonymous variants are P/LP, versus 979/1,388 other coding variants. A post-hoc rule scoring synonymous as 0 and the remaining coding classes as 1 has a descriptive pooled AUC of **0.8990**. This measures composition confounding; it is not an independently validated predictor or a statistical superiority claim. `evaluation_fold_composition.csv` and `evaluation_fold_rule_auc.csv` audit all 50 frozen outer-fold/seed combinations. Model comparisons must use the same split memberships and a prospectively specified baseline.

The historical **0.8144** is the bootstrap mean pooled AUC for the separately trained 3,121-row coding-nontruncating model (`results/summary.json`, historical key `pooled_missense`). Its pooled R90 recall/precision are 0.9099/0.5251. A different analysis evaluates the headline model within that same biological subgroup: AUC 0.7777, recall 0.530, precision 0.567 (`results/stratified_headline.csv`). These two analyses must not be conflated. Ablation/SHAP in the mixed cohort may reflect proxies for these consequence classes and do not establish missense-specific biology.

Keys are compared exactly without new normalization. Fallback is used only when the primary key is missing, not malformed; a malformed primary key stays explicitly invalid. Two valid but different representations are flagged as conflicts. Original selected SO terms remain available for start/stop distinctions and unsupported future terms.

References: [Ensembl consequence definitions](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html).
