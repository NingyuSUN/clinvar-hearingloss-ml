# Protein language model scores: results (`PLM-ESM-2026-10`)

This report follows the [prespecified protocol](plm_protocol.md) and does not
change it. Every number below comes from
[`results/plm_esm_2026-10/`](../results/plm_esm_2026-10/), whose manifest
hashes the inputs, the code, the sealed predictions and the outputs.

The later [missense-only training follow-up](../experiments/esm_missense_2x2/README.md)
uses newly fitted controls and a common missense validation target. Its results
are reported separately; they do not replace the historical numbers below.

## Findings

- **The primary rule is met, but an increment beyond scored status is not demonstrated.** On the
  3,121 coding non-truncating variants, adding ESM-1v raises the pooled
  gene-held-out AUC from 0.823 to 0.988. The difference is +0.165, with a paired
  gene-group 95% interval of +0.067 to +0.273.
- **A control reproduces almost all of that gain.** The control model receives
  only a 0/1 flag marking variants that *have* an ESM score. It reaches 0.985.
  ESM-1v versus this control (S3) is +0.002, interval −0.001 to +0.007.
  The prespecified verdict is therefore
  `improvement_demonstrated_not_attributable_beyond_scored_status_marker`.
  ORIGINAL18 has no feature that separates missense from synonymous variants.
  Among the 1,793 unscored coding non-truncating variants, 64 are P/LP (3.6%).
  Among the 1,328 scored missense variants, 925 are P/LP (69.7%).
- **The ESM-augmented model ranks scored missense variants better
  (secondary, not multiplicity adjusted).** Training remained mixed.
  AUC rises from 0.967 to 0.984, a difference of +0.017 with interval +0.005 to
  +0.033; 0.3% of bootstrap draws are at or below zero. ESM-2 650M gives a
  similar difference: +0.018, interval +0.008 to +0.033. This historical
  contrast does not isolate a score increment beyond availability: the
  availability-only model reaches 0.978396 on the same scored missense subset,
  versus 0.966674 for ORIGINAL18 and 0.983764 for ORIGINAL18 + ESM-1v.
- **Zero-shot ESM-1v alone is weaker than the existing features.** On scored
  missense its AUC is 0.875 (0.829–0.915), against 0.967 for ORIGINAL18.
  The difference is −0.092, interval −0.142 to −0.048.

The ESM-augmented model improves ranking on held-out scored missense records,
but this mixed-training experiment does not by itself separate score information
from availability information learned during training. Its apparent large effect on
the mixed coding non-truncating subgroup is largely reproduced by the
scored-status control, which also separates consequence classes. This does not
isolate a unique causal mechanism.

## What was run

| Step | Commit |
|---|---|
| Protocol, label-blind mapping and code frozen | `7e7a456` |
| Scoring run 1 failed (GPU memory); fix and execution-log entry | `f7987a8` |
| ESM scores and label-free features committed before fitting | `5f4323d` |
| Out-of-fold predictions sealed (SHA256) before any metric | `a661139` |

The scoring used the five pinned ESM-1v models and ESM-2 650M with
masked-marginal log-likelihood ratios. It ran on Kaggle (two Tesla T4s, torch
2.11.0+cu128, transformers 5.16.1). All 1,189 positions were scored by each
model, and every weight file matched its Hugging Face LFS SHA256.

The fits used the unchanged pipeline: seeds 101–110, five outer gene-group
folds, an inner gene-group holdout, XGBoost 3.4.1 and the R90 rule. There were
600 fold fits (2 tasks × 6 feature sets × 50).

## Denominators

| Evaluation set | Variants | P/LP | B/LB | Gene groups |
|---|---:|---:|---:|---:|
| Coding non-truncating (Task A) | 3,121 | 989 | 2,132 | 148 |
| Scored missense (subset of Task A) | 1,328 | 925 | 403 | 135 |
| Headline (Task B) | 7,125 | 3,576 | 3,549 | 167 |

1,329 frozen-cache missense variants were candidates. One, a multi-nucleotide
change with two altered residues, is unscored by rule. All 1,328 scored
variants are in the coding non-truncating subset.

## Prespecified contrasts

Delta = AUC(b) − AUC(a) on the same variants, from seed-averaged OOF scores.
Intervals are percentile 95% from 2,000 paired whole-gene-group draws (seed 0;
all draws valid). "LB > 0" is the superiority rule.

| ID | Tier | Evaluation set | a → b | AUC a | AUC b | Delta [95% interval] | LB > 0 |
|---|---|---|---|---:|---:|---|---|
| **P1** | **primary** | coding non-truncating | ORIGINAL18 → +ESM-1v | 0.823 | 0.988 | **+0.165 [+0.067, +0.273]** | yes |
| S1 | secondary | headline | ORIGINAL18 → +ESM-1v | 0.931 | 0.992 | +0.061 [+0.022, +0.095] | yes |
| S2 | secondary | scored missense | ORIGINAL18 → +ESM-1v | 0.967 | 0.984 | +0.017 [+0.005, +0.033] | yes |
| S3 | secondary | coding non-truncating | +scored flag → +ESM-1v | 0.985 | 0.988 | +0.002 [−0.001, +0.007] | **no** |
| D1 | descriptive | coding non-truncating | ORIGINAL18 → +ESM-2 | 0.823 | 0.987 | +0.164 [+0.066, +0.271] | — |
| D2 | descriptive | headline | ORIGINAL18 → +ESM-2 | 0.931 | 0.991 | +0.060 [+0.021, +0.093] | — |
| D3 | descriptive | scored missense | ORIGINAL18 → +ESM-2 | 0.967 | 0.985 | +0.018 [+0.008, +0.033] | — |
| D4 | descriptive | coding non-truncating | ORIGINAL18 → +scored flag | 0.823 | 0.985 | +0.163 [+0.064, +0.270] | — |
| D5 | descriptive | coding non-truncating | ORIGINAL18 → ESM-1v only | 0.823 | 0.929 | +0.106 [−0.003, +0.209] | — |
| D6 | descriptive | scored missense | ORIGINAL18 → ESM-1v only | 0.967 | 0.852 | −0.115 [−0.173, −0.066] | — |
| D7 | descriptive | headline | +scored flag → +ESM-1v | 0.991 | 0.992 | +0.001 [−0.000, +0.002] | — |

The rule reads the prespecified wording table as follows. P1's lower bound is
above zero, and S3's is not, so an increment beyond the marker is not
demonstrated. S1 has the same structure as P1: D7, the headline equivalent of
S3, is +0.001. S2 and D3 have lower bounds above zero, but a constant scored
flag in the evaluation subset does not exclude availability effects learned
from mixed training. S1–S3 are secondary and not multiplicity adjusted.
The later C/D contrast restricts both training and evaluation to scored
missense records and directly addresses this limitation.

Matched-fold mean deltas (descriptive; folds share training data):

| Contrast | Mean delta (50 folds) | Range of seed means |
|---|---:|---|
| P1 | +0.128 | +0.117 to +0.140 |
| S1 | +0.042 | +0.037 to +0.046 |
| S3 | +0.001 | +0.000 to +0.002 |

## Model metrics

Pooled AUC, plus R90 recall and precision from the fold-threshold vote. All
intervals are conditional gene-group bootstrap intervals.

| Evaluation set | Model | AUC [95% interval] | R90 recall | R90 precision |
|---|---|---|---:|---:|
| coding non-truncating | ORIGINAL18 | 0.823 [0.712, 0.924] | 0.908 | 0.531 |
| | ESM-1v only | 0.929 [0.896, 0.950] | 0.760 | 0.861 |
| | ORIGINAL18 + ESM-1v | 0.988 [0.982, 0.992] | 0.843 | 0.986 |
| | ORIGINAL18 + scored flag | 0.985 [0.977, 0.991] | 0.877 | 0.975 |
| | ESM-2 only | 0.930 [0.900, 0.950] | 0.799 | 0.856 |
| | ORIGINAL18 + ESM-2 | 0.987 [0.978, 0.992] | 0.814 | 0.985 |
| scored missense | ORIGINAL18 | 0.967 [0.948, 0.980] | 0.915 | 0.969 |
| | ORIGINAL18 + ESM-1v | 0.984 [0.974, 0.991] | 0.889 | 0.987 |
| | ORIGINAL18 + ESM-2 | 0.985 [0.976, 0.992] | 0.858 | 0.986 |
| | ESM-1v only (trained) | 0.852 [0.798, 0.899] | 0.813 | 0.861 |
| headline | ORIGINAL18 | 0.931 [0.894, 0.974] | 0.846 | 0.852 |
| | ORIGINAL18 + ESM-1v | 0.992 [0.986, 0.996] | 0.895 | 0.996 |
| | ESM-1v only | 0.571 [0.446, 0.706] | 0.628 | 0.529 |

The full table, with recall and precision intervals and every model on every
set, is in `model_metrics.csv`.

Adding ESM-1v moves the operating point. On coding non-truncating variants,
precision rises from 0.531 to 0.986 while recall falls from 0.908 to 0.843,
below the 0.90 validation target. On the headline cohort, PLM-only models are
near chance (0.571): only 19% of headline variants have a score.

## Zero-shot scores on scored missense

These scores use no training: they rank by −LLR. There are 1,328 variants,
925 of them P/LP, in 135 gene groups.

| Score | AUC [95% interval] |
|---|---|
| ESM-1v five-model mean | 0.875 [0.829, 0.915] |
| ESM-2 650M | 0.879 [0.840, 0.907] |
| ESM-1v single members | 0.856–0.875 |

Zero-shot ESM-1v minus ORIGINAL18 OOF on the same variants is −0.092
[−0.142, −0.048].

## ORIGINAL18 reproduction check

| Cohort | Published AUC | Re-run AUC | Published TP/FP/TN/FN | Re-run TP/FP/TN/FN |
|---|---:|---:|---|---|
| coding non-truncating | 0.8148 | 0.8226 | 900/813/1,319/89 | 898/794/1,338/91 |
| headline | 0.9301 | 0.9309 | 3,022/513/3,036/554 | 3,025/525/3,024/551 |

The re-run does not exactly reproduce the published predictions; this is
reported, as the protocol requires. This repository's pipeline is
deterministic in the same environment:
- 1, 8 and 10 XGBoost threads give identical predictions;
- an independent re-fit of seed 101 matches the sealed predictions to within
  3 × 10⁻⁸ (CSV rounding).

The published predictions come from the archived 2026-09-10 run, which predates
this repository's package (2026-09-11). Its code is not public, so the source
of the 0.008 difference cannot be identified here. All contrasts compare models
fitted by the same code in the same run, as the protocol specifies.

## Interpretation limits

- **Availability control.** Score availability is associated with missense
  status in this dataset. On mixed subgroups a scored-status control reproduces
  most of the gain (S3, D7), while a score increment beyond that control is not
  demonstrated. ORIGINAL18 lacks basic consequence resolution within coding
  non-truncating variants. This comparison does not prove that all improvement
  has a single causal explanation.
- **Label circularity.** ClinVar submitters may use computational evidence
  (ACMG PP3/BP4) and allele frequency (BA1/BS1) when assigning labels. Part of
  any predictor's apparent within-missense performance can reflect evidence
  that helped create the label. This applies to ESM scores and to ORIGINAL18's
  conservation and frequency features.
- **Conditional intervals.** Intervals resample gene groups with predictions
  held fixed. They exclude refitting, model selection and full cross-validation
  dependence ([statistical revision](statistical_revision.md)). Secondary
  results are unadjusted.
- **Small absolute effect at a high baseline.** Within missense, ORIGINAL18
  already reaches 0.967. The +0.017 gain is about half of its remaining
  ranking error (1 − AUC falls from 0.033 to 0.016). It does not establish
  clinical utility.
- **Coverage and windows.** 47 of 135 scored proteins exceed 1,022 residues
  and are scored in a centred window. Non-missense coding variants receive no
  score by design.
- **Not done in this historical run.** No calibration, no external validation,
  no model trained only on missense variants, and no comparison with AlphaMissense
  or other supervised predictors. Missense-only training was subsequently
  completed in the linked 2026-10-21 follow-up; the other gaps remain.

## Reproduce

```bash
python scripts/plm_map_variants.py --out data/plm            # label-blind, Ensembl REST
python kaggle/push_esm.py push                                # private Kaggle GPU kernel
python kaggle/push_esm.py fetch --dest data/plm/esm_scores
python scripts/plm_build_features.py
HLPATH_NTHREAD=10 python scripts/run_plm_evaluation.py --out runs/plm_esm_2026-10
python scripts/analyze_plm.py --runs runs/plm_esm_2026-10 --out results/plm_esm_2026-10
```

Ensembl REST serves only the current release, which is 116 at the time of
writing; after it moves on, a re-run of the mapping may differ. The committed
mapping files are pinned by SHA256, so re-run the mapping only to audit them.
