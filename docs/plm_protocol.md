# Prespecified protocol: protein language model scores (`PLM-ESM-2026-10`)

**Status: frozen.** This document was committed before any ESM score was
computed with the protocol's models, before any model containing a PLM feature
was fitted and before ORIGINAL18 was re-run. Work done before the freeze, all
of it label-blind:

- the variant-to-protein mapping (`scripts/plm_map_variants.py`; the label
  column is not read, and `data/plm/mapping_manifest.json` records the counts below);
- an engineering smoke test of the scoring kernel with the 8M-parameter
  `facebook/esm2_t6_8M_UR50D` model on 45 positions. No variant labels were
  involved, and the outputs were discarded.

`src/hlpath/plm.py` and `src/hlpath/compare.py` are the executable form of
this protocol. `scripts/analyze_plm.py` hard-codes the contrasts listed below.

## Question

Do zero-shot protein language model scores improve the gene-held-out
discrimination of the existing 18-feature model (ORIGINAL18 = `L6_full`)?
This is the setting where published work suggests PLMs help: missense
pathogenicity.

## Fixed denominators (published cohort, unchanged)

| Set | Variants | Gene groups | Source |
|---|---:|---:|---|
| Headline cohort | 7,125 | 167 | `results/statistics_20260921` |
| Coding non-truncating subset | 3,121 | 148 | same |
| Scored missense (subset of coding non-truncating) | 1,328 | reported at analysis | `data/plm/mapping_manifest.json` |

Composition of the coding non-truncating subset by frozen cache variant type:
1,733 synonymous, 1,329 missense, 28 in-frame indel, 25 start/stop-altering and
6 other coding. **Only missense variants can receive a score, so 57% of the
primary evaluation set is unscored by construction.** Label prevalence within the
scored missense set is not examined before fitting.

## Feature definition

**Transcript.** The frozen cached transcript of each variant
(`predict/bundle/annotation_cache.jsonl`, field `transcript_ids`; 7,122 of
7,125 are MANE Select). No new transcript choice is made.

**Candidates.** Headline variants whose frozen `consequence_class` is
`coding_nontruncating` or `canonical_splice`. Truncating and non-coding classes
cannot carry a missense change on the cached transcript and are unscored.

**Annotation.** Ensembl VEP REST (`POST /vep/human/region`, Ensembl release
116, GRCh38) with `protein=1, mane=1, canonical=1, transcript_version=1`. The
entry for the cached transcript's stable ID is used.

**Scored if all of:**
1. VEP `consequence_terms` include `missense_variant`;
2. `amino_acids` is a single-residue substitution `X/Y` with
   `protein_start == protein_end`;
3. X and Y are distinct standard residues (`ACDEFGHIKLMNPQRSTVWY`);
4. the Ensembl release-116 translation of the cached transcript
   (`POST /sequence/id`, `type=protein`) has X at `protein_start`.

Otherwise the variant gets a recorded status (`not_missense`,
`not_single_residue_substitution`, `nonstandard_residue`,
`reference_residue_mismatch`, …) and no score. Mapping outcome before
the freeze: 1,328 scored (all frozen-cache missense, all coding non-truncating),
1 missense multi-nucleotide change with two residues changed (unscored),
0 reference mismatches, 0 transcript-version changes; 1,189 unique positions on
135 proteins, 47 of them longer than 1,022 residues.

**Models** (Hugging Face, pinned commit):

| Role | Model | Revision |
|---|---|---|
| primary (5-model ensemble) | `facebook/esm1v_t33_650M_UR90S_1` … `_5` | `8bfdb18…`, `3c1e9e6…`, `0b00fd1…`, `443968f…`, `fb2e51c…` |
| sensitivity | `facebook/esm2_t33_650M_UR50D` | `08e4846…` |

Full hashes are in `src/hlpath/plm.py::MODEL_REVISIONS`. The kernel refuses
other revisions and checks downloaded weights against the hub LFS SHA256.

**Window.** Proteins of at most 1,022 residues are scored whole. Longer proteins
use a 1,022-residue window centred on the variant (`start = p − 1 − 511`),
clamped to the protein ends. The tokens are `<cls>` + window + `<eos>`.

**Score (masked marginal).** Mask the variant residue only. Run one fp32 forward
pass, take the log-softmax over the full vocabulary at the masked token, then
compute `LLR = log p(alt) − log p(ref)`.
- `esm1v_llr` = mean of the five ESM-1v LLRs (missing unless all five exist).
- `esm2_llr` = the ESM-2 650M LLR.
- `plm_scored` = 1 when `esm1v_llr` exists, else 0.

No further transformation is applied. Trees are invariant to monotone
transforms. Zero-shot AUCs use `−LLR`, a fixed direction, because a lower LLR
is predicted as more damaging.

**Unscored variants.** Their scores stay `NaN` and are imputed **inside every
split with the inner-train median**, the existing rule for continuous features.
Synonymous variants are *not* set to LLR 0. Imputing a single constant still
makes "unscored" visible to the trees, which is consequence information that
ORIGINAL18 lacks. Contrasts S2 and S3 below exist to separate this marker from
the score itself.

## Feature sets (models)

| Name | Features |
|---|---|
| `ORIGINAL18` | the 18 existing features (identical to `L6_full`) |
| `ESM1V_ONLY` | `esm1v_llr` (the PLM-score-only model) |
| `ORIGINAL18_ESM1V` | 18 features + `esm1v_llr` |
| `ORIGINAL18_SCOREDFLAG` | 18 features + `plm_scored` (scored-status control) |
| `ESM2_ONLY` | `esm2_llr` (sensitivity) |
| `ORIGINAL18_ESM2` | 18 features + `esm2_llr` (sensitivity) |

XGBoost parameters, the 700-round cap, post hoc argmin-log-loss iteration
selection, the R90 threshold and the imputation rule are unchanged
([methods](methods.md)).

## Evaluation (unchanged splits and seeds)

The repository's `run_experiment` is used with seeds 101–110, five outer
gene-group folds and the same balanced fold and inner-holdout functions. Fold
assignment depends only on group sizes and prevalence, so the folds are the
published ones.

- **Task A** (`plm_codingnt`): train and evaluate within the 3,121-variant
  coding non-truncating subset, as the published 0.815 subset model did.
- **Task B** (`plm_headline`): the 7,125-variant headline cohort.
- **Scored-missense evaluation set**: Task A variants with `plm_scored = 1`,
  evaluated with Task A predictions (no separate training).

Aggregation follows the 2026-09-21 revision: per-variant seed-averaged OOF
scores (`analysis._pooled`), with fold thresholds applied before voting. The
headline thread count is fixed at `HLPATH_NTHREAD=10`, matching the original
run.

Execution order:
1. commit this protocol;
2. score on Kaggle;
3. build features and commit their hashes;
4. fit, then seal predictions in `prediction_seal.json`;
5. analyse;
6. report.

No metric is printed during fitting.

## Endpoints

AUC is the variant-weighted pooled ROC-AUC of the seed-averaged OOF scores.
Each delta is `AUC(b) − AUC(a)` on identical variants.

| ID | Tier | Task / evaluation set | a → b |
|---|---|---|---|
| **P1** | **primary** | A / all coding non-truncating (3,121) | ORIGINAL18 → ORIGINAL18_ESM1V |
| S1 | secondary | B / headline (7,125) | ORIGINAL18 → ORIGINAL18_ESM1V |
| S2 | secondary | A / scored missense (1,328) | ORIGINAL18 → ORIGINAL18_ESM1V |
| S3 | secondary | A / all coding non-truncating | ORIGINAL18_SCOREDFLAG → ORIGINAL18_ESM1V |
| D1–D3 | descriptive | A all, B, A scored missense | ORIGINAL18 → ORIGINAL18_ESM2 |
| D4 | descriptive | A / all | ORIGINAL18 → ORIGINAL18_SCOREDFLAG |
| D5, D6 | descriptive | A all, A scored missense | ORIGINAL18 → ESM1V_ONLY |
| D7 | descriptive | B / headline | ORIGINAL18_SCOREDFLAG → ORIGINAL18_ESM1V |

Also reported descriptively:
- the AUC, R90 recall and R90 precision of every model on every evaluation set;
- zero-shot AUCs (`−esm1v_llr`, `−esm2_llr` and each ESM-1v member) on scored
  missense, plus zero-shot ESM-1v versus ORIGINAL18 there;
- matched-fold mean deltas (`metrics.paired_delta`, no significance).

## Uncertainty and decision rule

**Paired whole-gene-group bootstrap** (`compare.paired_group_bootstrap_auc`):
2,000 draws, seed 0. Each draw samples G gene groups (`p4_gene_group`) of the
evaluation set with replacement and keeps all of their variants. Both models'
AUCs and their difference are computed on the same draw, and the result is a
percentile 95% interval. A draw lacking a label class is undefined and counted;
an interval is withheld if fewer than 90% of draws are valid. Intervals are
conditional on the fitted models and fixed predictions. They exclude refitting
and full CV dependence, as in the [statistical revision](statistical_revision.md).

**Superiority rule:** a contrast shows superiority only if the 2.5th percentile
of its paired delta is **above zero**. Only P1 decides the headline claim. S1–S3
use the same rule, are reported as secondary and are not multiplicity adjusted.
They cannot rescue a negative P1.

Prespecified wording:

| P1 lower bound | S3 lower bound | Statement |
|---|---|---|
| ≤ 0 | any | "Improvement not demonstrated." This is never stated as equivalence, whatever the point estimate. |
| > 0 | > 0 | "ESM-1v adds discrimination beyond ORIGINAL18 and beyond the scored-status marker." |
| > 0 | ≤ 0 | "Improvement demonstrated, but not attributable to the score beyond the missense-status marker." |

S2 is reported in every case. It is the cleanest test of the score itself,
because `plm_scored` is constant within it.

**ORIGINAL18 reproduction check.** The re-run ORIGINAL18 is compared with the
published empirical AUCs (0.8147647 subset, 0.9300532 headline) and confusion
matrices. A mismatch is reported with its size. It does not stop the analysis.
All contrasts use the re-run ORIGINAL18 from the same code and environment,
never the published numbers.

## Not allowed after the freeze

- Changing hyperparameters, windows, models, revisions, imputation, feature
  transforms, contrasts, evaluation sets, bootstrap settings or the rule after
  seeing any fold result. Any extra analysis is labelled *post hoc exploratory*.
- Promoting ESM-2 to primary or replacing P1 by a secondary contrast.
- Any technical scoring failure is recorded as an amendment before fitting,
  with the affected variants left unscored, never after.

## Execution log (amendments after the freeze)

1. **Scoring run 1 (Kaggle kernel version 1, 2026-10-20) failed.** All five
   ESM-1v models scored the 1,189 positions (about 7.3 min each; weights matched
   the hub LFS SHA256). Loading ESM-2 then ran out of GPU memory, because the
   earlier models had not been freed: HF modules hold reference cycles that
   `del` does not break. The kernel writes its output only at the end, so no
   scores were saved, and none were seen or joined to labels.

   Engineering fix:
   - explicit `gc.collect()` before `empty_cache()`;
   - per-model partial output files;
   - one worker per available GPU pulling models from a queue.

   The scoring definition, models, revisions, windows and inputs are unchanged.
2. **Scoring run 2 (kernel version 2) succeeded.** All six models scored all
   1,189 positions on two T4s. Scores and label-free features were committed
   (`5f4323d`) before fitting.
3. **Fits sealed.** 600 fold fits finished in 18 minutes. The SHA256 of all
   14 prediction files was committed (`a661139`) before any metric was computed.
4. **ORIGINAL18 reproduction check: not exact.** The subset AUC is 0.8226
   against the published 0.8148; the headline AUC is 0.9309 against 0.9301.
   In a post hoc diagnostic, 1, 8 and 10 threads gave identical predictions, and
   a re-fit matched the sealed predictions. The published run used the archived
   pre-package implementation. As prespecified, all contrasts use the re-run. See
   [results](plm_results.md).

## Software

Python 3.12, xgboost 3.4.1, numpy 2.5.3, pandas 3.0.1 (matching the original
run's xgboost/numpy). The Kaggle kernel records its own torch, transformers and
CUDA versions, the GPU, and the weight, input and output SHA256
(`data/plm/esm_scores/esm_scoring_manifest.json`).
