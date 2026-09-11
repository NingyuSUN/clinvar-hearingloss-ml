# Results

> **Stratification correction (2026-09-11):** The historical 3,121-row `coding_nontruncating` subset includes 1,329 missense, 1,733 synonymous and 59 other coding variants. Its AUC 0.814, R90, ablation and SHAP results are **not strict-missense results**. Historical filenames containing `missense` are retained for reproducibility. See [stratification details](stratification.md).

Gene-held-out grouped 5-fold cross-validation, 10 frozen seeds. XGBoost.
Raw tables are in `results/`; `run_manifest.json` records the XGBoost version and
the matrix hash.

Three ways to summarise AUC:

| | what it measures |
|---|---|
| per-fold mean ± SD | average fold-model AUC on its own held-out genes — SD is **not** a CI |
| pooled + bootstrap 95% CI | the seed-ensemble ranking over the whole cohort, resampling variants |
| gene-group-macro | unweighted over gene groups; dominated by the many small single-class groups |

## Headline cohort (N = 7,125)

| Full model (18 features) | value |
|---|---|
| AUC, per-fold mean ± SD | 0.950 ± 0.016 |
| AUC, pooled + bootstrap 95% CI | **0.930 [0.924, 0.936]** |
| AUC, gene-group-macro | 0.978 |
| R90 recall (pooled) | 0.845 [0.833, 0.858] |
| R90 precision (pooled) | 0.855 [0.843, 0.866] |

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

| Step | Headline Δ (t) | coding non-truncating Δ (t) |
|---|---|---|
| Base → + consequence | +0.128 (23.7) ✓ | +0.006 (3.9) ✓ |
| + consequence → + constraint | +0.001 (0.2) | +0.025 (1.5) |
| + constraint → + frequency | +0.020 (3.1) ✓ | +0.146 (8.0) ✓ |
| + frequency → Full | +0.004 (4.4) ✓ | +0.002 (0.7) |

✓ = significant at p < 0.05. Note that + domain adds only 0.004 AUC on the full
cohort but is a significant *paired* difference — comparing the mean delta to the
fold SD (0.016) would call it noise.

### Leave-one-group-out from the full model

| Remove | Headline Δ (t) | coding non-truncating Δ (t) |
|---|---|---|
| consequence flags | −0.085 (−11.6) ✓ | +0.001 (0.2) |
| gene constraint | **+0.011 (+4.4)** ✓ | **−0.054 (−4.0)** ✓ |
| conservation | −0.012 (−4.0) ✓ | −0.040 (−4.6) ✓ |
| allele frequency | −0.010 (−1.8) | −0.136 (−6.8) ✓ |
| protein domain | −0.004 (−3.6) ✓ | −0.004 (−1.1) |
| allele shape | +0.001 (0.9) | +0.003 (0.9) |

Gene constraint slightly *hurts* the near-trivial full-cohort AUC (it lets the
model over-fit gene identity) but is the second most important feature for
coding non-truncating. It is kept.

## coding non-truncating subset (coding non-truncating, N = 3,121, prevalence 0.32)

| Full model | value |
|---|---|
| AUC, per-fold mean ± SD | 0.856 ± 0.038 |
| AUC, pooled + bootstrap 95% CI | **0.814 [0.800, 0.829]** |
| AUC, gene-group-macro | 0.941 |
| R90 recall (pooled) | 0.910 [0.891, 0.927] |
| R90 precision (pooled) | 0.525 [0.502, 0.548] |

Frequency-only baseline: 0.772. The full model adds ~0.08 over frequency alone,
from conservation, gene constraint and protein domain.

## Consequence-stratified recall (headline model, global R90 threshold)

| Consequence class | N | Prevalence | AUC | R90 recall | R90 precision |
|---|---:|---:|---:|---:|---:|
| truncating | 1,952 | 0.999 | — | 0.997 | 0.999 |
| canonical splice | 563 | 0.996 | — | 0.966 | 0.998 |
| coding non-truncating | 3,121 | 0.317 | 0.778 | **0.530** | 0.567 |
| non-coding / other | 1,489 | 0.051 | 0.833 | 0.158 | 0.098 |

(`AUC` blank for the single-class LoF strata. Recall is the fraction of pathogenic
variants in each stratum recovered at the one global R90 threshold.)

## Diagnostics

### Gene-held-out vs random row split (full model, pooled AUC)

| Split | AUC |
|---|---:|
| gene-held-out | 0.927 |
| random row split | 0.996 |
| gap | **0.069** |

### Review-status sensitivity (full model, pooled + bootstrap AUC)

| Cohort | N | AUC [95% CI] |
|---|---:|---|
| all-tier (0-star included) | 7,736 | 0.930 [0.924, 0.936] |
| ≥ 1 star (headline) | 7,125 | 0.930 [0.924, 0.936] |
| ≥ 2 star | 3,554 | 0.967 [0.962, 0.971] |

AUC is stable across the review-status cuts and only tightens as the labels get
cleaner.

## Feature attribution

`docs/shap.md` has the full write-up. Headline result: TreeSHAP attribution on
the fitted models ranks the coding non-truncating feature groups in the same order as the
ablation above (frequency > gene constraint > conservation > domain, Spearman
0.94) — two independent methods agreeing on which features carry the coding non-truncating
signal.
