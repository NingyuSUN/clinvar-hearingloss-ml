# Hearing-loss variant pathogenicity — gene-held-out evaluation

[![Repository checks](https://github.com/NingyuSUN/clinvar-hearingloss-ml/actions/workflows/ci.yml/badge.svg)](https://github.com/NingyuSUN/clinvar-hearingloss-ml/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A machine-learning study of ClinVar hearing-loss variants. The question is not
"how high can the AUC go" but **how much of a pathogenicity model's apparent
performance survives when the test genes are genuinely unseen**, and where the
signal actually comes from.

> Research / portfolio project. Not a clinical diagnostic tool. See
> [`MODEL_CARD.md`](MODEL_CARD.md) for intended/out-of-scope use and
> [`PROJECT_STATUS.md`](PROJECT_STATUS.md) for what is and isn't finished.

## What this project does

- Builds a per-variant feature table from ClinVar, Ensembl VEP consequence terms,
  gnomAD allele frequency and gene constraint, protein-domain membership, and an
  evolutionary conservation score.
- **Freezes** the label rule, the cohort, the gene-group split, the operating
  point, and the feature list *before any model is trained* (`docs/methods.md`).
- Evaluates an XGBoost classifier with **grouped 5-fold cross-validation on gene
  groups** (connected components of gene co-occurrence), a nested gene-group
  hold-out for iteration and threshold selection, and 10 frozen seeds.
- Reports an additive feature ablation with **paired per-fold tests**, a
  **patient-level bootstrap CI** on the pooled out-of-fold predictions, and a
  **consequence-stratified** breakdown.

## Headline results

Headline cohort: 7,125 variants (feature-complete, ClinVar review status ≥ 1 star),
167 gene groups, roughly balanced labels.

| Model (18 features) | AUC | 95% CI |
|---|---:|---:|
| Full model, gene-held-out | **0.930** | 0.924 – 0.936 |
| Full model, ordinary random split | 0.996 | 0.995 – 0.997 |

A random split lets the model recognise the gene and inflates AUC by ~0.07. Every
number below is gene-held-out.

### The cohort is mostly decided by consequence type

In a hearing-loss gene panel, loss-of-function is near-deterministic for
pathogenicity (ACMG `PVS1`). That shows up as near-perfect separation before any
model runs:

| Consequence class | Variants | Pathogenic |
|---|---:|---:|
| Truncating (frameshift / stop) | 1,952 | 99.9 % |
| Canonical splice | 563 | 99.6 % |
| Non-coding / other | 1,489 | 5.1 % |
| **Coding non-truncating (missense)** | **3,121** | **31.7 %** |

Roughly 35 % of the cohort is one class. A model using only the four consequence
flags already scores **AUC 0.85** on the full cohort. The full-cohort 0.93 is
mostly loss-of-function identification.

### On the genuine problem — missense — the model is weaker

| Missense subset (N = 3,121) | AUC | 95% CI | R90 precision |
|---|---:|---:|---:|
| Full model | **0.814** | 0.800 – 0.829 | 0.53 |
| Frequency only | 0.772 | — | — |

Allele frequency is the dominant feature here (removing it costs 0.136 AUC), and
frequency also feeds ClinVar's own benign calls through `BA1` / `BS1`, so part of
that signal is the model re-deriving the labelling rule.

### The R90 operating point hides the missense gap

The threshold is chosen for recall ≥ 0.90 on the validation set. At that single
global threshold:

| Consequence class | Recall |
|---|---:|
| Truncating | 99.7 % |
| Canonical splice | 96.6 % |
| **Missense** | **53.0 %** |
| Non-coding / other | 15.8 % |

A headline "90 % recall" is carried by the trivial classes; about one in two
missense pathogenic variants is missed.

### Feature ablation (paired per-fold, full cohort → missense)

| Add to Base | Full cohort ΔAUC | Missense ΔAUC |
|---|---:|---:|
| consequence flags | **+0.128** (t 23.7) | +0.006 (t 3.9) |
| gene constraint | +0.001 (ns) | +0.025 (ns) |
| allele frequency | +0.020 (t 3.1) | **+0.146** (t 8.0) |
| protein domain | +0.004 (t 4.4) | +0.002 (ns) |

Consequence type carries the full cohort; frequency carries missense.
Under leave-one-group-out from the full model, **gene constraint is the second
most important feature for missense** (−0.054, t −4.0) even though it adds nothing
on the full cohort — so it is kept.

Full tables: `docs/results.md` and `results/`.

### Feature attribution agrees with the ablation

TreeSHAP (exact, out-of-fold, on the same frozen models) ranks the missense
feature groups the same way the ablation does — frequency > gene constraint >
conservation > domain (Spearman rank correlation 0.94). Conservation's
attribution shows a sharp threshold around a score of ~1–2, not a smooth
gradient. See `docs/shap.md`.

## Repository layout

```text
LICENSE  CITATION.cff  MODEL_CARD.md  PROJECT_STATUS.md  CHANGELOG.md
.github/workflows/ci.yml   CI: tests + public-artifact validation, no retraining
Makefile                    test / validate / syntax targets
tools/validate_public_artifacts.py   dependency-free checks on results/ + docs/
tests/                       pytest wrapper around the validator
src/hlpath/
  protocol.py       frozen label / cohort / gene-group / feature definitions
  pipeline.py       balanced grouped CV, nested inner hold-out, training, R90
  metrics.py        R90 threshold (with one-class guard), bootstrap CI, paired tests
  analysis.py       ladder / stratified / split-gap tables
  shap_analysis.py  TreeSHAP attribution on the frozen model
  features.py       load the modelling matrix
  config.py         paths, seeds, XGBoost params
scripts/
  run_evaluation.py    the full cross-validation run
  analyze_results.py   the reported tables
  run_shap.py          TreeSHAP tables + figures
  build_matrix.py      regenerate the matrix from a raw annotated feature table
data/
  modeling_matrix.csv.gz   per-variant features + label + gene group + cohort flags
docs/
  methods.md  data.md  results.md  shap.md  limitations.md
results/
  ladder_*.csv  stratified_headline.csv  split_gap ...  run_manifest.json
  shap/  importance / group_importance / by_consequence / conservation_bands + figures/
```

## Reproduce

```bash
conda env create -f environment.yml
conda activate hearingloss-variants
pip install -e .

python scripts/run_evaluation.py --out results/
python scripts/analyze_results.py --dir results/
```

To check the repository itself (fast, no retraining, no data required beyond
what's already committed):

```bash
make test      # pytest: JSON/CSV well-formedness, no CRLF regressions
make validate  # same checks, standalone script (what CI runs)
make syntax    # ast-parse every tracked .py file
```

`run_evaluation.py` reads `data/modeling_matrix.csv.gz` (override with
`HLPATH_DATA`). The fold assignments are a deterministic function of the 10 frozen
seeds in `src/hlpath/config.py`. To rebuild the matrix from raw annotations, see
`docs/data.md` and `scripts/build_matrix.py`.

## Methods in brief

- **Label** — ClinVar aggregate `ClinicalSignificance`: P/LP → 1, B/LB → 0,
  "Likely" kept, VUS / Conflicting excluded. (Not `ClinSigSimple`.)
- **Cohort** — feature-complete gate (reliable consequence, present gene
  constraint, definitive frequency status, coding call made, conservation present,
  single gene, label present) ∩ review status ≥ 1 star. All-tier and ≥ 2-star
  cohorts are reported as sensitivity analyses.
- **Split** — gene groups are connected components of `GeneSymbol` co-occurrence;
  a whole component is held out together. Outer grouped 5-fold, folds balanced on
  row count and label prevalence. Inner: a balanced gene-group hold-out for
  iteration selection and the R90 threshold.
- **Model** — XGBoost (`eta` 0.03, depth 5), trained to a 700-round cap with no
  in-loop early stopping; the best iteration is the argmin of the
  inner-validation log-loss curve. The R90 threshold is frozen on the inner
  validation and applied to the outer test with the same fitted model.
- **Statistics** — the mean ± SD across overlapping folds is not a confidence
  interval; a patient-level bootstrap on the pooled predictions is used instead,
  and ladder steps use paired per-fold differences.

## Limitations

See `docs/limitations.md`. In short: excluding VUS / Conflicting makes an easier,
selected task; consequence type near-determines the label for a third of the
cohort; allele frequency is partly circular with the ClinVar labels; the
conservation feature's provenance and the ClinVar phenotype-condition mapping are
not fully audited.

## Data sources

- NCBI ClinVar `variant_summary.txt.gz` (public domain).
- Ensembl VEP (REST, GRCh38) — consequence terms, canonical / MANE transcripts,
  protein-domain overlaps.
- gnomAD v4.1.1 — allele frequency (PASS); gnomAD v2.1.1 — gene constraint
  (`oe_lof` / `oe_mis` upper bound).
- Ensembl comparative-genomics conservation score.
