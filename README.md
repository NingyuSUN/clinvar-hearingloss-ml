# Hearing Loss Variant Pathogenicity Prediction

Machine learning project for predicting pathogenicity of hearing-loss-related variants using ClinVar-derived features, evolutionary conservation, strict gene-based validation, and SHAP explainability.

> **Status:** ongoing research / portfolio project. This repository is not a clinical diagnostic tool.

## Highlights

- Diagnosed why ordinary random splits can overestimate variant pathogenicity models.
- Re-evaluated performance using **gene-based splits**, where genes in the test set are unseen during training.
- Added **evolutionary conservation** as a cross-gene biological feature.
- Evaluated models with both AUC and a clinically motivated **R90** operating point: recall >= 0.90.
- Used SHAP to interpret global feature importance, conservation gray zones, subgroup behavior, and representative model errors.

## Results

### 1. Random split vs gene-based split

Random split produced stronger apparent performance than gene-based split, suggesting that random split may partly reward gene-associated shortcut learning.

| Evaluation setting | AUC |
|---|---:|
| Random split | 0.878 |
| Gene-based split | 0.723 |

### 2. Effect of conservation under strict gene-based evaluation

Adding `ensembl_conservation` improved unseen-gene generalization and reduced false positives under the R90 setting.

| Model | AUC | R90 FN | R90 FP | R90 Recall | R90 Precision |
|---|---:|---:|---:|---:|---:|
| Baseline, no conservation | 0.653 | 1 | 664 | 0.998 | 0.449 |
| + Conservation | 0.812 | 54 | 394 | 0.900 | 0.553 |
| Delta | +0.159 | +53 | -270 | -0.098 | +0.104 |

`R90` means the threshold is selected under the constraint that recall is at least 0.90, then FN, FP, and precision are evaluated at that threshold.

### 3. SHAP interpretation

The SHAP analyses suggest a hierarchical decision pattern:

```text
Strong conservation signal       -> conservation-dominated prediction
Ambiguous conservation gray zone -> structural features refine prediction
```

In the approximate conservation gray zone, `ensembl_conservation` between -1 and 1, model decisions become more dependent on variant structure features such as deletion, stop-gain, splice-site signal, and variant length.

## Repository structure

```text
.
├── README.md
├── requirements.txt
├── environment.yml
├── pyproject.toml
├── src/hlpath/
│   ├── clinvar.py           # ClinVar download, filtering, labeling
│   ├── data.py              # loading and feature engineering
│   ├── metrics.py           # AUC, confusion matrix, R90 metrics
│   ├── modeling.py          # XGBoost training and evaluation
│   ├── heterogeneity.py     # subgroup and conservation-bin analyses
│   └── shap_analysis.py     # SHAP plotting utilities
├── scripts/
│   ├── 00_download_prepare_clinvar.py
│   ├── 01_initial_random_forest_baseline.py
│   ├── 01_diagnose_random_vs_gene_split.py
│   ├── 02_evaluate_conservation.py
│   ├── 03_heterogeneity_analysis.py
│   └── 04_shap_explainability.py
├── docs/
│   ├── data_source.md
│   ├── project_summary.md
│   ├── methods.md
│   ├── results_summary.md
│   └── code_inventory.md
├── data/
│   ├── raw/
│   └── processed/
├── results/
│   ├── figures/
│   └── tables/
└── archive/original_colab_exports/
```

The `archive/` folder preserves the original Colab-exported scripts. The cleaned and reusable code is in `src/hlpath/` and `scripts/`.

## Data source

This project starts from the public NCBI ClinVar tab-delimited file:

```text
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
```

The data preparation script filters hearing-loss-related rows using MedGen ID `C0018784` in `PhenotypeIDS` or the phrase `hearing loss` in `PhenotypeList`. ClinVar `ClinicalSignificance` is converted into a binary label: pathogenic/likely pathogenic = 1, benign/likely benign = 0, and uncertain/conflicting labels are removed.

Generate the initial clean tables with:

```bash
python scripts/00_download_prepare_clinvar.py
```

This creates:

```text
data/raw/variant_summary.txt.gz
data/raw/hearing_loss_clinvar_clean.csv
data/processed/hearing_loss_df_clean.csv
data/processed/hearing_loss_grch38.csv
data/processed/hearing_loss_X.parquet
data/processed/hearing_loss_y.npy
```

See `docs/data_source.md` for details.

## Input data for model experiments

Expected local files:

```text
data/processed/hearing_loss_df_clean.csv
data/processed/hearing_loss_grch38_ensembl_conservation.csv
```

Large data files are intentionally not tracked. Place local CSV files under `data/raw/` or `data/processed/` before running the scripts.

Minimum required columns:

```text
label
GeneSymbol
Name or hgvs
Start
Stop
ensembl_conservation  # required for conservation and SHAP experiments
```

## Run analysis

Install:

```bash
conda env create -f environment.yml
conda activate hearingloss-variants
pip install -e .
```

Prepare ClinVar hearing-loss data:

```bash
python scripts/00_download_prepare_clinvar.py
```

Run the initial Random Forest random-split baseline:

```bash
python scripts/01_initial_random_forest_baseline.py \
  --x data/processed/hearing_loss_X.parquet \
  --y data/processed/hearing_loss_y.npy \
  --out results/tables/initial_random_forest_baseline.csv
```

Diagnose random split vs gene-based split:

```bash
python scripts/01_diagnose_random_vs_gene_split.py \
  --data data/processed/hearing_loss_df_clean.csv \
  --out results/tables/random_vs_gene_split.csv
```

Evaluate baseline vs conservation:

```bash
python scripts/02_evaluate_conservation.py \
  --data data/processed/hearing_loss_grch38_ensembl_conservation.csv \
  --out results/tables/conservation_gene_split_results.csv
```

Run heterogeneity analysis:

```bash
python scripts/03_heterogeneity_analysis.py \
  --data data/processed/hearing_loss_grch38_ensembl_conservation.csv \
  --outdir results/tables
```

Generate SHAP figures:

```bash
python scripts/04_shap_explainability.py \
  --data data/processed/hearing_loss_grch38_ensembl_conservation.csv \
  --figdir results/figures
```

## Methods summary

- **Model:** XGBoost binary classifier.
- **Features:** variant length, deletion/insertion/duplication flags, frameshift flag, stop-gain flag, splice-site flag, evolutionary conservation.
- **Validation:** gene-based train/test split; test genes are excluded from training.
- **Early stopping:** validation split is created only from training genes.
- **Metrics:** AUC and R90 threshold analysis.
- **Interpretability:** SHAP global summary, dependence plots, gray-zone SHAP, subgroup SHAP, and local error explanation.

## Limitations

- ClinVar labels can include uncertainty, submitter bias, and historical inconsistency.
- Conservation and gene constraint are informative features, not causal proof.
- Current feature engineering is coarse; future work should add transcript context, protein domains, NMD prediction, allele frequency, and tissue expression.
- Results are preliminary and intended for research/portfolio presentation.
