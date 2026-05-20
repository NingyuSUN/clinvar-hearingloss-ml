# Code Inventory

## Cleaned code

| Path | Purpose |
|---|---|
| `src/hlpath/data.py` | Load data and create variant-level features. |
| `src/hlpath/metrics.py` | Compute confusion metrics and R90 threshold metrics. |
| `src/hlpath/modeling.py` | Train and evaluate XGBoost under gene-based splitting. |
| `src/hlpath/heterogeneity.py` | Run subgroup and conservation-bin analyses. |
| `src/hlpath/shap_analysis.py` | Generate SHAP values and plots. |
| `scripts/02_evaluate_conservation.py` | Main baseline vs conservation experiment. |
| `scripts/03_heterogeneity_analysis.py` | Subgroup and gray-zone style result tables. |
| `scripts/04_shap_explainability.py` | SHAP global and gray-zone figure generation. |

## Original Colab exports

The original scripts are preserved under `archive/original_colab_exports/` for traceability. They include earlier exploratory, teaching-oriented, and repeated versions of the analysis.

| Original file | Role in project |
|---|---|
| `hl_step3 (1).py` | Early random-split model comparison with RF/XGBoost/LightGBM. |
| `hl_4.py` | Diagnostic comparison between random split and gene-based split. |
| `hl_5 (2).py` | GRCh38 filtering and conservation annotation workflow. |
| `hl_+cons.py` | Baseline vs conservation under strict gene-based evaluation. |
| `hearing_loss_heterogeneity.py` | Heterogeneity, error, and conservation-boundary analysis. |
| `hearingloss_shap*.py` | SHAP explainability; multiple repeated exports were retained. |
