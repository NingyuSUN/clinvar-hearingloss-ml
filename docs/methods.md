# Methods

## Data

The project uses a cleaned ClinVar-derived table of hearing-loss-related variants. The main modeling table contains variant identifiers, gene symbols, clinical labels, genomic positions, HGVS-like variant descriptions, and external conservation annotation.

Expected label encoding:

- `0`: benign / likely benign
- `1`: pathogenic / likely pathogenic

## Feature engineering

The cleaned pipeline uses variant-centric features only. `GeneSymbol` is used for splitting, not as a model feature.

Core features:

- `variant_length`
- `has_del`
- `has_ins`
- `has_dup`
- `has_fs`
- `has_stop`
- `has_splice`
- `ensembl_conservation`

Optional future features:

- gnomAD gene constraint
- allele frequency
- transcript-level annotation
- protein domain annotation
- NMD prediction
- tissue expression

## Validation design

The main evaluation uses gene-based splitting:

- Training genes and test genes are disjoint.
- Early stopping uses only a validation subset from training genes.
- Test genes are never used for model selection.

This design tests cross-gene generalization rather than random-split memorization.

## Model

The main model is an XGBoost binary classifier with early stopping. The repository keeps the training logic in `src/hlpath/modeling.py` and runnable entry points in `scripts/`.

## Metrics

The project reports:

- AUC for threshold-independent ranking performance.
- Confusion matrix metrics at fixed threshold 0.5.
- R90 metrics, where the threshold is selected to satisfy recall >= 0.90 and then FN, FP, precision, and recall are reported.

## Interpretability

SHAP analysis is used to examine:

- global feature importance
- conservation dependence
- gray-zone behavior
- splice vs non-splice subgroups
- representative false positives and false negatives
