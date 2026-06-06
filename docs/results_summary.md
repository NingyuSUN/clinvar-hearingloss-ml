# Results Summary

## Random split vs gene-based split

| Evaluation setting | AUC |
|---|---:|
| Random split | 0.878 |
| Gene-based split | 0.723 |

Random split produced higher apparent performance than gene-based split, suggesting that standard random evaluation can overestimate generalization.

## Conservation experiment

| Model | AUC | R90 FN | R90 FP | R90 Recall | R90 Precision |
|---|---:|---:|---:|---:|---:|
| Baseline | 0.653 | 1 | 664 | 0.998 | 0.449 |
| + Conservation | 0.812 | 54 | 394 | 0.900 | 0.553 |

Adding conservation improved AUC by 0.159 and reduced R90 false positives by 270.

## SHAP interpretation

Main qualitative findings:

- `ensembl_conservation` is a major global decision feature.
- High conservation tends to push predictions toward pathogenic.
- Low conservation tends to push predictions toward benign.
- In the approximate gray zone from -1 to 1, conservation becomes less decisive and structural features gain importance.
- Stop-gain effects can be context-dependent, suggesting future analysis with transcript position, NMD prediction, and gene-level LoF tolerance.

## Ablation Study Interpretation
| Model	| Cohort N | Genes	| AUC	| R90 Precision|	R90 Recall |
|---|---:|---:|---:|---:|---:|
| Base	| 7827	| 198	| 0.828 ± 0.053	| 0.650 ± 0.146	| 0.904 ± 0.005 |
| Base + LoF | 7827	| 198	| 0.786 ± 0.071	| 0.627 ± 0.122	| 0.901 ± 0.000 |
| Base + LoF + Missense	| 7827	| 198	| 0.786 ± 0.088	| 0.642 ± 0.124	| 0.901 ± 0.002 |
| Base + LoF + Missense + MAF	| 7827	| 198	| 0.782 ± 0.087	| 0.641 ± 0.125	| 0.902 ± 0.001 |
| Full Optimized	| 7827	| 198	| 0.787 ± 0.084	| 0.640 ± 0.130	| 0.901 ± 0.000 |

We utilized 3 kinds of ablation study: Previous feasure-specific dropna, Feature-complete subset, and Fixed full cohort+imputation. 
The study indicates that more features don't necessarily become a better model. 

Based on the 3 ablation study, we concluded that:
- The full optimized model has small increases in the feature-complete subset
- The full optimized model is lower than the base model when broader fixed full cohort+imputation is used
- The missingness of the gene constraint annotations is a significant limiting factor in deciding the efficacy of the model

## Current conclusion

Evolutionary conservation provides useful cross-gene signal for hearing-loss variant pathogenicity prediction under strict unseen-gene evaluation.
