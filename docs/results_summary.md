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

## Current conclusion

Evolutionary conservation provides useful cross-gene signal for hearing-loss variant pathogenicity prediction under strict unseen-gene evaluation.
