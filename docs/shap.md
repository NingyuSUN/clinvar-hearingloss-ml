# Feature attribution (SHAP)

`scripts/run_shap.py` explains the frozen full model (`L6_full`) with exact
TreeSHAP (XGBoost's native `pred_contribs`, identical to `shap.TreeExplainer`).
Every fold is retrained exactly as `pipeline.run_experiment` trains it — same
balanced split, same post-hoc log-loss iteration — so the explained model is the
one the reported AUCs come from. SHAP values are taken on each fold's
out-of-fold outer-test rows (never in-sample) and averaged over the 10 seeds per
variant.

```bash
pip install -e ".[interpret]"
python scripts/run_shap.py --out results/shap/
```

## Ablation and SHAP agree

Two independent methods — ablation (perturbation on held-out AUC) and SHAP
(attribution on the fitted model) — rank the missense feature groups the same
way:

| Group | Ablation ΔAUC (removed) | Ablation rank | SHAP share | SHAP rank |
|---|---:|---:|---:|---:|
| frequency | −0.136 | 1 | 49.9 % | 1 |
| gene constraint | −0.054 | 2 | 28.6 % | 2 |
| conservation | −0.040 | 3 | 18.4 % | 3 |
| protein domain | −0.004 | 4 | 2.7 % | 4 |
| consequence flags | +0.001 (ns) | 5 | 0.1 % | 6 |
| allele shape | +0.003 (ns) | 6 | 0.3 % | 5 |

Spearman rank correlation **0.943**. The bottom two swap, but both methods put
them at ≈ 0 either way.

## Headline cohort — SHAP is consistent with the ablation story there too

| Group | Share |
|---|---:|
| frequency | 36.2 % |
| consequence | 33.2 % |
| gene constraint | 17.0 % |
| conservation | 7.8 % |
| domain | 4.4 % |
| allele | 1.4 % |

By consequence class, the same model attributes very differently:

| Class | consequence | conservation | frequency | constraint |
|---|---:|---:|---:|---:|
| truncating | 57.9 % | 4.5 % | 20.4 % | 12.9 % |
| canonical splice | 49.5 % | 9.4 % | 22.7 % | 13.7 % |
| coding non-truncating (missense) | 18.1 % | 10.0 % | 45.4 % | 22.0 % |
| non-coding / other | 17.5 % | 7.8 % | 48.5 % | 15.6 % |

Even on missense rows — where the four consequence flags never vary — they still
carry 18 % of the attribution: telling the model "this is not a clear-cut
loss-of-function call" is itself informative. Frequency + gene constraint
(67 %) drive the missense-specific signal, matching the dedicated missense model.

## Conservation: a threshold, not a gradient

`figures/dependence_missense_ensembl_conservation.png` — below
`ensembl_conservation ≈ 1`, the SHAP contribution is a flat, mild pull toward
benign (≈ −0.1 to −0.3); above it, a sharp jump to a strong pull toward
pathogenic (≈ +0.1 to +0.65). It is a step, not a smooth gradient.

Binning missense variants into conservation terciles, the *relative* share
carried by each feature group barely moves (frequency 48–52 %, constraint
27–30 %, conservation 17–19 % in every band — `results/shap/conservation_bands.csv`).
So there is a real conservation threshold, but it is not accompanied by other
features "taking over" in a middle zone.

## A caveat: gene-constraint direction is not reliable

`p3_gene_constraint_oe_lof_upper` and `p3_gene_constraint_oe_mis_upper` are
gene-level constants (134 distinct values across the 3,121 missense variants)
and are themselves correlated (r = 0.62). Their raw value is positively
correlated with the label for *both* (r ≈ +0.17 to +0.20), but the SHAP-value
sign for `oe_lof_upper` points the other way — the classic effect of fitting two
collinear features together. **Report the group's importance (rank 2 for
missense), not a directional claim for either constraint feature individually**;
disentangling them would need a model with only one of the two, or a partial
dependence analysis. Frequency and conservation directions are clean and match
their raw label correlation.

## Reproduce

```bash
python scripts/run_shap.py --out results/shap/
```

Writes `importance_{headline,missense}.csv`, `group_importance_{headline,missense}.csv`,
`by_consequence.csv`, `conservation_bands.csv`, and `figures/`.
