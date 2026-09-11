# results/

Frozen-run outputs cited by `README.md` and `docs/results.md`. Re-running
`scripts/run_evaluation.py` then `scripts/analyze_results.py` overwrites these
with near-identical numbers (the bootstrap CI moves at the third decimal).

| File | Contents |
|---|---|
| `ladder_headline.csv` / `ladder_missense.csv` | additive feature ladder, per-fold AUC |
| `paired_ladder_headline.csv` / `paired_ladder_missense.csv` | paired per-fold Δ and t for each ladder step |
| `bootstrap_ci.csv` | patient-level bootstrap 95% CI, per cohort and feature set |
| `stratified_headline.csv` | headline-model recall / precision by consequence class |
| `split_gap.csv` | gene-held-out vs random-split pooled AUC |
| `pooled_oof_by_cohort.csv` | pooled out-of-fold confusion, per cohort |
| `run_manifest.json` | XGBoost version, seeds, matrix hash |
