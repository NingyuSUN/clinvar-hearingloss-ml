# data/

`modeling_matrix.csv.gz` — one row per ClinVar hearing-loss variant:

- `y` — label (1 = P/LP, 0 = B/LB)
- `p4_gene_group` — connected-component gene group id (the CV split unit)
- `p4_review_star` — ClinVar gold stars, 0–4
- `p4_headline_cohort`, `p4_in_modeling_cohort`, `p4_sens_ge2star`,
  `p4_expert_panel` — cohort membership flags
- `consequence_class` — `truncating` / `canonical_splice` / `coding_nontruncating`
  / `noncoding_or_other`
- `GeneSymbol`, `ReviewStatus` — reporting only, never model features
- 18 feature columns — see `docs/methods.md`

Built from ClinVar + Ensembl VEP + gnomAD annotations; see `docs/data.md`. To
rebuild it from a raw annotated table: `python scripts/build_matrix.py --raw ...`.

`hlpath.features.load_matrix()` reads the `.csv` or the `.csv.gz`; set `HLPATH_DATA`
to point elsewhere.
