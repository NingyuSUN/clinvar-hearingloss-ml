# Data — feature construction

The modelling matrix (`data/modeling_matrix.csv.gz`) is one row per ClinVar
variant with the label, the gene group, the review star, the cohort flags, the
consequence class, and 18 model features. This document describes how it is built;
`scripts/build_matrix.py` runs `hlpath.protocol.build_matrix` on a raw annotated
table.

## 1. Variant selection

Start from NCBI ClinVar `variant_summary.txt.gz`
(`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/`). Keep GRCh38 rows
whose `PhenotypeIDS` contains MedGen `C0018784` or whose `PhenotypeList` contains
"hearing loss". Upstream, VUS and Conflicting rows are removed by the label rule
(`docs/methods.md`).

## 2. Consequence terms — Ensembl VEP

Each variant is annotated through the Ensembl VEP REST API (GRCh38). The
transcript is chosen with a target-gene → MANE Select → canonical →
protein-coding → any preference, and the consequence terms from that transcript
give the four flags:

| Flag | From VEP terms |
|---|---|
| `p1_has_frameshift` | `frameshift_variant` |
| `p1_has_stop_gained` | `stop_gained` |
| `p1_has_canonical_splice` | `splice_donor_variant`, `splice_acceptor_variant` |
| `p1_has_splice_region` | `splice_region_variant` |

`p1_status` records whether a single reliable transcript was found
(`ok_single` / `ok_consensus`) or the variant key / transcript was unreliable.
A `c.*` 3′-UTR position is not a stop gain; positions off the reference transcript
do not contribute flags.

## 3. Allele frequency — gnomAD v4.1.1

Each variant is queried against the gnomAD v4.1.1 GraphQL API.

- `p2_af_status` = `LIVE_PASS` (observed, PASS filter), `NOT_IN_GNOMAD_R4`
  (confirmed absent), `LIVE_FILTERED`, or an unresolved state.
- `p4_freq_log10` = `log10(max PASS allele frequency + 1e-8)` for `LIVE_PASS`;
  a floor `log10(1 / 2N_max + 1e-8)` (≈ 6.25 × 10⁻⁷) for `NOT_IN_GNOMAD_R4`;
  `NaN` otherwise (imputed per fold).
- `p2_af_known_missing`, `p2_not_observed_in_gnomad_r4` are indicator flags.
- `p2_low_an_flag` marks an allele number below 2,000 (below the ACMG `BA1` /
  `BS1` minimum).

## 4. Protein domains — Ensembl VEP `domains`

Overlap of the variant's protein position with a curated domain from Pfam, SMART,
PROSITE (profiles / patterns), PRINTS, or CDD. Structure-model hits
(PANTHER, Gene3D, SUPERFAMILY, …) and non-functional tracks (mappings, mobidb,
phobius) are excluded.

- `p3_in_specific_domain` — inside a curated domain (0 for non-coding, not `NaN`).
- `p3_specific_domain_count` — number of overlapping curated domains.
- `p3_variant_is_coding` — 1 for coding, 0 for non-coding.

## 5. Gene constraint — gnomAD v2.1.1

`p3_gene_constraint_oe_lof_upper` and `p3_gene_constraint_oe_mis_upper` are the
upper bounds of the observed/expected ratio for LoF and missense variation, from
`gnomad.v2.1.1.lof_metrics.by_gene`, matched by current HGNC stable ID. gnomAD v2
and v4 constraint are never mixed. Multi-gene rows are not assigned an arbitrary
first value.

## 6. Conservation

`ensembl_conservation` is a single comparative-genomics conservation score at the
variant position. Its provenance (release, coordinate mapping) is **not fully
re-audited** in this cycle — see `docs/limitations.md`.

## Matrix schema

| Column | Meaning |
|---|---|
| `VariationID` | ClinVar variation id |
| `y` | label (1 = P/LP, 0 = B/LB) |
| `p4_gene_group` | connected-component gene group id |
| `p4_review_star` | 0–4 |
| `p4_headline_cohort` / `p4_in_modeling_cohort` / `p4_sens_ge2star` / `p4_expert_panel` | cohort membership flags |
| `consequence_class` | `truncating` / `canonical_splice` / `coding_nontruncating` / `noncoding_or_other` |
| `GeneSymbol`, `ReviewStatus` | kept for reporting, never used as features |
| 18 feature columns | see `docs/methods.md` |
