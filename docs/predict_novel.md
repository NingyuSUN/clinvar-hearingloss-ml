# Experimental live annotation

`predict/predict_novel.py` queries live sources for GRCh38 single-nucleotide
variants. The reproducible review entrypoint remains the
[frozen-cache demo](predict.md). Live annotation has no independent external
validation and produces uncalibrated research scores.

## Install and run

```bash
pip install -r predict/requirements-inference.txt -r predict/requirements-novel.txt
python predict/predict_novel.py --input my_variants.csv --bundle predict/bundle --output my_output/
python predict/render_variant_report.py --predictions my_output/predictions.jsonl
```

Use a fresh output directory. CSV columns:
`variant_id,assembly,chrom,pos,ref,alt`. Assembly must be `GRCh38`.
Optionally install `alphagenome` and set `ALPHAGENOME_API_KEY` to request
AVI/splice annotations. An API key does not resolve incompatible model features.

**Only GPN3 can currently score live inputs**, provided identity, supported
consequence class, gene mapping and all three GPN scores pass validation.
The four models requiring `ensembl_conservation` return
`missing_or_invalid_features`: the original conservation source is unresolved.
UCSC phyloP is retained as a separate reference annotation and never passed to
these trained models. No substitute or median bypass is provided.

## Data flow and contracts

| Source | Role | Failure behavior |
|---|---|---|
| Ensembl VEP | Consequence, gene, domains; MANE > canonical > coding > any transcript | Check response assembly, chromosome, position, strand and REF/ALT; reject mismatches and malformed transcripts. |
| gnomAD v4 frequency | Max genome/exome PASS AF | GraphQL errors and invalid payloads are failures, not absence. Require explicit filters, finite AF in [0,1] and positive AN. |
| gnomAD v2.1.1 constraint | LoF/missense observed/expected upper bounds | Missing values remain missing; reject malformed or nonfinite values. |
| UCSC phyloP100way | Reference annotation only | Original model feature remains unavailable even if phyloP succeeds. |
| Public GPN-Star Parquet | Three precomputed scores | Missing dependencies, failed lookups or ambiguous/nonfinite matches remain missing. No local foundation-model inference. |
| Optional AlphaGenome | AVI and splice annotations | Missing key or failed query remains missing. AVI is a composite, not an AlphaGenome-only model. |

Only an error-free explicit `data.variant: null` is treated as absent from
gnomAD. Its transform is `log10(1 / 1,600,000 + 1e-8)`, matching the frozen
protocol and 2,085 absent cache records. Observed AF zero instead maps to -8;
filtered/unavailable AF cannot use the cache-only median-imputation allowance.

Every input keeps its row index, including duplicates and failures. Known
cache membership comes from the verified bundle; membership outside that cache
is unknown (`null`), not proof of a new sample. Live consequence/gene context
conflicting with a cached record is rejected. Selected models must explicitly
exclude the known held-out group; a broken bundle fails the run.

## Boundaries

The original VEP target-gene preference tier is not reconstructed. A historical
spot check on `16-2496587-G-C` is not a cohort-level compatibility study. Source
versions and annotations may change. HTTP sources use request timeouts; GPN's
remote Parquet transport is managed by Polars and has no CLI-wide deadline.
Successful live annotation does not establish probability calibration,
independent external validity or clinical utility. The frozen demo is unchanged.

## Tests

Default offline tests exercise malformed responses, transport failures,
identity checks, source compatibility, frozen frequency parity, held-out-group
integrity and row-preserving CLI output:

```bash
pip install pytest requests
python -m pytest tests/test_live_sources_unit.py -q
```

Optional service checks require network access and can drift or fail with the
services. They are not part of required CI and were not rerun for this revision:

```bash
RUN_NETWORK_TESTS=1 python -m pytest tests/test_predict_novel.py -v
```

Response contracts: [Ensembl VEP](https://rest.ensembl.org/documentation/info/vep_region_get)
and [GraphQL response specification](https://spec.graphql.org/September2025/#sec-Response).
