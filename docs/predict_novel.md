# Scoring a genuinely novel variant

`docs/predict.md` covers `predict_variants.py`, which only replays the ~4,050
variants already in the frozen, hash-verified cache. This page covers
`predict/predict_novel.py`, which scores a variant that **isn't** in that
cache, by calling live public data sources (plus, optionally, the AlphaGenome
API with your own key) at request time.

**Read this whole page before trusting the output.** This path has much
weaker guarantees than the cache-only path, and it says so in its own
`run_manifest.json` every time it runs.

## Install

```bash
pip install -r predict/requirements-inference.txt
pip install -r predict/requirements-novel.txt   # requests, polars, and optionally the alphagenome SDK
```

`alphagenome` is listed but not required to run the script — omit an API key
and the script will still score three of the five models (see below).

## Run

```bash
# Without an AlphaGenome key: ORIGINAL18, GPN3, ORIGINAL18_GPN are scored;
# FULL_UNIFIED and FULL_STRATIFIED report missing_or_invalid_features.
python predict/predict_novel.py --input my_variants.csv --bundle predict/bundle --output my_output/

# With your own AlphaGenome key, apply at https://deepmind.google.com/science/alphagenome:
export ALPHAGENOME_API_KEY=...
python predict/predict_novel.py --input my_variants.csv --bundle predict/bundle --output my_output/

python predict/render_variant_report.py --predictions my_output/predictions.jsonl
```

Input CSV columns are the same as `predict_variants.py`:
`variant_id,assembly,chrom,pos,ref,alt`, GRCh38, single-nucleotide
substitutions only.

## What actually happens per variant

| Step | Source | Network? | Notes |
|---|---|---|---|
| Consequence class, gene, protein-domain overlap | Ensembl VEP REST | Yes (public, no key) | Transcript picked by MANE Select > canonical > protein-coding > any. |
| Allele frequency | gnomAD v4.1.1 GraphQL | Yes (public, no key) | Same max(genome, exome) PASS rule as `docs/data.md`. |
| Gene constraint (`oe_lof_upper`/`oe_mis_upper`) | gnomAD v2.1.1, via gnomAD's own API | Yes (public, no key) | |
| Conservation | UCSC phyloP (100-way), via UCSC's public REST API | Yes (public, no key) | **A deliberate substitute** — see below. |
| GPN scores (×3 checkpoints) | `songlab/gpn-star-scores` on Hugging Face, remote Parquet row lookup | Yes (public, no key) | No local model, no GPU — this is a precomputed public dataset, not live inference. |
| AVI / splice_sites | AlphaGenome API | Yes (needs **your own key**) | Skipped if no key given; the two models that need it report missing, not a guess. |

Every one of these fetch functions (except AlphaGenome, which needs a key I
don't have) was checked against a variant already in the frozen cache
(`16-2496587-G-C`, TBC1D24) before this was shipped: VEP's transcript choice,
the domain-overlap count, the allele frequency, the gene constraint values,
and all three GPN scores reproduced the cache's values **exactly**. The
resulting `ORIGINAL18`/`ORIGINAL18_GPN` scores differ very slightly from the
cache's for that same variant — that's the conservation substitution below
propagating through, not a bug.

## The one substitution you should know about: conservation

The frozen cache's `ensembl_conservation` feature is described in
`docs/data.md` as "a single comparative-genomics conservation score," and
`docs/limitations.md` already says its exact provenance was never fully
audited. Investigating this further while building this page: Ensembl's own
GERP conservation score (computed from EPO_EXTENDED/Pecan alignments) has
**no REST API** — the only ways to get it are the Ensembl Compara Perl API or
decoding packed binary data directly out of Ensembl's `conservation_score`
MySQL table, neither practical to ask a downloader to set up.

So `predict_novel.py` uses **UCSC phyloP (100-way vertebrate)** instead, via
UCSC's public `api.genome.ucsc.edu` REST endpoint. This is a related but
**different** conservation measure, on a different numeric scale, from a
different multiple-sequence alignment. It is not a reproduction of the
original feature — it's a documented, deliberate stand-in, chosen so the
other 17 hand-engineered features (which are exactly reproducible) don't all
have to sit unused behind one unfixable feature. If you have access to the
original conservation source, `docs/predict_novel.md`'s "Open" item in
`PROJECT_STATUS.md` tracks swapping it in properly.

## What's simplified versus the frozen pipeline

- **VEP transcript selection** drops the "target-gene" preference tier
  (`docs/data.md`'s first tier, before MANE Select) — the internal
  hearing-loss gene-panel list that tier depended on wasn't preserved
  anywhere this repository could find. MANE Select > canonical >
  protein-coding > any is used instead. This matched the frozen cache
  exactly on every case checked, but isn't a guarantee for every gene.
- **Only the four supported consequence classes score at all**
  (`missense`, `canonical_splice`, `splice_region`, `synonymous`) — same
  restriction as the cache-only path. A truncating or non-coding variant
  correctly comes back `unsupported_variant_type`.
- **No probability, ever.** Same as `predict_variants.py`: scores are
  uncalibrated classifier outputs, and this live path has had *zero*
  external validation — it is more experimental than the already-uncalibrated
  cache-only path, not less.

## Failure modes you'll actually see

- **Wrong REF allele**: Ensembl VEP will say so directly (`"request for
  consequence of [X] matches reference [X]"`), surfaced verbatim in
  `annotation_status` as `vep_failed:vep_http_400:...`. This is usually a
  typo in your input, not a pipeline bug.
- **Gene not resolvable to a known gene-group**: the variant still gets
  scored, but with `novel_gene_context` instead of `known_gene_held_out` —
  there is no held-out-fold guarantee for a gene the models never trained on,
  and the report says so.
- **No AlphaGenome key**: `FULL_UNIFIED`/`FULL_STRATIFIED` report
  `missing_or_invalid_features`, not a fabricated score.

## Testing this yourself

`tests/test_predict_novel.py` exercises these live sources for real (no
mocks) but is skipped unless you opt in, since it depends on external
services that can fail or change independently of this repository:

```bash
RUN_NETWORK_TESTS=1 python -m pytest tests/test_predict_novel.py -v
```
