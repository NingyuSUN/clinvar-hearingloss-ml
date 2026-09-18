# Changelog

## 2026-09-18 — Add live annotation for genuinely novel variants (`predict_novel.py`)

- Added `predict/live_sources.py` and `predict/predict_novel.py`: a second
  entry point that annotates a variant **not** in `predict/bundle`'s frozen
  cache by calling live sources at request time — Ensembl VEP REST (consequence,
  gene, protein-domain overlap), gnomAD v4.1.1 GraphQL (allele frequency),
  gnomAD v2.1.1 via gnomAD's own API (gene constraint), a remote Hugging
  Face-hosted Parquet lookup against `songlab/gpn-star-scores` for all three
  GPN checkpoints (no local model or GPU — GPN-Star is a precomputed public
  dataset), and optionally the AlphaGenome API (the user's own key) for
  `avi`/`splice_sites`. Scoring reuses `prediction_core.py`'s existing
  building blocks (`extract_features`, `select_fold`, `prediction_set`,
  `Bundle.model`) without touching or weakening its hash-verified,
  cache-only `predict_batch` path, which stays exactly as it was.
- **Verified against the frozen cache before shipping**: every one of VEP's
  transcript selection, the domain-curation count, the allele frequency, the
  gene-constraint values, and all three GPN scores reproduced a known cached
  variant's values exactly (tested against two variants of different
  consequence classes and genes). The one deliberate exception is
  `ensembl_conservation`, which uses UCSC phyloP (100-way) via UCSC's public
  REST API as a documented substitute — see "Known gaps" below and
  `docs/predict_novel.md`.
- Added `docs/predict_novel.md` (the tutorial), `predict/requirements-novel.txt`,
  and `tests/test_predict_novel.py` (real network calls against live services,
  skipped by default; opt in with `RUN_NETWORK_TESTS=1`, not part of required
  CI since external-service flakiness shouldn't block this repository's own
  checks).
- Updated `README.md`, `MODEL_CARD.md` and `PROJECT_STATUS.md`, including
  correcting `MODEL_CARD.md`'s previous "not a novel-variant annotation
  service" line, which this change makes no longer fully accurate — replaced
  with an explanation of the two entry points' very different guarantees.
- **Known gaps, stated rather than hidden**: this path has had no external
  validation of its own; `ensembl_conservation`'s live substitute is on a
  different scale from the original (unaudited) training-data source, so
  `ORIGINAL18`/`ORIGINAL18_GPN`/`FULL_UNIFIED`/`FULL_STRATIFIED` scores for a
  novel variant will differ somewhat from what an identical cached variant
  would have shown; the frozen pipeline's VEP "target-gene" transcript
  preference tier is not reproduced (the underlying gene-panel list wasn't
  preserved anywhere this repository could find); and this path needs live
  network access, unlike `predict_variants.py`.

## 2026-09-18 — Add downloadable five-model prediction CLI (`predict/`)

- Added `predict/`: a CLI (`predict_variants.py` + `prediction_core.py`) that
  scores a variant against five related models (`ORIGINAL18`, `GPN3`,
  `ORIGINAL18_GPN`, `FULL_UNIFIED`, `FULL_STRATIFIED`) and a new
  `render_variant_report.py` that turns the structured output into a
  plain-language Markdown report — the "detailed analysis" piece that an
  earlier engineering pass on this same prediction core had explicitly
  deferred. See `docs/predict.md` (new) and the `MODEL_CARD.md` update.
- The bundle (`predict/bundle/`: hash-verified model weights + a frozen,
  pre-verified annotation cache for ~4,050 variants) and the two inference
  scripts originate from a prior engineering pass (dated 2026-09-15) that had
  been run and verified on a separate compute server but never pushed to this
  repository, pending further validation work. Before folding it in here:
  - `prediction_core.py` and `predict_variants.py` were reformatted from a
    deliberately compact single-file style into normal multi-line Python for
    readability. The reformatted code was checked to produce **byte-identical**
    `predictions.jsonl` / `predictions.csv` / `predictions_wide.csv` against
    the original code, on the shipped examples and on the full 7,125-variant
    development cohort (35,625 model rows), before being committed.
  - Every file under `predict/bundle/` was copied byte-for-byte and re-verified
    against its recorded sha256 in `bundle_manifest.json` after the copy.
- Added `tests/test_predict.py` (needs `predict/requirements-inference.txt`;
  skipped otherwise) and `tests/test_render_variant_report.py`
  (dependency-free), including a regression test pinning known scores for one
  variant so a future change can't silently drift the numbers.
- Added a `predict-test` / `predict-demo` target to `Makefile`, and a CI step
  that runs the full CLI + report generator on every push.
- **Known gaps carried over, not fixed here** (tracked in `PROJECT_STATUS.md`):
  this only scores variants already in the frozen cache — a genuinely novel
  variant returns null on all five models by design, since an earlier version
  of this tool had a real bug where caller-supplied external annotations could
  be accepted without being verified. Scores are uncalibrated, not
  probabilities. There is no dedicated AlphaGenome-only model — the closest
  things (`FULL_UNIFIED`/`FULL_STRATIFIED`) use a composite AVI score that
  also includes AlphaMissense and related signals, not pure AlphaGenome.
- The prior engineering pass's own working notes said not to push this to
  GitHub yet, pending that further validation work. Publishing it now, with
  the scope limitations above stated plainly in `docs/predict.md` and
  `MODEL_CARD.md`, was a deliberate decision to revisit that hold rather than
  an oversight.

## 2026-09-18 — Repository hardening

- Added `LICENSE` (MIT), `CITATION.cff`, `MODEL_CARD.md`, `PROJECT_STATUS.md`.
- Added `.gitattributes` (`text=auto eol=lf`) so a text file re-saved with
  CRLF line endings by an external tool can no longer show up as a false
  "modified" diff — this closed a gap surfaced the same day, when a batch of
  `results/*.csv`/`*.json` files was found to differ from HEAD only in line
  endings, not content.
- Added a dependency-free artifact validator (`tools/validate_public_artifacts.py`)
  and a pytest suite (`tests/`) that check: no CRLF regressions, `results/summary.json`
  and `results/run_manifest.json` have their expected keys and a well-formed
  dataset hash, every `results/*.csv` is non-empty with a valid header, and no
  `__pycache__` directory is tracked by git.
- Added GitHub Actions CI (`.github/workflows/ci.yml`) running the test suite,
  the artifact validator, and a Python syntax check on every push/PR — none of
  it retrains a model, so CI stays fast and needs no private data.
- Added a `Makefile` (`test`, `validate`, `syntax` targets).

## 2026-09-18 — Revert stratification-correction documentation change

- Reverted commit `1235591` ("Add stratification correction: separate
  missense/synonymous/other in the coding non-truncating cohort"), which had
  been sitting as an un-reviewed local commit — authored by a collaborator,
  not yet reviewed by the repository owner — and had already been pushed to
  `main`. No model or result file was touched by either the original commit
  or the revert; only `README.md` and `docs/*.md` wording changed. Pending
  review, `README.md` again describes the 3,121-row coding non-truncating
  cohort as "missense"; see `MODEL_CARD.md` and `PROJECT_STATUS.md` for the
  accuracy caveat this leaves open.

## 2026-09-11 — Add stratification correction (reverted 2026-09-18, see above)

## 2026-09-11 — Add TreeSHAP feature attribution

- Added `scripts/run_shap.py` / `src/hlpath/shap_analysis.py` and
  `docs/shap.md`: exact TreeSHAP attribution on the frozen fitted models,
  cross-checked against the feature ablation.

## 2026-09-11 — Rewrite: frozen protocol, gene-held-out evaluation, honest stratified results

- Full rewrite of the evaluation around a frozen protocol (label rule, cohort
  gate, gene-group split, feature list, operating point fixed before
  training), replacing the earlier ad hoc results.

## 2026-06-06 — Update results_summary.md

## 2026-06-04 — Update README.md

## 2026-05-23 — Update project_summary.md

## 2026-05-20 — Initial upload

- First commit of data, scripts and README.
