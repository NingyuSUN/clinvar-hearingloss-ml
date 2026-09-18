# Changelog

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
