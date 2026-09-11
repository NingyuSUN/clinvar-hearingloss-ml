# Notebooks

There are no notebooks in this repository. The analysis is a small set of runnable
scripts (`scripts/`) over a reusable package (`src/hlpath/`):

```bash
python scripts/run_evaluation.py --out results/
python scripts/analyze_results.py --dir results/
```

The frozen protocol is `src/hlpath/protocol.py`; the cross-validation is
`src/hlpath/pipeline.py`. See `docs/methods.md`.
