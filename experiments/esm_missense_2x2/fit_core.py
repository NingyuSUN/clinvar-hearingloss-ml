'''Single-fold XGBoost fit: train-only imputation, post hoc iteration, validation-only R90.'''
from __future__ import annotations

import numpy as np
import pandas as pd
import xgboost as xgb

from hlpath.config import NUM_BOOST_ROUND, POSTHOC_MIN_ITER, xgb_params
from hlpath.metrics import r90_threshold


def _as_numeric(df: pd.DataFrame, name: str) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        try:
            out[c] = pd.to_numeric(out[c], errors='raise').astype('float64')
        except (TypeError, ValueError) as exc:
            raise ValueError(f'{name}: column {c!r} is not numeric') from exc
    if np.isinf(out.to_numpy(dtype='float64')).any():
        raise ValueError(f'{name}: infinite values present')
    return out


def impute_train_only(train: pd.DataFrame, *others: pd.DataFrame, continuous):
    '''Continuous columns get train medians; all other columns get 0.'''
    cols = list(train.columns)
    if len(set(cols)) != len(cols):
        raise ValueError('duplicate feature columns in train')
    cont = set(continuous)
    frames = [_as_numeric(train, 'train')]
    for i, other in enumerate(others):
        if set(other.columns) != set(cols) or len(other.columns) != len(cols):
            raise ValueError(f'frame {i}: columns differ from train')
        frames.append(_as_numeric(other.loc[:, cols], f'frame {i}'))

    medians: dict = {}
    for c in cols:
        if c in cont:
            med = frames[0][c].median(skipna=True)
            if pd.isna(med) or not np.isfinite(med):
                raise ValueError(f'continuous column {c!r} is entirely missing in train')
            medians[c] = float(med)

    fixed = []
    for i, df in enumerate(frames):
        for c in cols:
            df[c] = df[c].fillna(medians.get(c, 0.0))
        if df.isna().to_numpy().any():
            raise ValueError(f'frame {i}: NaN remains after imputation')
        fixed.append(df)
    return tuple(fixed), medians


def select_iteration(curve) -> int:
    '''Zero-indexed argmin of validation logloss, searching from POSTHOC_MIN_ITER onward.'''
    ll = np.asarray(curve, dtype='float64')
    if ll.ndim != 1 or len(ll) <= POSTHOC_MIN_ITER:
        raise ValueError(f'logloss curve needs more than {POSTHOC_MIN_ITER} entries')
    if not np.isfinite(ll).all():
        raise ValueError('non-finite validation logloss')
    return int(np.argmin(ll[POSTHOC_MIN_ITER:])) + POSTHOC_MIN_ITER


def _labels(y, n_rows: int, name: str) -> np.ndarray:
    arr = np.asarray(y)
    if arr.ndim != 1 or len(arr) != n_rows:
        raise ValueError(f'{name}: expected 1-D labels with {n_rows} rows, got shape {arr.shape}')
    try:
        vals = arr.astype('float64')
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name}: labels are not numeric') from exc
    if not np.isfinite(vals).all() or not np.isin(vals, (0.0, 1.0)).all():
        raise ValueError(f'{name}: labels must be binary 0/1')
    out = vals.astype(int)
    if len(np.unique(out)) < 2:
        raise ValueError(f'{name}: both classes required')
    return out


def _check_features(features, frames) -> list:
    feats = list(features)
    if not feats:
        raise ValueError('no features given')
    if len(set(feats)) != len(feats):
        raise ValueError('duplicate feature names')
    if not all(isinstance(f, str) for f in feats):
        raise ValueError('feature names must be strings')
    for name, df in frames:
        missing = [f for f in feats if f not in df.columns]
        if missing:
            raise ValueError(f'{name}: missing features {missing}')
    return feats


def fit_one(X_train, y_train, X_validation, y_validation, X_test, features, continuous,
            seed, threads) -> dict:
    feats = _check_features(features, [('X_train', X_train), ('X_validation', X_validation),
                                       ('X_test', X_test)])
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise ValueError('seed must be an integer')
    if isinstance(threads, bool) or not isinstance(threads, (int, np.integer)) or threads <= 0:
        raise ValueError('threads must be a positive integer')

    (xtr, xva, xte), medians = impute_train_only(
        X_train.loc[:, feats], X_validation.loc[:, feats], X_test.loc[:, feats],
        continuous=continuous)
    ytr = _labels(y_train, len(xtr), 'y_train')
    yva = _labels(y_validation, len(xva), 'y_validation')

    params = dict(xgb_params(int(seed)))
    params['nthread'] = int(threads)

    dtr = xgb.DMatrix(xtr.to_numpy(dtype='float64'), label=ytr, feature_names=feats)
    dva = xgb.DMatrix(xva.to_numpy(dtype='float64'), label=yva, feature_names=feats)
    dte = xgb.DMatrix(xte.to_numpy(dtype='float64'), feature_names=feats)

    hist: dict = {}
    booster = xgb.train(params, dtr, num_boost_round=NUM_BOOST_ROUND,
                        evals=[(dva, 'validation')], evals_result=hist, verbose_eval=False)
    curve = [float(v) for v in hist['validation']['logloss']]
    if len(curve) != NUM_BOOST_ROUND:
        raise RuntimeError(f'expected {NUM_BOOST_ROUND} logloss entries, got {len(curve)}')
    best_iter = select_iteration(curve)

    rng = (0, best_iter + 1)
    val_scores = np.asarray(booster.predict(dva, iteration_range=rng), dtype=np.float64)
    test_scores = np.asarray(booster.predict(dte, iteration_range=rng), dtype=np.float64)
    if not np.isfinite(test_scores).all():
        raise RuntimeError('non-finite test scores')

    thr, ok, reason = r90_threshold(yva, val_scores)
    if not ok or thr is None or not np.isfinite(thr):
        raise RuntimeError(f'R90 threshold unavailable: {reason}')

    return {
        'model': booster,
        'validation_scores': val_scores,
        'test_scores': test_scores,
        'imputation_medians': medians,
        'validation_logloss': curve,
        'best_iteration': int(best_iter),
        'threshold': float(thr),
        'threshold_reason': str(reason),
        'xgb_params': params,
    }
