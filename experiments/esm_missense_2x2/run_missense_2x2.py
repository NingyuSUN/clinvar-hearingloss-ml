#!/usr/bin/env python3
"""Frozen four-arm runner with matched missense validation and test records."""
from __future__ import annotations
import argparse
import contextlib
import fcntl
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback
import uuid
from datetime import datetime, timezone


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp.' + uuid.uuid4().hex)
    tmp.write_text(json.dumps(data, ensure_ascii=False, sort_keys=True,
                              indent=2, allow_nan=False) + '\n')
    os.replace(tmp, path)


def atomic_csv(path, frame):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp.' + uuid.uuid4().hex)
    frame.to_csv(tmp, index=False, lineterminator='\n', float_format='%.17g')
    os.replace(tmp, path)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def libs(repo_root=None):
    if repo_root:
        sys.path.insert(0, str(Path(repo_root).resolve() / 'src'))
    global np, pd, config, protocol, pipeline, plm, features, metrics, fit_core
    import numpy as np
    import pandas as pd
    from hlpath import config, protocol, pipeline, plm, features, metrics
    import fit_core


def counts(frame):
    return {'n': len(frame), 'positive': int((frame.y == 1).sum()),
            'negative': int((frame.y == 0).sum()),
            'gene_groups': int(frame.gene_group.nunique())}


def both(frame, name):
    need(len(frame) > 0 and set(frame.y.unique()) == {0, 1},
         f'{name}: require nonempty both classes')


def validate_ids(frame, name):
    need('VariationID' in frame, f'{name}: missing VariationID')
    need(not frame.VariationID.isna().any(), f'{name}: missing ID')
    need(not frame.VariationID.astype(str).str.strip().eq('').any(), f'{name}: blank ID')
    need(not frame.VariationID.duplicated().any(), f'{name}: duplicate ID')


def validate_scored(frame, name):
    flag = frame['plm_scored']
    score = pd.to_numeric(frame['esm1v_llr'], errors='raise')
    need(flag.isin([0, 1]).all(), f'{name}: nonbinary scored flag')
    need(not np.isinf(score.to_numpy(dtype=float)).any(), f'{name}: infinite ESM score')
    need((score.notna() == flag.eq(1)).all(), f'{name}: scored flag and score disagree')


def validate_cohort(data, expected=True):
    validate_ids(data, 'cohort')
    need(data.y.isin([0, 1]).all(), 'nonbinary cohort label')
    need(not data.gene_group.isna().any(), 'missing gene group')
    need(not data.gene_group.astype(str).str.strip().str.lower().isin(
        ['', 'na', 'nan', 'none', 'null']).any(), 'invalid gene group')
    validate_scored(data, 'cohort')
    if expected:
        m = data[data.plm_scored == 1]
        need(counts(data) == {'n': 3121, 'positive': 989, 'negative': 2132,
                             'gene_groups': 148}, 'mixed cohort denominator mismatch')
        need(counts(m) == {'n': 1328, 'positive': 925, 'negative': 403,
                          'gene_groups': 135}, 'missense denominator mismatch')


def load_data(matrix_path, plm_path, expected=True):
    matrix = features.load_matrix(matrix_path)
    validate_ids(matrix, 'matrix')
    pf = pd.read_csv(plm_path, float_precision='round_trip')
    validate_ids(pf, 'PLM features')
    validate_scored(pf, 'PLM features')
    data = plm.attach_plm_features(matrix, pf)
    data = data[(data.p4_headline_cohort == 1) &
                (data.consequence_class == 'coding_nontruncating')].copy()
    data['VariationID'] = data.VariationID.astype(str)
    data['gene_group'] = data.p4_gene_group
    data['y'] = data.y.astype(int)
    validate_cohort(data, expected)
    data['gene_group'] = data.gene_group.astype(str)
    return data.reset_index(drop=True)


def arms():
    base = list(protocol.ALL_FEATURES)
    return {'A': {'training': 'mixed', 'features': base},
            'B': {'training': 'mixed', 'features': base + ['esm1v_llr']},
            'C': {'training': 'scored_missense', 'features': base},
            'D': {'training': 'scored_missense', 'features': base + ['esm1v_llr']}}


def scientific_config(threads):
    need(threads > 0, 'threads must be positive')
    need(config.K_OUTER == 5 and config.SEEDS == list(range(101, 111)) and
         config.NUM_BOOST_ROUND == 700 and config.POSTHOC_MIN_ITER == 10 and
         config.INNER_VALID_FRAC == .2, 'frozen source config mismatch')
    params = config.xgb_params(101)
    params['nthread'] = threads
    return {'seeds': list(config.SEEDS), 'k_outer': config.K_OUTER,
            'num_boost_round': config.NUM_BOOST_ROUND,
            'posthoc_min_iter': config.POSTHOC_MIN_ITER,
            'inner_valid_fraction': config.INNER_VALID_FRAC,
            'target_recall': config.TARGET_RECALL, 'threads': threads,
            'xgb_params_seed101': params,
            'validation_population': 'scored_missense', 'refit': False}


def provenance(matrix_path, plm_path, threads):
    modules = {m.__name__: sha(m.__file__) for m in
               (config, protocol, pipeline, plm, features, metrics, fit_core)}
    modules['run_missense_2x2.py'] = sha(__file__)
    return {'inputs_sha256': {'matrix': sha(matrix_path), 'plm_features': sha(plm_path)},
            'code_sha256': modules,
            'versions': {'python': platform.python_version(), **{
                name: importlib.metadata.version(name) for name in
                ['numpy', 'pandas', 'xgboost', 'scikit-learn', 'scipy']}},
            'config': scientific_config(threads), 'arms': arms()}


def build_splits(data, seeds=None, k=None):
    """Split mixed cohort first; only then restrict training by arm."""
    seeds = config.SEEDS if seeds is None else seeds
    k = config.K_OUTER if k is None else k
    sizes = data.groupby('gene_group').size().to_dict()
    prevalence = data.groupby('gene_group').y.mean().to_dict()
    rows, outer_all, inner_all = [], {}, {}
    for seed in seeds:
        outer = pipeline.balanced_group_folds(sizes, prevalence, k, seed)
        outer_all[str(seed)] = outer
        for fold in range(k):
            train_groups = [g for g, f in outer.items() if f != fold]
            sz = {g: sizes[g] for g in train_groups}
            pr = {g: prevalence[g] for g in train_groups}
            inner_fold = pipeline.balanced_group_folds(
                sz, pr, max(2, round(1 / config.INNER_VALID_FRAC)), seed * 97 + fold)
            inner = {g: int(inner_fold[g] == 0) for g in train_groups}
            inner_all[f'{seed}:{fold}'] = inner
            te = data[data.gene_group.map(outer) == fold]
            tr = data[data.gene_group.map(outer) != fold]
            iv = tr[tr.gene_group.map(inner) == 1]
            it = tr[tr.gene_group.map(inner) == 0]
            gs = [set(x.gene_group) for x in [it, iv, te]]
            need(not (gs[0] & gs[1] or gs[0] & gs[2] or gs[1] & gs[2]),
                 'original mixed partition gene overlap')
            for arm in arms():
                used = {'inner_train': it if arm in 'AB' else it[it.plm_scored == 1],
                        'inner_val': iv[iv.plm_scored == 1],
                        'outer_test': te[te.plm_scored == 1]}
                for role, sub in used.items():
                    both(sub, f'{seed}/{fold}/{arm}/{role}')
                    part = sub[['VariationID', 'y', 'gene_group', 'plm_scored']].copy()
                    part['seed'], part['outer_fold'], part['arm'], part['role'] = seed, fold, arm, role
                    rows.append(part)
    split = pd.concat(rows, ignore_index=True)
    validate_splits(data, split, seeds, k)
    return split, outer_all, inner_all


def validate_splits(data, split, seeds, k):
    keys = ['seed', 'outer_fold', 'arm', 'role', 'VariationID']
    need(not split.duplicated(keys).any(), 'duplicate split membership')
    need(set(split.arm) == set(arms()), 'split arms mismatch')
    need(set(split.seed) == set(seeds), 'split seeds mismatch')
    need(set(split.outer_fold) == set(range(k)), 'split folds mismatch')
    need(set(split.role) == {'inner_train', 'inner_val', 'outer_test'}, 'split roles mismatch')
    identity = data.set_index('VariationID')
    need(split.VariationID.isin(identity.index).all(), 'unknown split ID')
    for col in ['y', 'gene_group', 'plm_scored']:
        wanted = split.VariationID.map(identity[col])
        need((wanted.to_numpy() == split[col].to_numpy()).all(), 'split identity mismatch: '+col)
    target = set(data.loc[data.plm_scored == 1, 'VariationID'])
    for seed in seeds:
        ss = split[split.seed == seed]
        for arm in arms():
            test = ss[(ss.arm == arm) & (ss.role == 'outer_test')]
            need(not test.VariationID.duplicated().any() and set(test.VariationID) == target,
                 'OOF target coverage mismatch')
        for fold in range(k):
            sf = ss[ss.outer_fold == fold]
            common = {}
            for arm in arms():
                sa = sf[sf.arm == arm]
                groups = []
                for role in ['inner_train', 'inner_val', 'outer_test']:
                    sub = sa[sa.role == role]
                    both(sub, f'validate {seed}/{fold}/{arm}/{role}')
                    groups.append(set(sub.gene_group))
                    if role != 'inner_train' or arm in 'CD':
                        need(sub.plm_scored.eq(1).all(), 'unscored restricted population')
                    if role != 'inner_train':
                        ids = list(sub.VariationID)
                        if role in common:
                            need(ids == common[role], 'cross-arm validation/test order differs')
                        common[role] = ids
                need(not (groups[0] & groups[1] or groups[0] & groups[2] or groups[1] & groups[2]),
                     'train/validation/test group overlap')
            for pair in [('A', 'B'), ('C', 'D')]:
                got = [list(sf[(sf.arm == a) & (sf.role == 'inner_train')].VariationID) for a in pair]
                need(got[0] == got[1], 'paired arm training identity mismatch')
            a_ids = set(sf[(sf.arm == 'A') & (sf.role == 'inner_train')].VariationID)
            c_ids = set(sf[(sf.arm == 'C') & (sf.role == 'inner_train')].VariationID)
            need(c_ids == {v for v in a_ids if identity.loc[v, 'plm_scored'] == 1},
                 'missense training is not exact mixed subset')


def prepare(args):
    out = args.out
    need(not out.exists() or not any(out.iterdir()), 'prepare needs new/empty output directory')
    out.mkdir(parents=True, exist_ok=True)
    data = load_data(args.matrix, args.plm_features)
    prov = provenance(args.matrix, args.plm_features, args.threads)
    split, outer, inner = build_splits(data)
    atomic_csv(out / 'split_manifest.csv', split)
    manifest = {**prov, 'protocol': 'ESM-MISSENSE-2X2-20261021', 'created_utc': utc(),
                'source_commit': 'd27ff142fed63a110114f8efa928cdd26bedd2df',
                'cohort_counts': {'mixed': counts(data), 'missense': counts(data[data.plm_scored == 1])},
                'outer_assignments': outer, 'inner_assignments': inner,
                'split_manifest_sha256': sha(out / 'split_manifest.csv')}
    atomic_json(out / 'frozen_manifest.json', manifest)
    print(json.dumps({'stage': 'prepared', 'planned_fits': 200, 'manifest_sha256': sha(out / 'frozen_manifest.json')}), flush=True)


@contextlib.contextmanager
def run_lock(out):
    with (out / 'fit.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


CHECKPOINT_FILES = {'test_predictions.csv', 'validation_predictions.csv', 'fold.json',
                    'validation_logloss.json', 'model.ubj'}


def validate_checkpoint(path, binding, identity):
    complete = json.loads((path / 'COMPLETE.json').read_text())
    need(complete['frozen_manifest_sha256'] == binding, 'checkpoint manifest mismatch')
    need(complete['identity'] == identity, 'checkpoint identity mismatch')
    need(set(complete['files_sha256']) == CHECKPOINT_FILES, 'incomplete checkpoint artifact list')
    for name, expected in complete['files_sha256'].items():
        need((path / name).is_file() and sha(path / name) == expected,
             'checkpoint hash mismatch: '+name)
    row = json.loads((path / 'fold.json').read_text())
    need({key: row[key] for key in identity} == identity, 'fold metadata mismatch')
    return row


def read_predictions(path):
    return pd.read_csv(path, dtype={'VariationID': str, 'gene_group': str},
                       float_precision='round_trip')


def predictions(frame, scores, threshold, identity):
    scores = np.asarray(scores, dtype=np.float64)
    need(len(scores) == len(frame) and np.isfinite(scores).all(), 'invalid prediction scores')
    pred = frame[['VariationID', 'y', 'gene_group']].copy()
    pred['score'] = scores
    pred['threshold'] = float(threshold)
    pred['decision'] = (scores >= threshold).astype(int)
    for key, value in identity.items():
        pred[key] = value
    return pred


def fit(args):
    out = args.out
    need((out / 'frozen_manifest.json').is_file(), 'missing prepared manifest')
    with run_lock(out):
        manifest = json.loads((out / 'frozen_manifest.json').read_text())
        binding = sha(out / 'frozen_manifest.json')
        for key, value in provenance(args.matrix, args.plm_features, args.threads).items():
            need(manifest[key] == value, 'frozen provenance differs: '+key)
        need(sha(out / 'split_manifest.csv') == manifest['split_manifest_sha256'], 'split hash mismatch')
        data = load_data(args.matrix, args.plm_features)
        split = pd.read_csv(out / 'split_manifest.csv', dtype={'VariationID': str, 'gene_group': str})
        validate_splits(data, split, config.SEEDS, config.K_OUTER)
        # Verify full assignments and order against source; fit STILL reads saved ID manifests.
        rebuilt, outer, inner = build_splits(data)
        need(split.equals(rebuilt), 'saved split does not match frozen source partition')
        need(outer == manifest['outer_assignments'] and inner == manifest['inner_assignments'],
             'frozen assignments differ from source')
        data = data.set_index('VariationID', drop=False)
        checkpoints = out / 'checkpoints'
        checkpoints.mkdir(exist_ok=True)
        all_test, fold_rows, cp_hashes = [], [], {}
        for seed in config.SEEDS:
            for fold in range(config.K_OUTER):
                for arm, armcfg in arms().items():
                    identity = {'seed': seed, 'outer_fold': fold, 'arm': arm}
                    name = f's{seed}_f{fold}_{arm}'
                    cp = checkpoints / name
                    sm = split[(split.seed == seed) & (split.outer_fold == fold) & (split.arm == arm)]
                    datasets = {role: data.loc[list(sm[sm.role == role].VariationID)].copy()
                                for role in ['inner_train', 'inner_val', 'outer_test']}
                    if cp.exists():
                        row = validate_checkpoint(cp, binding, identity)
                    else:
                        started = time.monotonic()
                        tmp = checkpoints / ('.tmp.' + name + '.' + uuid.uuid4().hex)
                        tmp.mkdir()
                        train, val, test = (datasets[r] for r in ['inner_train', 'inner_val', 'outer_test'])
                        fs = armcfg['features']
                        result = fit_core.fit_one(train[fs], train.y.to_numpy(), val[fs], val.y.to_numpy(),
                                                  test[fs], fs, plm.PLM_CONTINUOUS, seed, args.threads)
                        threshold = float(result['threshold'])
                        need(np.isfinite(threshold), 'invalid threshold')
                        atomic_csv(tmp / 'test_predictions.csv', predictions(test, result['test_scores'], threshold, identity))
                        atomic_csv(tmp / 'validation_predictions.csv', predictions(val, result['validation_scores'], threshold, identity))
                        result['model'].save_model(str(tmp / 'model.ubj'))
                        atomic_json(tmp / 'validation_logloss.json', result['validation_logloss'])
                        row = {**identity, 'roles_counts': {r: counts(d) for r, d in datasets.items()},
                               'best_iteration': int(result['best_iteration']), 'threshold': threshold,
                               'threshold_reason': result['threshold_reason'],
                               'imputation_medians': result['imputation_medians'], 'xgb_params': result['xgb_params'],
                               'features': fs, 'elapsed_seconds': time.monotonic() - started, 'completed_utc': utc()}
                        atomic_json(tmp / 'fold.json', row)
                        atomic_json(tmp / 'COMPLETE.json', {'identity': identity, 'frozen_manifest_sha256': binding,
                                    'files_sha256': {f: sha(tmp / f) for f in sorted(CHECKPOINT_FILES)}})
                        need(not cp.exists(), 'checkpoint unexpectedly appeared')
                        os.rename(tmp, cp)
                        row = validate_checkpoint(cp, binding, identity)
                    for filename, role in [('test_predictions.csv', 'outer_test'), ('validation_predictions.csv', 'inner_val')]:
                        pr = read_predictions(cp / filename)
                        ref = datasets[role]
                        need(list(pr.VariationID) == list(ref.VariationID), 'checkpoint IDs differ from manifest')
                        need(list(pr.y) == list(ref.y) and list(pr.gene_group) == list(ref.gene_group), 'checkpoint labels/groups mismatch')
                        need(np.isfinite(pr.score).all() and np.isfinite(pr.threshold).all(), 'nonfinite saved predictions')
                        need((pr.threshold == row['threshold']).all(), 'saved threshold mismatch')
                        need((pr.decision == (pr.score >= pr.threshold).astype(int)).all(), 'saved decisions inconsistent')
                        for key, value in identity.items():
                            need((pr[key] == value).all(), 'saved prediction identity mismatch')
                    all_test.append(read_predictions(cp / 'test_predictions.csv'))
                    fold_rows.append(row)
                    cp_hashes[name] = sha(cp / 'COMPLETE.json')
                    print(json.dumps({'stage': 'fit_complete', 'completed': len(fold_rows), 'total': 200}), flush=True)
        oof = pd.concat(all_test, ignore_index=True)
        need(len(oof) == 53120 and not oof.duplicated(['arm', 'seed', 'VariationID']).any(), 'final OOF rowcount/duplicate failure')
        for (_, _), group in oof.groupby(['arm', 'seed']):
            need(len(group) == 1328 and int(group.y.sum()) == 925, 'final OOF denominator failure')
        pf = pd.DataFrame([{k: v for k, v in row.items() if k not in {'roles_counts', 'imputation_medians', 'xgb_params', 'features'}}
                           for row in fold_rows])
        seal_path = out / 'prediction_seal.json'
        if seal_path.exists():
            seal = json.loads(seal_path.read_text())
            need(set(seal['files_sha256']) == {'oof_predictions.csv', 'per_fold.csv'} and seal['total_fits'] == 200 and seal['oof_rows'] == 53120, 'invalid final seal schema')
            need(seal['frozen_manifest_sha256'] == binding and seal['checkpoint_complete_sha256'] == cp_hashes,
                 'final seal binding mismatch')
            for name, value in seal['files_sha256'].items():
                need(sha(out / name) == value, 'final artifact hash mismatch')
        else:
            atomic_csv(out / 'oof_predictions.csv', oof.rename(columns={'gene_group': 'group', 'outer_fold': 'fold'}))
            atomic_csv(out / 'per_fold.csv', pf)
            atomic_json(seal_path, {'generated_utc': utc(), 'total_fits': 200, 'oof_rows': len(oof),
                        'frozen_manifest_sha256': binding, 'checkpoint_complete_sha256': cp_hashes,
                        'files_sha256': {n: sha(out / n) for n in ['oof_predictions.csv', 'per_fold.csv']}})
        print(json.dumps({'stage': 'sealed', 'total_fits': 200, 'seal_sha256': sha(seal_path)}), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['prepare', 'fit'])
    ap.add_argument('--matrix', type=Path, required=True)
    ap.add_argument('--plm-features', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--repo-root', type=Path)
    ap.add_argument('--threads', type=int, default=10)
    args = ap.parse_args()
    try:
        libs(args.repo_root)
        (prepare if args.command == 'prepare' else fit)(args)
    except Exception:
        if args.out.is_dir():
            atomic_json(args.out / ('failure.' + uuid.uuid4().hex + '.json'),
                        {'utc': utc(), 'command': args.command, 'traceback': traceback.format_exc()})
        raise

if __name__ == '__main__':
    main()
