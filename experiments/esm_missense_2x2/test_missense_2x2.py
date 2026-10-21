"""Synthetic engineering tests; no real variants are used."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_missense_2x2 as r
r.libs()
import numpy as np
import pandas as pd
import fit_core as core


def synthetic():
    rows = []
    for gene in range(40):
        for i in range(6):
            scored = int(i < 4)
            rows.append({'VariationID': f'v{gene:03d}_{i}', 'y': i % 2,
                         'gene_group': f'G{gene:03d}', 'plm_scored': scored,
                         'esm1v_llr': float(i) if scored else np.nan})
    return pd.DataFrame(rows)


class SplitTests(unittest.TestCase):
    def test_common_targets_subset_training_and_gene_isolation(self):
        data = synthetic()
        split, outer, inner = r.build_splits(data, seeds=[101, 102], k=5)
        target = set(data.loc[data.plm_scored == 1, 'VariationID'])
        for seed in [101, 102]:
            for arm in 'ABCD':
                ss = split[(split.seed == seed) & (split.arm == arm)]
                test = ss[ss.role == 'outer_test']
                self.assertEqual(set(test.VariationID), target)
                self.assertFalse(test.VariationID.duplicated().any())
                for fold in range(5):
                    sf = ss[ss.outer_fold == fold]
                    train = set(sf[sf.role == 'inner_train'].gene_group)
                    val = set(sf[sf.role == 'inner_val'].gene_group)
                    heldout = set(sf[sf.role == 'outer_test'].gene_group)
                    self.assertFalse(train & val or train & heldout or val & heldout)
                    # Every held-out group must follow original mixed fold assignment.
                    self.assertTrue(all(outer[str(seed)][g] == fold for g in heldout))
                    self.assertTrue(all(inner[f'{seed}:{fold}'][g] == 1 for g in val))
        for (seed, fold), sf in split.groupby(['seed', 'outer_fold']):
            a = sf[(sf.arm == 'A') & (sf.role == 'inner_train')]
            c = sf[(sf.arm == 'C') & (sf.role == 'inner_train')]
            self.assertEqual(list(c.VariationID), list(a.loc[a.plm_scored == 1, 'VariationID']))
            for role in ['inner_val', 'outer_test']:
                vals = [list(sf[(sf.arm == arm) & (sf.role == role)].VariationID) for arm in 'ABCD']
                self.assertTrue(all(x == vals[0] for x in vals))

    def test_missing_gene_and_duplicate_id_fail(self):
        data = synthetic()
        r.validate_cohort(data, expected=False)
        bad = data.copy(); bad.loc[0, 'gene_group'] = None
        with self.assertRaises(ValueError):
            r.validate_cohort(bad, expected=False)
        bad = data.copy(); bad.loc[1, 'VariationID'] = bad.loc[0, 'VariationID']
        with self.assertRaises(ValueError):
            r.validate_cohort(bad, expected=False)

    def test_score_flag_mismatch_fails(self):
        data = synthetic()
        data.loc[0, 'esm1v_llr'] = np.nan
        with self.assertRaises(ValueError):
            r.validate_scored(data, 'test')
        data = synthetic(); data.loc[4, 'esm1v_llr'] = .2
        with self.assertRaises(ValueError):
            r.validate_scored(data, 'test')

    def test_altered_manifest_identity_or_overlap_fails(self):
        data = synthetic()
        split, _, _ = r.build_splits(data, seeds=[101], k=5)
        bad = split.copy(); bad.loc[0, 'y'] = 1 - bad.loc[0, 'y']
        with self.assertRaises(ValueError):
            r.validate_splits(data, bad, [101], 5)
        bad = pd.concat([split, split.iloc[[0]]], ignore_index=True)
        with self.assertRaises(ValueError):
            r.validate_splits(data, bad, [101], 5)


class CoreTests(unittest.TestCase):
    def test_training_only_imputation_and_binary_zero(self):
        train = pd.DataFrame({'continuous': [0., 2., np.nan], 'binary': [1., 0., np.nan]})
        val = pd.DataFrame({'continuous': [10000., np.nan], 'binary': [np.nan, 1.]})
        (fixed_train, fixed_val), med = core.impute_train_only(train, val, continuous={'continuous'})
        self.assertEqual(med, {'continuous': 1.})
        self.assertEqual(fixed_val.loc[1, 'continuous'], 1.)
        self.assertEqual(fixed_train.loc[2, 'binary'], 0.)
        self.assertTrue(np.isnan(train.loc[2, 'continuous']))

    def test_all_missing_continuous_and_infinite_fail(self):
        with self.assertRaises(ValueError):
            core.impute_train_only(pd.DataFrame({'x': [np.nan, np.nan]}), continuous={'x'})
        with self.assertRaises(ValueError):
            core.impute_train_only(pd.DataFrame({'x': [1., np.inf]}), continuous={'x'})

    def test_best_iteration_ignores_initial_rounds(self):
        curve = [0.] * 10 + [2., .4, .4, 1.]
        self.assertEqual(core.select_iteration(curve), 11)
        with self.assertRaises(ValueError):
            core.select_iteration([.1] * 10)

    def test_fit_uses_only_train_and_validation_labels(self):
        captured = []
        class Matrix:
            def __init__(self, data, label=None, feature_names=None):
                self.data = np.asarray(data); self.label = label
                captured.append(self)
        class Booster:
            def predict(self, mat, iteration_range):
                self_range = (0, 24)
                if iteration_range != self_range:
                    raise AssertionError('unexpected best iteration')
                return np.asarray([.9, .8, .2, .1] if mat.label is not None else [.8, .01], dtype=np.float32)
        def train(params, dtrain, num_boost_round, evals, evals_result, verbose_eval):
            self.assertEqual(num_boost_round, 700)
            self.assertEqual(params['nthread'], 10)
            self.assertEqual(len(evals), 1)
            self.assertEqual(list(dtrain.data[:, 0]), [1., 2., 3., 4.])
            self.assertEqual(list(evals[0][0].data[:, 0]), [10., 11., 12., 13.])
            self.assertFalse(verbose_eval)
            curve = [1.] * 700; curve[0] = 0.; curve[23] = .2
            evals_result['validation'] = {'logloss': curve}
            return Booster()
        with patch.object(core.xgb, 'DMatrix', Matrix), patch.object(core.xgb, 'train', train):
            got = core.fit_one(pd.DataFrame({'x': [1., 2., 3., 4.]}), [0, 1, 0, 1],
                               pd.DataFrame({'x': [10., 11., 12., 13.]}), [1, 1, 0, 0],
                               pd.DataFrame({'x': [90., 91.]}), ['x'], {'x'}, 101, 10)
        self.assertIsNone(captured[2].label)
        self.assertEqual(got['best_iteration'], 23)
        self.assertEqual(got['test_scores'].dtype, np.float64)
        self.assertEqual(got['threshold'], got['test_scores'][0])
        self.assertEqual(got['imputation_medians'], {'x': 2.5})


class ArtifactTests(unittest.TestCase):
    def test_float32_threshold_tie_survives_csv(self):
        threshold = float(np.float32(.1))
        frame = pd.DataFrame({'VariationID': ['a', 'b'], 'y': [1, 0], 'gene_group': ['G1', 'G2']})
        score = np.array([np.float32(.1), np.nextafter(np.float32(.1), np.float32(0))])
        got = r.predictions(frame, score, threshold, {'seed': 101, 'outer_fold': 0, 'arm': 'D'})
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'p.csv'; r.atomic_csv(p, got)
            loaded = r.read_predictions(p)
        self.assertEqual(list(loaded.decision), [1, 0])
        self.assertEqual(loaded.score.iloc[0], loaded.threshold.iloc[0])
        self.assertTrue((loaded.decision == (loaded.score >= loaded.threshold).astype(int)).all())

    def checkpoint(self, path):
        identity = {'seed': 101, 'outer_fold': 0, 'arm': 'D'}
        for f in r.CHECKPOINT_FILES:
            (path / f).write_text('synthetic bytes')
        r.atomic_json(path / 'fold.json', identity)
        r.atomic_json(path / 'COMPLETE.json', {'identity': identity, 'frozen_manifest_sha256': 'frozen',
                      'files_sha256': {f: r.sha(path / f) for f in r.CHECKPOINT_FILES}})
        return identity

    def test_resume_rejects_corrupt_or_different_checkpoint(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d); ident = self.checkpoint(path)
            r.validate_checkpoint(path, 'frozen', ident)
            with self.assertRaises(ValueError):
                r.validate_checkpoint(path, 'different manifest', ident)
            with self.assertRaises(ValueError):
                r.validate_checkpoint(path, 'frozen', dict(ident, arm='C'))
            (path / 'model.ubj').write_text('corrupt')
            with self.assertRaises(ValueError):
                r.validate_checkpoint(path, 'frozen', ident)

    def test_resume_requires_every_artifact(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d); ident = self.checkpoint(path)
            seal = json.loads((path / 'COMPLETE.json').read_text())
            del seal['files_sha256']['validation_predictions.csv']
            r.atomic_json(path / 'COMPLETE.json', seal)
            with self.assertRaises(ValueError):
                r.validate_checkpoint(path, 'frozen', ident)

if __name__ == '__main__':
    unittest.main()
