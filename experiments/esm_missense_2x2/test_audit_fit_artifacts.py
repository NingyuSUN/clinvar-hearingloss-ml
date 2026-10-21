"""Synthetic checks for the independent fit-artifact audit."""
import unittest

import numpy as np
import pandas as pd

from audit_fit_artifacts import (compare_split_tables, expected_medians, imputed_matrix,
                                 independent_r90, reconstruct_splits, validate_prediction_table)


class IndependentAuditTests(unittest.TestCase):
    def test_positive_quantile_matches_exhaustive_r90_with_ties(self):
        rng = np.random.default_rng(17)
        for n_pos in (7, 10, 11, 20, 149):
            for _ in range(10):
                y = np.array([1] * n_pos + [0] * 17)
                scores = rng.integers(0, 12, len(y)) / 12
                choices = []
                for threshold in np.unique(scores):
                    predicted = scores >= threshold
                    tp = int(((y == 1) & predicted).sum())
                    fp = int(((y == 0) & predicted).sum())
                    if tp / n_pos >= 0.9:
                        choices.append((fp, -threshold, threshold))
                self.assertEqual(independent_r90(y, scores), min(choices)[2])

    def test_negative_scores_do_not_change_fixed_positive_quantile(self):
        y = np.array([1] * 10 + [0] * 5)
        fixed = np.linspace(0.1, 1, 10)
        low = independent_r90(y, np.r_[fixed, np.zeros(5)])
        high = independent_r90(y, np.r_[fixed, np.ones(5)])
        self.assertEqual(low, high)

    def test_training_medians_and_binary_zero_ignore_validation_extremes(self):
        train = pd.DataFrame({"continuous": [1, 3, np.nan], "binary": [1, np.nan, 0]})
        validation = pd.DataFrame({"continuous": [np.nan, 10000], "binary": [np.nan, 1]})
        med = expected_medians(train, ["continuous", "binary"], {"continuous"})
        self.assertEqual(med, {"continuous": 2.0})
        np.testing.assert_array_equal(imputed_matrix(validation, ["continuous", "binary"], med),
                                      np.array([[2, 0], [10000, 1]]))

    def test_reconstruction_group_isolation_and_missense_filter_after_splitting(self):
        rows = []
        for group in range(10):
            for record in range(3):
                rows.append({"VariationID": f"synthetic_{group}_{record}", "gene_group": f"g{group}",
                             "y": int(record == 1), "plm_scored": int(record != 2)})
        data = pd.DataFrame(rows)
        # A deterministic substitute for the base allocator tests audit control
        # flow; the live audit separately binds the actual base Git function.
        def assign(sizes, prevalence, k, seed):
            return {g: (i + seed) % k for i, g in enumerate(sorted(sizes))}
        roles, split, outer, inner = reconstruct_splits(data, assign)
        self.assertEqual(len(roles), 200)
        self.assertEqual(len(outer), 10)
        self.assertEqual(len(inner), 50)
        for seed, fold, arm in roles:
            used = roles[seed, fold, arm]
            groups = [set(frame.gene_group) for frame in used.values()]
            self.assertFalse(groups[0] & groups[1] or groups[0] & groups[2] or groups[1] & groups[2])
            self.assertTrue(used["inner_val"].plm_scored.eq(1).all())
            self.assertTrue(used["outer_test"].plm_scored.eq(1).all())
            if arm in "CD":
                expected = roles[seed, fold, "A"]["inner_train"]
                self.assertEqual(list(used["inner_train"].VariationID),
                                 list(expected.loc[expected.plm_scored == 1, "VariationID"]))
        compare_split_tables(split.copy(), split)
        changed = split.copy()
        changed.loc[0, "role"] = "outer_test"
        with self.assertRaises((AssertionError, ValueError)):
            compare_split_tables(changed, split)

    def test_prediction_table_rejects_changed_identity_and_threshold_decisions(self):
        reference = pd.DataFrame({"VariationID": ["s1", "s2"], "y": [1, 0], "gene_group": ["g1", "g2"]})
        identity = {"seed": 101, "outer_fold": 0, "arm": "C"}
        frame = reference.assign(score=[0.5, 0.1], threshold=0.5, decision=[1, 0], **identity)
        validate_prediction_table(frame, reference, identity, 0.5)
        corrupted = frame.copy()
        corrupted.loc[0, "decision"] = 0
        with self.assertRaisesRegex(ValueError, "decision mismatch"):
            validate_prediction_table(corrupted, reference, identity, 0.5)
        corrupted = frame.copy()
        corrupted.loc[1, "gene_group"] = "wrong"
        with self.assertRaisesRegex(ValueError, "alignment mismatch"):
            validate_prediction_table(corrupted, reference, identity, 0.5)


if __name__ == "__main__":
    unittest.main(verbosity=2)
