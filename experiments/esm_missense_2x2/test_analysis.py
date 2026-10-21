"""Synthetic tests only; no real variant table is needed."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from analysis import (ARMS, SEEDS, aggregate, bootstrap, discordance, metrics,
                      validate_rows, weighted_auc, write_csv)


def synthetic_rows():
    rows = []
    for arm in ARMS:
        for seed in SEEDS:
            for i in range(10):
                y = i % 2
                score = 0.8 if y else 0.2
                if arm == "C" and i == 1:
                    score = 0.2
                if arm == "D" and i == 3:
                    score = 0.2
                # A 5/5 vote tie for this benign item must classify positive.
                if i == 0 and seed <= 105:
                    score = 0.8
                rows.append(dict(VariationID=f"synthetic_{i}", y=y, group=f"g{i//2}",
                                 score=score, threshold=0.5, decision=int(score >= 0.5),
                                 seed=seed, fold=i//2, arm=arm))
    return rows


def panel():
    return validate_rows(synthetic_rows(), 10, 5, 5)


class AnalysisTests(unittest.TestCase):
    def test_auc_ties(self):
        self.assertEqual(weighted_auc([0, 1], [0.5, 0.5]), 0.5)
        self.assertEqual(weighted_auc([0, 1], [0.1, 0.9]), 1.0)
        self.assertEqual(weighted_auc([0, 1], [0.9, 0.1]), 0.0)
        self.assertIsNone(weighted_auc([0, 1], [0.1, 0.9], [1, 0]))

    def test_cluster_multiplicity_matches_explicit_replication(self):
        y = np.array([1, 0, 1, 0, 1, 0])
        score = np.array([0.6, 0.6, 0.8, 0.2, 0.3, 0.4])
        decision = (score >= 0.5).astype(int)
        # Groups (0,1), (2,3), (4,5) drawn twice, zero, once.
        weight = np.array([2, 2, 0, 0, 1, 1])
        repeated = np.repeat(np.arange(6), weight)
        self.assertEqual(metrics(y, score, decision, weight),
                         metrics(y[repeated], score[repeated], decision[repeated]))

    def test_duplicate_input_rejected(self):
        rows = synthetic_rows()
        rows.append(rows[0].copy())
        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_rows(rows, 10, 5, 5)

    def test_missing_input_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            validate_rows(synthetic_rows()[:-1], 10, 5, 5)

    def test_truth_and_group_inconsistency_rejected(self):
        rows = synthetic_rows()
        rows[10]["y"] = 1
        with self.assertRaisesRegex(ValueError, "label/group"):
            validate_rows(rows, 10, 5, 5)

    def test_cross_arm_fold_mismatch_rejected(self):
        rows = synthetic_rows()
        rows[100]["fold"] = 1
        with self.assertRaisesRegex(ValueError, "outer fold"):
            validate_rows(rows, 10, 5, 5)

    def test_group_leakage_rejected(self):
        rows = synthetic_rows()
        rows[1]["fold"] = 1
        with self.assertRaisesRegex(ValueError, "group split"):
            validate_rows(rows, 10, 5, 5)

    def test_threshold_inconsistency_rejected(self):
        rows = synthetic_rows()
        rows[1]["threshold"] = 0.6
        with self.assertRaisesRegex(ValueError, "threshold differs"):
            validate_rows(rows, 10, 5, 5)

    def test_decision_threshold_mismatch_rejected(self):
        rows = synthetic_rows()
        rows[0]["decision"] = 0
        with self.assertRaisesRegex(ValueError, "decision !="):
            validate_rows(rows, 10, 5, 5)

    def test_invalid_threshold_and_missing_group_rejected(self):
        rows = synthetic_rows()
        rows[0]["threshold"] = -0.1
        with self.assertRaisesRegex(ValueError, "invalid threshold"):
            validate_rows(rows, 10, 5, 5)
        rows = synthetic_rows()
        rows[0]["group"] = "nan"
        with self.assertRaisesRegex(ValueError, "empty ID/group"):
            validate_rows(rows, 10, 5, 5)

    def test_truncated_row_whitespace_and_single_class_expectation_rejected(self):
        rows = synthetic_rows()
        rows[0]["fold"] = None
        with self.assertRaisesRegex(ValueError, "missing required"):
            validate_rows(rows, 10, 5, 5)
        rows = synthetic_rows()
        rows[0]["group"] = "g0 "
        with self.assertRaisesRegex(ValueError, "whitespace"):
            validate_rows(rows, 10, 5, 5)
        rows[0]["group"] = "g0"
        rows[0]["fold"] = "NA"
        with self.assertRaisesRegex(ValueError, "Line 2: fold must be an integer"):
            validate_rows(rows, 10, 5, 5)
        with self.assertRaisesRegex(ValueError, "both expected classes"):
            validate_rows(synthetic_rows(), 10, 10, 5)

    def test_aggregation_vote_tie_and_paired_fn_flips(self):
        p = panel()
        scores, decisions, summary = aggregate(p)
        self.assertAlmostEqual(scores["A"][0], 0.5)
        self.assertEqual(decisions["A"][0], 1)
        counts, changed = discordance(p, scores, decisions)
        self.assertEqual(counts["pathogenic"]["C0_D1"], 1)
        self.assertEqual(counts["pathogenic"]["C1_D0"], 1)
        self.assertEqual(summary["C"]["FN"], 1)
        self.assertEqual(summary["D"]["FN"], 1)
        self.assertEqual(len(changed), 2)
        self.assertEqual({r["change"] for r in changed}, {"corrected", "new_error"})

    def test_shared_draws_and_reproducibility(self):
        p = panel()
        scores, decisions, summary = aggregate(p)
        c1, meta1 = bootstrap(p, scores, decisions, summary, 200, 0)
        c2, meta2 = bootstrap(p, scores, decisions, summary, 200, 0)
        self.assertEqual(c1, c2)
        self.assertEqual(meta1, meta2)
        # Identical arms A/B yield exactly zero every paired interval.
        for row in c1:
            if row["contrast"] == "B-A":
                self.assertEqual((row["difference"], row["ci95_low"], row["ci95_high"]), (0, 0, 0))

    def test_intervals_withheld_with_insufficient_draws(self):
        p = panel()
        scores, decisions, summary = aggregate(p)
        contrasts, meta = bootstrap(p, scores, decisions, summary, 20, 0)
        self.assertFalse(meta["intervals_available"])
        self.assertTrue(all(r["ci95_low"] is None and r["ci95_high"] is None for r in contrasts))

    def test_actual_single_class_group_draws_withhold_intervals(self):
        p = panel()
        # Each cluster is now pure-class. Drawing two clusters with replacement
        # yields one-class data for approximately half the draws, not zero.
        p["groups"] = np.where(p["y"] == 1, "positive_cluster", "negative_cluster")
        scores, decisions, summary = aggregate(p)
        contrasts, meta = bootstrap(p, scores, decisions, summary, 200, 0)
        self.assertGreater(meta["invalid_single_class_draws"], 20)
        self.assertEqual(meta["valid_draws"] + meta["invalid_single_class_draws"], 200)
        self.assertFalse(meta["intervals_available"])
        self.assertTrue(all(r["difference"] is not None and r["ci95_low"] is None for r in contrasts))

    def test_csv_explicit_order_and_empty_schema(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.csv"
            write_csv(path, [{"a": 1, "b": 2}], ["b", "a"])
            self.assertEqual(path.read_text().splitlines()[0], "b,a")
            write_csv(path, [], ["a", "b"])
            self.assertEqual(path.read_text(), "a,b\n")
            with self.assertRaisesRegex(ValueError, "empty CSV without fieldnames"):
                write_csv(path, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
