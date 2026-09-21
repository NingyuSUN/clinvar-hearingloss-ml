"""Deterministic oracles for grouped inference and repeated-OOF validation."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hlpath.metrics import bootstrap_ci, paired_delta
from hlpath.analysis import _pooled


def test_unequal_clusters_match_explicit_resampling_oracle():
    y = np.array([0, 1, 1, 0, 1, 0, 0])
    score = np.array([0.2, 0.2, 0.8, 0.4, 0.9, 0.6, 0.1])
    groups = np.array(["a", "a", "b", "c", "c", "c", "c"])
    pred = (score >= 0.5).astype(int)
    actual = bootstrap_ci(y, score, pred, n_boot=500, seed=17, groups=groups)
    rng = np.random.default_rng(17)
    members = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    auc, recall, precision = [], [], []
    for _ in range(500):
        ix = np.concatenate([members[g] for g in rng.integers(0, 3, 3)])
        yy, pp = y[ix], pred[ix]
        if np.unique(yy).size == 2:
            auc.append(roc_auc_score(yy, score[ix]))
        if yy.sum():
            recall.append((yy * pp).sum() / yy.sum())
        if pp.sum():
            precision.append((yy * pp).sum() / pp.sum())
    assert actual["AUC"]["estimate"] == pytest.approx(roc_auc_score(y, score))
    for key, values in [("AUC", auc), ("R90_recall", recall), ("R90_precision", precision)]:
        assert actual[key]["n_valid"] == len(values)
        assert actual[key]["bootstrap_mean"] == pytest.approx(np.mean(values))
        assert [actual[key]["lo"], actual[key]["hi"]] == pytest.approx(np.percentile(values, [2.5, 97.5]))


def test_undefined_draws_counted_per_metric():
    out = bootstrap_ci([0, 1], [0.1, 0.9], [0, 0], groups=["a", "b"], n_boot=100)
    assert out["R90_precision"]["estimate"] is None
    assert out["R90_precision"]["n_valid"] == 0
    assert out["AUC"]["n_valid"] < out["R90_recall"]["n_valid"]
    assert out["AUC"]["lo"] is None
    assert out["AUC"]["status"] == "insufficient_valid_draws"


@pytest.mark.parametrize("change", ["missing_group", "row_group", "one_group", "bad_label", "infinite", "bad_pred"])
def test_bootstrap_rejects_invalid_input(change):
    y, p, g, pred = [0, 1], [0.1, 0.9], ["a", "b"], [0, 1]
    if change == "missing_group": g = [None, "b"]
    if change == "row_group": g = ["row0", "row1"]
    if change == "one_group": g = ["a", "a"]
    if change == "bad_label": y = [0, 2]
    if change == "infinite": p = [float("inf"), 0.9]
    if change == "bad_pred": pred = [0, 2]
    with pytest.raises(ValueError): bootstrap_ci(y, p, pred, groups=g, n_boot=100)


def fixture_oof():
    pf = pd.DataFrame({"seed": [1, 2], "outer_fold": [0, 0], "threshold": [0.5, 0.5], "threshold_ok": [True, True]})
    oof = pd.DataFrame({"seed": [1, 1, 2, 2], "outer_fold": [0]*4, "VariationID": [10, 11, 10, 11],
                        "y": [0, 1, 0, 1], "prob": [0.2, 0.7, 0.6, 0.9],
                        "gene_group": ["a", "b", "a", "b"], "consequence_class": ["coding"]*4})
    return pf, oof


def test_oof_averages_seeds_and_keeps_tie_positive_contract():
    pf, oof = fixture_oof()
    pooled = _pooled(pf, oof)
    assert pooled["prob"].tolist() == pytest.approx([0.4, 0.8])
    assert pooled["hit"].tolist() == [0.5, 1.0]
    assert pooled["gene_group"].tolist() == ["a", "b"]


@pytest.mark.parametrize("change", ["duplicate_prediction", "missing_seed", "label_conflict", "group_conflict", "missing_threshold", "duplicate_threshold", "failed_threshold", "missing_fold"])
def test_oof_fails_closed(change):
    pf, oof = fixture_oof()
    if change == "duplicate_prediction": oof = pd.concat([oof, oof.iloc[:1]])
    if change == "missing_seed": oof = oof.iloc[:-1]
    if change == "label_conflict": oof.loc[2, "y"] = 1
    if change == "group_conflict": oof.loc[2, "gene_group"] = "c"
    if change == "missing_threshold": pf.loc[0, "threshold"] = np.nan
    if change == "duplicate_threshold": pf = pd.concat([pf, pf.iloc[:1]])
    if change == "failed_threshold": pf.loc[0, "threshold_ok"] = False
    if change == "missing_fold": pf = pf.iloc[:1]
    with pytest.raises(ValueError): _pooled(pf, oof)


def test_paired_differences_are_descriptive_and_require_matching_folds():
    pf = pd.DataFrame({"seed": [1, 1, 1, 1], "outer_fold": [0, 1, 0, 1],
                       "feature_set": ["a", "a", "b", "b"], "auc": [0.5, 0.6, 0.7, 0.7]})
    result = paired_delta(pf, "a", "b")
    assert result["mean_delta"] == pytest.approx(0.15)
    assert result["inference"] == "descriptive_only"
    assert not {"paired_se", "paired_t", "significant_p05"} & result.keys()
    with pytest.raises(ValueError): paired_delta(pf.iloc[:-1], "a", "b")
    with pytest.raises(ValueError): paired_delta(pd.concat([pf, pf.iloc[:1]]), "a", "b")


def test_oof_checks_recorded_test_sizes():
    pf, oof = fixture_oof()
    pf["n_test"] = [2, 3]
    with pytest.raises(ValueError, match="test sizes"):
        _pooled(pf, oof)


def test_random_split_retains_biological_groups(monkeypatch):
    # Fake the booster; exercise the real OOF-output construction without fitting.
    import types
    fake_xgb = types.ModuleType("xgboost")
    fake_xgb.DMatrix = lambda x, **kw: x
    class Booster:
        def predict(self, x, **kw):
            return np.linspace(0.1, 0.9, len(x))
    def train(*args, **kwargs):
        kwargs["evals_result"]["iv"] = {"logloss": [0.5]}
        return Booster()
    fake_xgb.train = train
    from hlpath import pipeline
    monkeypatch.setattr(pipeline, "xgb", fake_xgb)
    data = pd.DataFrame({"VariationID": range(6), "y": [0, 1, 0, 1, 0, 1],
                         "x": range(6), "GeneSymbol": ["A"]*6, "p4_gene_group": ["GG01"]*6,
                         "_grp": [f"row{i}" for i in range(6)], "consequence_class": ["coding"]*6,
                         "p4_review_star": [1]*6})
    outer = {f"row{i}": 0 if i < 2 else 1 for i in range(6)}
    inner = {f"row{i}": 1 if i < 4 else 0 for i in range(2, 6)}
    _, oof = pipeline._one_fold(data, ["x"], outer, inner, 0, 101)
    assert oof["gene_group"].tolist() == ["GG01", "GG01"]
    assert oof["split_group"].tolist() == ["row0", "row1"]
