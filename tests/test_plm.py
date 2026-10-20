"""Deterministic checks for the PLM-ESM-2026-10 feature rules and comparisons."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hlpath import plm  # noqa: E402
from hlpath.compare import paired_group_bootstrap_auc  # noqa: E402
from hlpath.pipeline import _impute, run_experiment  # noqa: E402
from hlpath.protocol import ALL_FEATURES, FEATURE_SETS  # noqa: E402


# --- windows and scores ------------------------------------------------------
@pytest.mark.parametrize("length,position,expected", [
    (300, 1, (0, 300)), (1022, 1022, (0, 1022)),
    (2000, 1, (0, 1022)), (2000, 600, (88, 1110)), (2000, 2000, (978, 2000)),
    (1023, 512, (0, 1022)), (1023, 513, (1, 1023)),
])
def test_esm_window(length, position, expected):
    start, end = plm.esm_window(length, position)
    assert (start, end) == expected
    assert start < position <= end and end - start <= plm.MAX_RESIDUES


@pytest.mark.parametrize("position", [0, 301])
def test_esm_window_rejects_outside_positions(position):
    with pytest.raises(ValueError):
        plm.esm_window(300, position)


def test_llr_sign_and_finiteness():
    lp = {a: -3.0 for a in plm.AA20} | {"R": -0.5, "W": -6.0}
    assert plm.masked_marginal_llr(lp, "R", "W") == pytest.approx(-5.5)
    with pytest.raises(ValueError):
        plm.masked_marginal_llr(lp | {"W": float("nan")}, "R", "W")


def test_vcf_input():
    assert plm.vcf_input("16-2496587-G-C") == "16 2496587 . G C . . ."
    with pytest.raises(ValueError):
        plm.vcf_input("16-x-G-C")


# --- mapping rule ----------------------------------------------------------
def entry(**kw):
    base = {"transcript_id": "ENST1.2", "protein_id": "ENSP1", "consequence_terms": ["missense_variant"],
            "amino_acids": "R/W", "protein_start": 3, "protein_end": 3}
    return base | kw


@pytest.mark.parametrize("kw,seq,status", [
    ({}, "MAR", plm.STATUS_SCORED),
    ({}, "MAK", plm.STATUS_REF_MISMATCH),
    ({}, "MA", plm.STATUS_REF_MISMATCH),
    ({}, None, plm.STATUS_NO_SEQUENCE),
    ({"consequence_terms": ["synonymous_variant"], "amino_acids": "R"}, "MAR", plm.STATUS_NOT_MISSENSE),
    ({"amino_acids": "RG/WS", "protein_end": 4}, "MARG", plm.STATUS_MULTI_RESIDUE),
    ({"amino_acids": "R/X"}, "MAR", plm.STATUS_NONSTANDARD),
    ({"consequence_terms": ["missense_variant", "splice_region_variant"]}, "MAR", plm.STATUS_SCORED),
])
def test_map_consequence(kw, seq, status):
    assert plm.map_consequence(entry(**kw), seq)["map_status"] == status


def test_map_consequence_absent_transcript():
    assert plm.map_consequence(None, "MAR")["map_status"] == plm.STATUS_TRANSCRIPT_ABSENT


def test_select_transcript_by_stable_id_ignores_version():
    rec = {"transcript_consequences": [entry(transcript_id="ENST9.1"), entry(transcript_id="ENST1.5")]}
    assert plm.select_transcript_consequence(rec, "ENST1.2")["transcript_id"] == "ENST1.5"
    assert plm.select_transcript_consequence(rec, "ENST7.1") is None
    assert plm.select_transcript_consequence(None, "ENST1.2") is None
    dup = {"transcript_consequences": [entry(), entry()]}
    with pytest.raises(ValueError):
        plm.select_transcript_consequence(dup, "ENST1")


# --- features ------------------------------------------------------------
def logprob_rows(tid, pos, models, lp_ref, lp_alt, ref="R", alt="W"):
    rows = []
    for m, r, a in zip(models, lp_ref, lp_alt):
        lp = {f"lp_{x}": -4.0 for x in plm.AA20} | {f"lp_{ref}": r, f"lp_{alt}": a}
        rows.append({"transcript_stable_id": tid, "protein_position": pos, "model": m} | lp)
    return rows


def mapping_frame():
    return pd.DataFrame([
        {"VariationID": 1, "map_status": "scored", "cached_transcript_id": "ENST1.2",
         "protein_position": 3, "ref_aa": "R", "alt_aa": "W"},
        {"VariationID": 2, "map_status": "not_missense", "cached_transcript_id": "ENST1.2",
         "protein_position": np.nan, "ref_aa": None, "alt_aa": None},
        {"VariationID": 3, "map_status": "scored", "cached_transcript_id": "ENST1.2",
         "protein_position": 3, "ref_aa": "R", "alt_aa": "C"},
    ])


def test_variant_llrs_ensemble_mean_and_unscored():
    models = [*plm.ESM1V_MODELS, plm.ESM2_MODEL]
    lp = pd.DataFrame(logprob_rows("ENST1", 3, models, [-1, -1, -1, -1, -1, -2], [-2, -3, -4, -5, -6, -7]))
    out = plm.variant_llrs(mapping_frame(), lp).set_index("VariationID")
    assert out.loc[1, "esm1v_llr"] == pytest.approx(np.mean([-1, -2, -3, -4, -5]))
    assert out.loc[1, "esm2_llr"] == pytest.approx(-5)
    assert out.loc[3, "esm1v_llr"] == pytest.approx(-3.0)  # C at -4 vs R at -1
    assert np.isnan(out.loc[2, "esm1v_llr"]) and out.loc[2, plm.SCORED_FLAG] == 0
    assert out.loc[1, plm.SCORED_FLAG] == 1


def test_variant_llrs_accepts_kernel_output_columns():
    models = [*plm.ESM1V_MODELS, plm.ESM2_MODEL]
    lp = pd.DataFrame(logprob_rows("ENST1", 3, models, [-1] * 6, [-2] * 6))
    lp = lp.assign(ref_aa="R", model_revision="abc", window_start=0, window_end=3)
    out = plm.variant_llrs(mapping_frame(), lp).set_index("VariationID")
    assert out.loc[1, "esm1v_llr"] == pytest.approx(-1.0)


def test_variant_llrs_requires_all_esm1v_members():
    models = list(plm.ESM1V_MODELS[:4])
    lp = pd.DataFrame(logprob_rows("ENST1", 3, models, [-1] * 4, [-2] * 4))
    out = plm.variant_llrs(mapping_frame(), lp).set_index("VariationID")
    assert np.isnan(out.loc[1, "esm1v_llr"]) and out.loc[1, plm.SCORED_FLAG] == 0


def test_variant_llrs_rejects_duplicates_and_nonfinite():
    models = [plm.ESM2_MODEL, plm.ESM2_MODEL]
    with pytest.raises(ValueError):
        plm.variant_llrs(mapping_frame(), pd.DataFrame(logprob_rows("ENST1", 3, models, [-1, -1], [-2, -2])))
    bad = pd.DataFrame(logprob_rows("ENST1", 3, [plm.ESM2_MODEL], [float("inf")], [-2]))
    with pytest.raises(ValueError):
        plm.variant_llrs(mapping_frame(), bad)


def test_attach_plm_features_left_join():
    matrix = pd.DataFrame({"VariationID": [1, 2, 4], "y": [1, 0, 0]})
    feats = pd.DataFrame({"VariationID": [1, 2], "esm1v_llr": [-5.0, np.nan], "esm2_llr": [-4.0, np.nan],
                          "plm_scored": [1, 0]})
    out = plm.attach_plm_features(matrix, feats)
    assert out["plm_scored"].tolist() == [1, 0, 0]
    assert out["esm1v_llr"].isna().tolist() == [False, True, True]
    with pytest.raises(ValueError):
        plm.attach_plm_features(out, feats)


def test_feature_sets_extend_original18_without_changing_it():
    assert plm.PLM_FEATURE_SETS["ORIGINAL18"] == ALL_FEATURES == FEATURE_SETS["L6_full"]
    for name, cols in plm.PLM_FEATURE_SETS.items():
        assert len(cols) == len(set(cols)), name
    assert plm.PLM_FEATURE_SETS["ORIGINAL18_ESM1V"][:18] == ALL_FEATURES
    assert set(plm.SCORE_COLUMNS) <= plm.PLM_CONTINUOUS and plm.SCORED_FLAG not in plm.PLM_CONTINUOUS
    assert set(plm.MODEL_REVISIONS) == {*plm.ESM1V_MODELS, plm.ESM2_MODEL}


def test_impute_uses_inner_train_median_for_plm_scores():
    train = pd.DataFrame({"esm1v_llr": [-1.0, -3.0, np.nan, -9.0], "plm_scored": [1, 1, 0, 1]})
    test = pd.DataFrame({"esm1v_llr": [np.nan], "plm_scored": [0]})
    tr, te = _impute(train, test, continuous=plm.PLM_CONTINUOUS)
    assert te["esm1v_llr"].item() == pytest.approx(-3.0) and tr["esm1v_llr"].iloc[2] == pytest.approx(-3.0)
    _, te_default = _impute(train, test)  # frozen rule: not a continuous column -> 0
    assert te_default["esm1v_llr"].item() == 0


def synthetic_matrix(n_groups=30, per_group=8, seed=3):
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_groups):
        for i in range(per_group):
            y = int(rng.random() < 0.4)
            row = {"VariationID": g * 100 + i, "y": y, "p4_gene_group": f"GG{g:03d}", "GeneSymbol": f"G{g}",
                   "consequence_class": "coding_nontruncating", "p4_review_star": 1, "cohort": 1}
            row |= {c: float(rng.normal(y, 1.0)) for c in ALL_FEATURES}
            row["esm1v_llr"] = float(rng.normal(-3 * y, 1.0)) if i % 2 else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def test_run_experiment_default_feature_map_is_unchanged():
    m = synthetic_matrix()
    pf_a, oof_a = run_experiment(m, "cohort", ["L6_full"], seeds=[101])
    pf_b, oof_b = run_experiment(m, "cohort", ["L6_full"], seeds=[101], feature_map=FEATURE_SETS)
    pd.testing.assert_frame_equal(pf_a, pf_b)
    pd.testing.assert_frame_equal(oof_a["L6_full"], oof_b["L6_full"])
    pf_c, oof_c = run_experiment(m, "cohort", ["ORIGINAL18"], seeds=[101],
                                 feature_map=plm.PLM_FEATURE_SETS, continuous=plm.PLM_CONTINUOUS)
    np.testing.assert_array_equal(oof_a["L6_full"]["prob"], oof_c["ORIGINAL18"]["prob"])
    pf_d, _ = run_experiment(m, "cohort", ["ESM1V_ONLY"], seeds=[101],
                             feature_map=plm.PLM_FEATURE_SETS, continuous=plm.PLM_CONTINUOUS)
    assert pf_d["n_features"].eq(1).all()


# --- paired gene-group bootstrap ---------------------------------------------
def test_paired_bootstrap_matches_explicit_resampling_oracle():
    y = np.array([0, 1, 1, 0, 1, 0, 0, 1])
    a = np.array([0.2, 0.2, 0.8, 0.4, 0.9, 0.6, 0.1, 0.3])
    b = np.array([0.1, 0.7, 0.8, 0.2, 0.9, 0.3, 0.1, 0.6])
    g = np.array(["a", "a", "b", "c", "c", "c", "c", "d"])
    out = paired_group_bootstrap_auc(y, a, b, groups=g, n_boot=400, seed=11)
    rng = np.random.default_rng(11)
    members = [np.flatnonzero(g == k) for k in np.unique(g)]
    deltas = []
    for _ in range(400):
        ix = np.concatenate([members[k] for k in rng.integers(0, 4, 4)])
        if np.unique(y[ix]).size == 2:
            deltas.append(roc_auc_score(y[ix], b[ix]) - roc_auc_score(y[ix], a[ix]))
    d = out["delta_b_minus_a"]
    assert d["estimate"] == pytest.approx(roc_auc_score(y, b) - roc_auc_score(y, a))
    assert d["n_valid"] == len(deltas)
    assert d["bootstrap_mean"] == pytest.approx(np.mean(deltas))
    assert [d["lo"], d["hi"]] == pytest.approx(np.percentile(deltas, [2.5, 97.5]))
    assert out["superior"] == (d["lo"] > 0)


def test_identical_scores_are_never_superior():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 200)
    s = rng.random(200)
    g = np.repeat([f"g{i}" for i in range(20)], 10)
    out = paired_group_bootstrap_auc(y, s, s, groups=g, n_boot=200)
    assert out["delta_b_minus_a"]["estimate"] == 0
    assert out["delta_b_minus_a"]["lo"] == 0 and out["delta_b_minus_a"]["hi"] == 0
    assert out["superior"] is False


def test_paired_bootstrap_detects_a_clear_improvement():
    rng = np.random.default_rng(1)
    y = rng.integers(0, 2, 600)
    g = np.repeat([f"g{i}" for i in range(60)], 10)
    weak = rng.random(600)
    strong = y + rng.normal(0, 0.3, 600)
    out = paired_group_bootstrap_auc(y, weak, strong, groups=g, n_boot=300)
    assert out["superior"] is True and out["delta_b_minus_a"]["lo"] > 0


@pytest.mark.parametrize("change", ["missing_group", "row_group", "one_group", "bad_label",
                                    "infinite", "misaligned", "one_class", "few_draws"])
def test_paired_bootstrap_rejects_invalid_input(change):
    y, a, b, g, n = [0, 1, 0, 1], [0.1, 0.9, 0.2, 0.8], [0.2, 0.8, 0.1, 0.9], ["a", "a", "b", "b"], 100
    if change == "missing_group": g = [None, "a", "b", "b"]
    if change == "row_group": g = ["row0", "row1", "row2", "row3"]
    if change == "one_group": g = ["a"] * 4
    if change == "bad_label": y = [0, 2, 0, 1]
    if change == "infinite": a = [float("inf"), 0.9, 0.2, 0.8]
    if change == "misaligned": b = [0.2, 0.8]
    if change == "one_class": y = [1, 1, 1, 1]
    if change == "few_draws": n = 50
    with pytest.raises(ValueError):
        paired_group_bootstrap_auc(y, a, b, groups=g, n_boot=n)
